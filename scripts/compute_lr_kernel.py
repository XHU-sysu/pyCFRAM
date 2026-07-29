#!/usr/bin/env python3
"""Compute the radiative-kernel Lapse-Rate / Planck decomposition (native,
xesmf-free implementation) for a pyCFRAM case.

Reads `cases/<case>/output/cfram_result.nc::dT_observed` (the model's
total temperature response -- same quantity ClimKern's calc_T_feedbacks
would be fed, see docs/plan.md §5.1) plus `perturbed_surf.nc::ps`, applies
each configured kernel (default: CloudSat=Kramer, GFDL), and writes
`cases/<case>/output/lr_kernel.nc`.

Usage:
    python scripts/compute_lr_kernel.py cesm2_4xco2_official
    python scripts/compute_lr_kernel.py cesm2_4xco2_official --kernels CloudSat GFDL
"""
import os
import sys
import argparse
import numpy as np
from netCDF4 import Dataset

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.config import load_case, PROJECT_ROOT
from core.kernels import KernelSet
from core.lr_kernel import extract_and_apply
from data.kernel_source import get_kernel_path

FILL = -999.0


def _mask_fill(arr):
    return np.where(np.abs(arr) > 900.0, np.nan, arr)


def check_freshness(nc_path):
    """WP-M2.R / §2.1 新鲜度门: ocndyn+lhflx+shflx == sfcdyn to machine
    precision. Raises if the file is a pre-fedcf4d stale output."""
    with Dataset(nc_path) as d:
        ocndyn = _mask_fill(d.variables['dT_ocndyn'][:])
        lhflx = _mask_fill(d.variables['dT_lhflx'][:])
        shflx = _mask_fill(d.variables['dT_shflx'][:])
        sfcdyn = _mask_fill(d.variables['dT_sfcdyn'][:])
    resid = np.nanmax(np.abs(ocndyn + lhflx + shflx - sfcdyn))
    print("  freshness gate: identity residual max = %.3e K" % resid)
    if resid >= 1e-10:
        raise RuntimeError(
            "STALE cfram_result.nc (identity residual %.3e K >= 1e-10). "
            "Rerun the case per docs/plan.md WP-M2.R before using it for "
            "the lapse-rate kernel module." % resid)


def load_case_fields(cfg):
    out_dir = cfg['_output_dir']
    result_nc = os.path.join(out_dir, 'cfram_result.nc')
    check_freshness(result_nc)

    with Dataset(result_nc) as d:
        lev = np.array(d.variables['lev'][:], dtype=np.float64)  # hPa
        lat = np.array(d.variables['lat'][:], dtype=np.float64)
        lon = np.array(d.variables['lon'][:], dtype=np.float64)
        dT_observed = _mask_fill(np.array(d.variables['dT_observed'][:]))

    skin_idx = np.argmax(lev > 1005.0) if np.any(lev > 1005.0) else lev.size - 1
    atm_mask = np.arange(lev.size) != skin_idx
    plev_native_pa = lev[atm_mask] * 100.0
    dT_atm = dT_observed[atm_mask]
    dT_skin = dT_observed[skin_idx]

    surf_path = cfg['input']['perturbed_surf']
    with Dataset(surf_path) as d:
        ps = np.array(d.variables['ps'][:], dtype=np.float64)
        if ps.ndim == 3:
            ps = ps[0]

    return dict(lat=lat, lon=lon, plev_native_pa=plev_native_pa,
                dT_atm=dT_atm, dT_skin=dT_skin, ps=ps)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('case')
    parser.add_argument('--kernels', nargs='+', default=None,
                         help='Override case.yaml lapse_rate.kernels')
    parser.add_argument('--month', default=None,
                         help="Override case.yaml lapse_rate.kernel_months "
                              "('annual' or an int 1-12)")
    parser.add_argument('--sky', default=None, choices=['all-sky', 'clear-sky'])
    args = parser.parse_args()

    cfg = load_case(args.case)
    lr_cfg = cfg.get('lapse_rate', {})
    kernel_names = args.kernels or lr_cfg.get('kernels', ['CloudSat', 'GFDL'])
    month_cfg = lr_cfg.get('kernel_months', 'annual')
    month = args.month if args.month is not None else month_cfg
    if isinstance(month, list):
        month = month[0] if len(month) == 1 else 'annual'
    sky = args.sky or lr_cfg.get('sky', 'all-sky')

    print("=" * 60)
    print("compute_lr_kernel: case=%s kernels=%s month=%s sky=%s"
          % (args.case, kernel_names, month, sky))
    print("=" * 60)

    fields = load_case_fields(cfg)
    lat, lon = fields['lat'], fields['lon']

    out_path = os.path.join(cfg['_output_dir'], 'lr_kernel.nc')
    os.makedirs(cfg['_output_dir'], exist_ok=True)
    nc_out = Dataset(out_path, 'w')
    nc_out.createDimension('lat', lat.size)
    nc_out.createDimension('lon', lon.size)
    v = nc_out.createVariable('lat', 'f8', ('lat',)); v[:] = lat
    v = nc_out.createVariable('lon', 'f8', ('lon',)); v[:] = lon
    nc_out.setncattr('month', str(month))
    nc_out.setncattr('sky', sky)
    nc_out.setncattr('source', 'core.lr_kernel.extract_and_apply (native, xesmf-free)')

    for kname in kernel_names:
        kpath = get_kernel_path(kname)
        ks_native = KernelSet(kpath, name=kname)
        ks = ks_native.regrid_to(lat, lon)

        result = extract_and_apply(
            fields['dT_atm'], fields['dT_skin'], fields['ps'],
            fields['plev_native_pa'], lat, ks, month=month, sky=sky)

        dR_lr, dR_pl = result['dR_lr'], result['dR_pl']
        print("  %-10s dR_lr: mean=%.4f W/m^2  dR_pl: mean=%.4f W/m^2"
              % (kname, np.nanmean(dR_lr), np.nanmean(dR_pl)))

        v = nc_out.createVariable('dR_lr_%s' % kname, 'f8', ('lat', 'lon'),
                                   fill_value=FILL)
        v[:] = np.where(np.isnan(dR_lr), FILL, dR_lr)
        v.units = 'W m-2'
        v.long_name = '%s kernel lapse-rate TOA LW radiative perturbation' % kname

        v = nc_out.createVariable('dR_pl_%s' % kname, 'f8', ('lat', 'lon'),
                                   fill_value=FILL)
        v[:] = np.where(np.isnan(dR_pl), FILL, dR_pl)
        v.units = 'W m-2'
        v.long_name = '%s kernel Planck TOA LW radiative perturbation' % kname

    nc_out.close()
    print("Wrote %s" % out_path)


if __name__ == '__main__':
    main()
