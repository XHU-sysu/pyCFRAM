#!/usr/bin/env python3
"""M3 route i: compute per-physical-process Lapse-Rate attribution.

Reads cases/<case>/output/cfram_result.nc (each dT_<term>) and
cases/<case>/output/lr_kernel.nc (the M2 total ΔR_LR from dT_observed,
used as the additivity-residual reference), and writes
cases/<case>/output/lr_attribution.nc with dR_lr_from_<term>_<kernel> for
every configured kernel + term.

Usage:
    python scripts/compute_lr_attribution.py cesm2_4xco2_official
    python scripts/compute_lr_attribution.py cesm2_4xco2_official --terms q atmdyn sfcdyn
"""
import os
import sys
import argparse
import numpy as np
from netCDF4 import Dataset

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.config import load_case
from core.kernels import KernelSet
from core.lr_attribution import attribute_lr, additivity_residual, DEFAULT_TERMS
from data.kernel_source import get_kernel_path

FILL = -999.0


def _mask_fill(arr):
    return np.where(np.abs(arr) > 900.0, np.nan, arr)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('case')
    parser.add_argument('--terms', nargs='+', default=None,
                         help='Override default term list: %s' % DEFAULT_TERMS)
    parser.add_argument('--kernels', nargs='+', default=None)
    parser.add_argument('--month', default=None)
    parser.add_argument('--sky', default=None, choices=['all-sky', 'clear-sky'])
    args = parser.parse_args()

    cfg = load_case(args.case)
    lr_cfg = cfg.get('lapse_rate', {})
    kernel_names = args.kernels or lr_cfg.get('kernels', ['CloudSat', 'GFDL'])
    month = args.month or lr_cfg.get('kernel_months', 'annual')
    if isinstance(month, list):
        month = month[0] if len(month) == 1 else 'annual'
    sky = args.sky or lr_cfg.get('sky', 'all-sky')
    terms = args.terms or DEFAULT_TERMS

    result_nc = os.path.join(cfg['_output_dir'], 'cfram_result.nc')
    lr_kernel_nc = os.path.join(cfg['_output_dir'], 'lr_kernel.nc')
    if not os.path.exists(lr_kernel_nc):
        sys.exit("ERROR: %s missing -- run scripts/compute_lr_kernel.py %s first "
                  "(M3 route i needs the M2 total for additivity comparison)."
                  % (lr_kernel_nc, args.case))

    print("=" * 60)
    print("compute_lr_attribution: case=%s terms=%s kernels=%s"
          % (args.case, terms, kernel_names))
    print("=" * 60)

    with Dataset(result_nc) as d:
        lev = np.array(d.variables['lev'][:], dtype=np.float64)
        lat = np.array(d.variables['lat'][:], dtype=np.float64)
        lon = np.array(d.variables['lon'][:], dtype=np.float64)
        skin_idx = np.argmax(lev > 1005.0) if np.any(lev > 1005.0) else lev.size - 1
        atm_mask = np.arange(lev.size) != skin_idx
        plev_native_pa = lev[atm_mask] * 100.0

        dT_terms = {}
        available = []
        for term in terms:
            vname = 'dT_%s' % term
            if vname not in d.variables:
                print("  (skip %s: %s not in cfram_result.nc)" % (term, vname))
                continue
            arr = _mask_fill(np.array(d.variables[vname][:]))
            dT_terms[term] = (arr[atm_mask], arr[skin_idx])
            available.append(term)

    surf_path = cfg['input']['perturbed_surf']
    with Dataset(surf_path) as d:
        ps = np.array(d.variables['ps'][:], dtype=np.float64)
        if ps.ndim == 3:
            ps = ps[0]

    with Dataset(lr_kernel_nc) as d:
        total_by_kernel = {k: _mask_fill(np.array(d.variables['dR_lr_%s' % k][:]))
                            for k in kernel_names}

    out_path = os.path.join(cfg['_output_dir'], 'lr_attribution.nc')
    nc_out = Dataset(out_path, 'w')
    nc_out.createDimension('lat', lat.size)
    nc_out.createDimension('lon', lon.size)
    nc_out.createVariable('lat', 'f8', ('lat',))[:] = lat
    nc_out.createVariable('lon', 'f8', ('lon',))[:] = lon
    nc_out.setncattr('terms', ','.join(available))
    nc_out.setncattr('month', str(month))
    nc_out.setncattr('sky', sky)

    for kname in kernel_names:
        kpath = get_kernel_path(kname)
        ks = KernelSet(kpath, name=kname).regrid_to(lat, lon)

        by_term = attribute_lr(dT_terms, ps, plev_native_pa, lat, ks,
                                month=month, sky=sky)

        dR_lr_by_term = {}
        for term, result in by_term.items():
            dR_lr_by_term[term] = result['dR_lr']
            v = nc_out.createVariable('dR_lr_from_%s_%s' % (term, kname), 'f8',
                                       ('lat', 'lon'), fill_value=FILL)
            v[:] = np.where(np.isnan(result['dR_lr']), FILL, result['dR_lr'])
            v.units = 'W m-2'
            v.long_name = ("%s-attributed lapse-rate TOA LW perturbation (%s kernel)"
                            % (term, kname))

        resid = additivity_residual(total_by_kernel[kname], dR_lr_by_term)
        resid_mean = np.nanmean(np.abs(resid))
        total_mean = np.nanmean(np.abs(total_by_kernel[kname]))
        pct = 100.0 * resid_mean / total_mean if total_mean else np.nan
        print("  %-10s additivity |residual| mean = %.4f W/m^2 (%.1f%% of |total| mean)"
              % (kname, resid_mean, pct))

        v = nc_out.createVariable('dR_lr_additivity_residual_%s' % kname, 'f8',
                                   ('lat', 'lon'), fill_value=FILL)
        v[:] = np.where(np.isnan(resid), FILL, resid)
        v.units = 'W m-2'
        v.long_name = ("Total ΔR_LR (from dT_observed) minus Σ per-process "
                        "ΔR_LR (%s kernel) -- CFRAM first-order non-additivity, "
                        "not expected to be zero (docs/plan.md R7)" % kname)

    nc_out.close()
    print("Wrote %s" % out_path)


if __name__ == '__main__':
    main()
