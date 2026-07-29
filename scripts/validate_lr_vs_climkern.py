#!/usr/bin/env python3
"""Cross-validate the native lapse-rate kernel module (core/lr_kernel.py)
against ClimKern's own calc_T_feedbacks, fed the *same* base/perturbed
temperature response (docs/plan.md §5.1 "同一输入" control-variable
design).

Requires the `pycfram-kern` conda env (climkern + xesmf); this is the
*only* script in the M2 module that imports climkern/xesmf -- core/
stays dependency-free (docs/plan.md §3.3 "原则").

Usage (inside `pycfram-kern` env):
    python scripts/validate_lr_vs_climkern.py cesm2_4xco2_official
    python scripts/validate_lr_vs_climkern.py cesm2_4xco2_official --kernels CloudSat GFDL
"""
import os
import sys
import argparse
import numpy as np
import pandas as pd
import xarray as xr
from netCDF4 import Dataset

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.config import load_case

try:
    import climkern as ck
except ImportError as e:
    sys.exit(
        "ERROR: climkern (and its xesmf/esmpy dependency) is required for "
        "this script. Run from inside the `pycfram-kern` conda env "
        "(see docs/plan.md §3.2). Original error: %r" % e)

FILL = -999.0
MONTHS = pd.date_range('2000-01-15', periods=12, freq='MS') + pd.Timedelta(days=14)


def _mask_fill(arr):
    return np.where(np.abs(arr) > 900.0, np.nan, arr)


def tile_to_12_months(da):
    """Tile a (lat,lon) or (plev,lat,lon) DataArray (single state, no time
    dim) into 12 identical months with a datetime64 time coordinate
    (docs/plan.md §5.2 -- climkern's tile_data requires len(time)%12==0,
    and groupby(time.dt.month) requires datetime64, not integer months)."""
    tiled = xr.concat([da] * 12, dim='time')
    tiled = tiled.assign_coords(time=('time', MONTHS))
    return tiled


def load_state_for_climkern(pres_path, surf_path):
    """Read pyCFRAM's *_pres.nc / *_surf.nc into ClimKern-ready, 12-month
    tiled DataArrays (ta, ts, ps)."""
    ds_pres = xr.open_dataset(pres_path)
    ds_surf = xr.open_dataset(surf_path)

    ta = ds_pres['ta_lay'].isel(time=0, drop=True)  # (lev, lat, lon)
    ta = ta.rename({'lev': 'plev'})
    ta['plev'].attrs['units'] = 'hPa'
    ta.attrs['units'] = 'K'

    ts = ds_surf['ts'].isel(time=0, drop=True)
    ts.attrs['units'] = 'K'

    ps = ds_surf['ps'].isel(time=0, drop=True)
    ps.attrs['units'] = 'Pa'

    return tile_to_12_months(ta), tile_to_12_months(ts), tile_to_12_months(ps)


def area_weighted_mean(field, lat):
    w = np.cos(np.deg2rad(lat))
    w2d = np.broadcast_to(w[:, None], field.shape)
    mask = ~np.isnan(field)
    if not np.any(mask):
        return np.nan
    return np.sum(np.where(mask, field * w2d, 0.0)) / np.sum(np.where(mask, w2d, 0.0))


def area_weighted_corr(a, b, lat):
    w = np.cos(np.deg2rad(lat))
    w2d = np.broadcast_to(w[:, None], a.shape)
    mask = ~np.isnan(a) & ~np.isnan(b)
    aw, bw, ww = a[mask], b[mask], w2d[mask]
    ww = ww / ww.sum()
    ma = np.sum(aw * ww)
    mb = np.sum(bw * ww)
    cov = np.sum(ww * (aw - ma) * (bw - mb))
    va = np.sum(ww * (aw - ma) ** 2)
    vb = np.sum(ww * (bw - mb) ** 2)
    return cov / np.sqrt(va * vb)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('case')
    parser.add_argument('--kernels', nargs='+', default=None)
    parser.add_argument('--sky', default=None, choices=['all-sky', 'clear-sky'])
    args = parser.parse_args()

    cfg = load_case(args.case)
    lr_cfg = cfg.get('lapse_rate', {})
    kernel_names = args.kernels or lr_cfg.get('kernels', ['CloudSat', 'GFDL'])
    sky = args.sky or lr_cfg.get('sky', 'all-sky')

    lr_kernel_nc = os.path.join(cfg['_output_dir'], 'lr_kernel.nc')
    if not os.path.exists(lr_kernel_nc):
        sys.exit("ERROR: %s not found. Run scripts/compute_lr_kernel.py %s first."
                  % (lr_kernel_nc, args.case))

    print("=" * 60)
    print("validate_lr_vs_climkern: case=%s kernels=%s sky=%s"
          % (args.case, kernel_names, sky))
    print("=" * 60)

    with Dataset(lr_kernel_nc) as d:
        lat = np.array(d.variables['lat'][:], dtype=np.float64)
        lon = np.array(d.variables['lon'][:], dtype=np.float64)
        native_lr = {k: _mask_fill(np.array(d.variables['dR_lr_%s' % k][:]))
                     for k in kernel_names}

    print("Loading base/perturbed state for ClimKern reference ...")
    ctrl_ta, ctrl_ts, ctrl_ps = load_state_for_climkern(
        cfg['input']['base_pres'], cfg['input']['base_surf'])
    pert_ta, pert_ts, pert_ps = load_state_for_climkern(
        cfg['input']['perturbed_pres'], cfg['input']['perturbed_surf'])

    report_lines = []
    report_lines.append("Lapse-Rate kernel cross-validation report")
    report_lines.append("case=%s sky=%s" % (args.case, sky))
    report_lines.append("")

    ref_path = os.path.join(cfg['_output_dir'], 'lr_climkern_ref.nc')
    ref_nc = Dataset(ref_path, 'w')
    ref_nc.createDimension('lat', lat.size)
    ref_nc.createDimension('lon', lon.size)
    ref_nc.createVariable('lat', 'f8', ('lat',))[:] = lat
    ref_nc.createVariable('lon', 'f8', ('lon',))[:] = lon
    ref_nc.setncattr('source', 'climkern.calc_T_feedbacks (xesmf reference)')

    all_pass = True
    for kname in kernel_names:
        print("\n--- kernel: %s ---" % kname)
        lr_ck, pl_ck = ck.calc_T_feedbacks(
            ctrl_ta, ctrl_ts, ctrl_ps, pert_ta, pert_ts, pert_ps,
            pert_trop=None, kern=kname, sky=sky, fixRH=False)
        lr_ck_annual = lr_ck.mean(dim='time').values  # (lat,lon)

        native = native_lr[kname]
        corr = area_weighted_corr(native, lr_ck_annual, lat)
        mean_native = area_weighted_mean(native, lat)
        mean_ck = area_weighted_mean(lr_ck_annual, lat)
        rel_diff = abs(mean_native - mean_ck) / abs(mean_ck) if mean_ck != 0 else np.nan

        gate = (corr >= 0.85) and (rel_diff <= 0.15)
        all_pass = all_pass and gate

        line = ("%-10s  corr=%.4f  mean_native=%.4f W/m2  mean_ck=%.4f W/m2  "
                "rel_diff=%.2f%%  -> %s"
                % (kname, corr, mean_native, mean_ck, rel_diff * 100,
                   "PASS" if gate else "FAIL"))
        print(line)
        report_lines.append(line)

        v = ref_nc.createVariable('dR_lr_%s' % kname, 'f8', ('lat', 'lon'), fill_value=FILL)
        v[:] = np.where(np.isnan(lr_ck_annual), FILL, lr_ck_annual)
        v.units = 'W m-2'
        v.long_name = '%s ClimKern reference lapse-rate TOA LW radiative perturbation' % kname

    ref_nc.close()
    print("Wrote %s" % ref_path)

    report_lines.append("")
    report_lines.append("OVERALL: %s" % ("PASS" if all_pass else "FAIL"))

    report_path = os.path.join(cfg['_output_dir'], 'lr_validation_report.txt')
    with open(report_path, 'w') as f:
        f.write('\n'.join(report_lines) + '\n')
    print("\nWrote %s" % report_path)

    if not all_pass:
        sys.exit(1)


if __name__ == '__main__':
    main()
