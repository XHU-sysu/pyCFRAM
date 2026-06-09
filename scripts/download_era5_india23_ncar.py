#!/usr/bin/env python3
"""Build pyCFRAM-ready ERA5 files for india_wb23 from NCAR ERA5 on AWS.

This avoids CDS request queues by reading NCAR's public NetCDF4 ERA5 files
directly from S3 and writing only the India-Bangladesh regional subset.

Output layout matches data/era5_source.py:
  OUTDIR/
    era5_pl_{t,q,o3,cc,clwc,ciwc}_{YYYY04}.nc
    era5_sl_{YYYY04}/
      data_stream-oper_stepType-instant.nc   # skt, sp
      data_stream-oper_stepType-accum.nc     # ssrd, ssr, tisr
      data_stream-oper_stepType-accum_flux.nc # slhf, sshf
      data_stream-oper_stepType-max.nc       # mx2t
"""

import argparse
import calendar
import os
from datetime import datetime, timedelta

import numpy as np
import pandas as pd
import xarray as xr

try:
    import s3fs
except ImportError:
    s3fs = None


BUCKET = "nsf-ncar-era5"
OUTDIR = "/disk/r074/lzhenn/workspace/ust-jumper/pyCFRAM/era5_data/india_wb23_ncar"

MONTH = 4
CLIM_YEARS = list(range(1991, 2021))
EVENT_YEAR = 2023

LAT_NORTH = 35.0
LAT_SOUTH = 10.0
LON_WEST = 65.0
LON_EAST = 110.0

PL_VARS = {
    "t": ("128_130_t", "T"),
    "q": ("128_133_q", "Q"),
    "o3": ("128_203_o3", "O3"),
    "cc": ("128_248_cc", "CC"),
    "clwc": ("128_246_clwc", "CLWC"),
    "ciwc": ("128_247_ciwc", "CIWC"),
}

DIAG_PL_VARS = {
    "u": ("128_131_u", "U"),
    "v": ("128_132_v", "V"),
}

INSTANT_VARS = {
    "skt": ("128_235_skt", "SKT"),
    "sp": ("128_134_sp", "SP"),
}

DIAG_INSTANT_VARS = {
    "t2m": ("128_167_2t", "VAR_2T"),
    "d2m": ("128_168_2d", "VAR_2D"),
}

ACCUM_VARS = {
    "ssrd": ("128_169_ssrd", "SSRD"),
    "ssr": ("128_176_ssr", "SSR"),
    "tisr": ("128_212_tisr", "TISR"),
    "slhf": ("128_147_slhf", "SLHF"),
    "sshf": ("128_146_sshf", "SSHF"),
}

MAX_VARS = {
    "mx2t": ("128_201_mx2t", "MX2T"),
}

ENCODING = {
    "zlib": True,
    "complevel": 1,
    "shuffle": True,
    "dtype": "float32",
}


def s3_path(*parts):
    return "/".join([BUCKET, *parts])


def yyyymm(year, month=None):
    if month is None:
        month = MONTH
    return f"{year}{month:02d}"


def _wanted_times(year):
    last_day = calendar.monthrange(year, MONTH)[1]
    return pd.date_range(
        f"{year}-{MONTH:02d}-01 00:00",
        f"{year}-{MONTH:02d}-{last_day:02d} 18:00",
        freq="6h",
    )


def _open_dataset(fs, path):
    handle = fs.open(path, "rb")
    try:
        ds = xr.open_dataset(handle, engine="h5netcdf", chunks=None)
        return handle, ds
    except Exception:
        handle.close()
        raise


def _close(handle, ds):
    ds.close()
    handle.close()


def _sel_region(da):
    return da.sel(latitude=slice(LAT_NORTH, LAT_SOUTH),
                  longitude=slice(LON_WEST, LON_EAST))


def _write_dataset(path, data_vars, coords, attrs=None, overwrite=False):
    if os.path.exists(path) and not overwrite:
        size = os.path.getsize(path) / 1e6
        print(f"  SKIP  {os.path.basename(path)} ({size:.1f} MB)")
        return
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = f"{path}.tmp"
    if os.path.exists(tmp):
        os.remove(tmp)
    ds = xr.Dataset(data_vars=data_vars, coords=coords, attrs=attrs or {})
    enc = {name: ENCODING for name in data_vars}
    ds.to_netcdf(tmp, engine="h5netcdf", encoding=enc)
    ds.close()
    os.replace(tmp, path)
    size = os.path.getsize(path) / 1e6
    print(f"  wrote {path} ({size:.1f} MB)", flush=True)


def download_pl_var(fs, year, var_short, overwrite=False, dry_run=False,
                    var_table=None):
    ym = yyyymm(year)
    var_table = var_table or PL_VARS
    code, ncar_name = var_table[var_short]
    grid = "ll025uv" if var_short in ("u", "v") else "ll025sc"
    out = os.path.join(OUTDIR, f"era5_pl_{var_short}_{ym}.nc")
    if os.path.exists(out) and not overwrite:
        print(f"  SKIP  PL {var_short} {ym}")
        return
    print(f"  PL {var_short} {ym}", flush=True)
    if dry_run:
        return

    parts = []
    last_day = calendar.monthrange(year, MONTH)[1]
    for day in range(1, last_day + 1):
        d0 = f"{ym}{day:02d}00"
        d1 = f"{ym}{day:02d}23"
        path = s3_path("e5.oper.an.pl", ym,
                       f"e5.oper.an.pl.{code}.{grid}.{d0}_{d1}.nc")
        handle, ds = _open_dataset(fs, path)
        try:
            da = _sel_region(ds[ncar_name].isel(time=[0, 6, 12, 18])).load()
            da = da.rename({"level": "pressure_level"})
            parts.append(da.astype("float32"))
        finally:
            _close(handle, ds)

    all_da = xr.concat(parts, dim="time")
    all_da.name = var_short
    _write_dataset(
        out,
        {var_short: all_da},
        {
            "time": all_da.time,
            "pressure_level": all_da.pressure_level,
            "latitude": all_da.latitude,
            "longitude": all_da.longitude,
        },
        attrs={"source": "NCAR ERA5 AWS", "created": datetime.utcnow().isoformat()},
        overwrite=True,
    )


def download_instant_sl(fs, year, overwrite=False, dry_run=False,
                        var_table=None, out_name="data_stream-oper_stepType-instant.nc"):
    ym = yyyymm(year)
    sl_dir = os.path.join(OUTDIR, f"era5_sl_{ym}")
    out = os.path.join(sl_dir, out_name)
    if os.path.exists(out) and not overwrite:
        print(f"  SKIP  SL instant {ym} {out_name}")
        return
    print(f"  SL instant {ym} {out_name}", flush=True)
    if dry_run:
        return

    var_table = var_table or INSTANT_VARS
    data_vars = {}
    coords = None
    last_day = calendar.monthrange(year, MONTH)[1]
    for out_var, (code, ncar_name) in var_table.items():
        path = s3_path("e5.oper.an.sfc", ym,
                       f"e5.oper.an.sfc.{code}.ll025sc."
                       f"{ym}0100_{ym}{last_day:02d}23.nc")
        handle, ds = _open_dataset(fs, path)
        try:
            if ncar_name not in ds.variables:
                candidates = [v for v in ds.variables
                              if v not in ("time", "latitude", "longitude")]
                if len(candidates) == 1:
                    ncar_name = candidates[0]
                else:
                    raise KeyError(
                        f"{path}: cannot find {ncar_name}; candidates={candidates}")
            da = _sel_region(ds[ncar_name].isel(time=slice(0, None, 6))).load()
            da = da.astype("float32")
            data_vars[out_var] = da
            if coords is None:
                coords = {
                    "time": da.time,
                    "latitude": da.latitude,
                    "longitude": da.longitude,
                }
        finally:
            _close(handle, ds)

    _write_dataset(out, data_vars, coords,
                   attrs={"source": "NCAR ERA5 AWS"}, overwrite=True)


def _boundary_month_paths(stream, year, code):
    ym = yyyymm(year)
    first_day = datetime(year, MONTH, 1)
    prev_day = first_day - timedelta(days=1)
    next_month = (first_day.replace(day=28) + timedelta(days=4)).replace(day=1)
    prev = yyyymm(prev_day.year, prev_day.month)
    nxt = yyyymm(next_month.year, next_month.month)
    return [
        s3_path(stream, prev,
                f"{stream}.{code}.ll025sc.{prev}1606_{ym}0106.nc"),
        s3_path(stream, ym,
                f"{stream}.{code}.ll025sc.{ym}0106_{ym}1606.nc"),
        s3_path(stream, ym,
                f"{stream}.{code}.ll025sc.{ym}1606_{nxt}0106.nc"),
    ]


def _load_forecast_series(fs, paths, ncar_name, year, accum_to_hourly=False):
    wanted = _wanted_times(year)
    wanted_set = set(wanted)
    pieces = []
    times = []

    for path in paths:
        handle, ds = _open_dataset(fs, path)
        try:
            da = _sel_region(ds[ncar_name])
            init_times = pd.to_datetime(ds.forecast_initial_time.values)
            hours = ds.forecast_hour.values.astype(int)
            hour_to_index = {int(h): i for i, h in enumerate(hours)}

            for i, init_time in enumerate(init_times):
                for hour in (6, 12):
                    valid_time = init_time + pd.Timedelta(hours=hour)
                    if valid_time not in wanted_set:
                        continue
                    j = hour_to_index[hour]
                    value = da.isel(forecast_initial_time=i, forecast_hour=j)
                    if accum_to_hourly:
                        if hour == 12:
                            j6 = hour_to_index[6]
                            value = value - da.isel(
                                forecast_initial_time=i, forecast_hour=j6)
                        value = value / 6.0
                    value = value.load().astype("float32")
                    value = xr.DataArray(
                        np.asarray(value.data, dtype=np.float32),
                        dims=("latitude", "longitude"),
                        coords={
                            "latitude": value.latitude,
                            "longitude": value.longitude,
                        },
                        name=ncar_name,
                    )
                    pieces.append(value)
                    times.append(valid_time.to_datetime64())
        finally:
            _close(handle, ds)

    if len(times) != len(wanted):
        missing = sorted(set(wanted.values) - set(times))
        raise RuntimeError(f"{ncar_name}: expected {len(wanted)} times, "
                           f"got {len(times)}, missing={missing[:4]}")

    order = np.argsort(np.array(times))
    sorted_pieces = [pieces[i] for i in order]
    sorted_times = np.array(times)[order]
    out = xr.concat(sorted_pieces, dim="time")
    out = out.assign_coords(time=sorted_times)
    return out


def download_accum_sl(fs, year, overwrite=False, dry_run=False):
    ym = yyyymm(year)
    sl_dir = os.path.join(OUTDIR, f"era5_sl_{ym}")
    out_main = os.path.join(sl_dir, "data_stream-oper_stepType-accum.nc")
    out_flux = os.path.join(sl_dir, "data_stream-oper_stepType-accum_flux.nc")
    if (os.path.exists(out_main) and os.path.exists(out_flux)
            and not overwrite):
        print(f"  SKIP  SL accum {ym}")
        return
    print(f"  SL accum/flux {ym}", flush=True)
    if dry_run:
        return

    all_vars = {}
    coords = None
    for out_name, (code, ncar_name) in ACCUM_VARS.items():
        paths = _boundary_month_paths("e5.oper.fc.sfc.accumu", year, code)
        da = _load_forecast_series(
            fs, paths, ncar_name, year, accum_to_hourly=True)
        da.name = out_name
        all_vars[out_name] = da
        if coords is None:
            coords = {
                "time": da.time,
                "latitude": da.latitude,
                "longitude": da.longitude,
            }

    _write_dataset(
        out_main,
        {k: all_vars[k] for k in ("ssrd", "ssr", "tisr")},
        coords,
        attrs={"source": "NCAR ERA5 AWS; 6h accum converted to hourly equivalent"},
        overwrite=True,
    )
    _write_dataset(
        out_flux,
        {k: all_vars[k] for k in ("slhf", "sshf")},
        coords,
        attrs={"source": "NCAR ERA5 AWS; 6h accum converted to hourly equivalent"},
        overwrite=True,
    )


def download_max_sl(fs, year, overwrite=False, dry_run=False):
    ym = yyyymm(year)
    sl_dir = os.path.join(OUTDIR, f"era5_sl_{ym}")
    out = os.path.join(sl_dir, "data_stream-oper_stepType-max.nc")
    if os.path.exists(out) and not overwrite:
        print(f"  SKIP  SL max {ym}")
        return
    print(f"  SL max {ym}", flush=True)
    if dry_run:
        return

    data_vars = {}
    coords = None
    for out_name, (code, ncar_name) in MAX_VARS.items():
        paths = _boundary_month_paths("e5.oper.fc.sfc.minmax", year, code)
        da = _load_forecast_series(
            fs, paths, ncar_name, year, accum_to_hourly=False)
        da.name = out_name
        data_vars[out_name] = da
        coords = {
            "time": da.time,
            "latitude": da.latitude,
            "longitude": da.longitude,
        }

    _write_dataset(out, data_vars, coords,
                   attrs={"source": "NCAR ERA5 AWS"}, overwrite=True)


def download_year(fs, year, overwrite=False, dry_run=False):
    print(f"\n=== {year} ===", flush=True)
    for var_short in PL_VARS:
        download_pl_var(fs, year, var_short, overwrite, dry_run)
    download_instant_sl(fs, year, overwrite, dry_run)
    download_accum_sl(fs, year, overwrite, dry_run)
    download_max_sl(fs, year, overwrite, dry_run)


def download_diag_year(fs, year, overwrite=False, dry_run=False):
    print(f"\n=== diagnostics {year} ===", flush=True)
    for var_short in DIAG_PL_VARS:
        download_pl_var(fs, year, var_short, overwrite, dry_run,
                        var_table=DIAG_PL_VARS)
    download_instant_sl(
        fs, year, overwrite, dry_run,
        var_table=DIAG_INSTANT_VARS,
        out_name="data_stream-oper_stepType-instant_diag.nc")


def main():
    global OUTDIR

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--start-year", type=int, default=min(CLIM_YEARS))
    ap.add_argument("--end-year", type=int, default=EVENT_YEAR)
    ap.add_argument("--outdir", default=OUTDIR)
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--mode", choices=("cfram", "diag", "all"), default="cfram",
                    help="cfram: original pyCFRAM inputs; diag: u/v + t2m/d2m; all: both")
    args = ap.parse_args()

    OUTDIR = args.outdir

    years = [y for y in CLIM_YEARS + [EVENT_YEAR]
             if args.start_year <= y <= args.end_year]
    print(f"Output dir: {OUTDIR}")
    print(f"Years: {years}")
    print("Source: NCAR ERA5 on AWS")

    if s3fs is None and not args.dry_run:
        raise RuntimeError("s3fs is required unless --dry-run is used")

    fs = None if args.dry_run else s3fs.S3FileSystem(
        anon=True, default_fill_cache=False, default_cache_type="none")
    for year in years:
        if args.mode in ("cfram", "all"):
            download_year(fs, year, args.overwrite, args.dry_run)
        if args.mode in ("diag", "all"):
            download_diag_year(fs, year, args.overwrite, args.dry_run)


if __name__ == "__main__":
    main()
