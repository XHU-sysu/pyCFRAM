#!/usr/bin/env python3
"""Download buffered May-July ERA5 inputs for the Alaska wildfire case.

The downloader reads NCAR's public ERA5 archive on AWS and writes the same
regional NetCDF layout consumed by pyCFRAM's ERA5 interface.
"""

import argparse

import download_era5_india23_ncar as ncar


OUTDIR = (
    "/disk/r074/lzhenn/workspace/ust-jumper/pyCFRAM/"
    "era5_data/alaska_wf22_ncar"
)
YEARS = list(range(1991, 2021)) + [2022]
MONTHS = [5, 6, 7]

# Buffer all documented Alaska case regions, including the Arctic response box.
LAT_NORTH = 78.0
LAT_SOUTH = 55.0
LON_WEST = 180.0
LON_EAST = 245.0

T2M_VAR = {"t2m": ncar.DIAG_INSTANT_VARS["t2m"]}


def parse_int_list(value):
    return [int(item) for item in value.split(",") if item.strip()]


def configure(month, outdir):
    ncar.MONTH = month
    ncar.OUTDIR = outdir
    ncar.LAT_NORTH = LAT_NORTH
    ncar.LAT_SOUTH = LAT_SOUTH
    ncar.LON_WEST = LON_WEST
    ncar.LON_EAST = LON_EAST


def download_month(fs, year, month, outdir, overwrite=False, dry_run=False):
    configure(month, outdir)
    print(f"\n=== Alaska ERA5 {year}-{month:02d} ===", flush=True)
    for var_short in ncar.PL_VARS:
        ncar.download_pl_var(fs, year, var_short, overwrite, dry_run)
    ncar.download_instant_sl(fs, year, overwrite, dry_run)
    ncar.download_accum_sl(fs, year, overwrite, dry_run)
    ncar.download_instant_sl(
        fs,
        year,
        overwrite,
        dry_run,
        var_table=T2M_VAR,
        out_name="data_stream-oper_stepType-instant_diag.nc",
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--years", type=parse_int_list, default=YEARS)
    ap.add_argument("--months", type=parse_int_list, default=MONTHS)
    ap.add_argument("--outdir", default=OUTDIR)
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    invalid_years = sorted(set(args.years) - set(YEARS))
    invalid_months = sorted(set(args.months) - set(MONTHS))
    if invalid_years or invalid_months:
        ap.error(
            f"outside Alaska case window: years={invalid_years}, "
            f"months={invalid_months}"
        )

    print(f"Output dir: {args.outdir}")
    print(f"Years: {args.years}")
    print(f"Months: {args.months}")
    print(
        f"Buffered domain: {LAT_SOUTH:g}-{LAT_NORTH:g}N, "
        f"{LON_WEST:g}-{LON_EAST:g}E"
    )
    print("Variables: CFRAM PL/SL inputs + t2m")
    print("Source: NCAR ERA5 on AWS")

    if ncar.s3fs is None and not args.dry_run:
        raise RuntimeError("s3fs is required unless --dry-run is used")
    fs = None if args.dry_run else ncar.s3fs.S3FileSystem(
        anon=True, default_fill_cache=False, default_cache_type="none"
    )

    for year in args.years:
        for month in args.months:
            download_month(
                fs, year, month, args.outdir, args.overwrite, args.dry_run
            )


if __name__ == "__main__":
    main()
