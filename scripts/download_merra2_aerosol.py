#!/usr/bin/env python3
"""Download MERRA-2 3D aerosol mixing ratios (M2I3NVAER) via OPeNDAP.

Downloads regional subset for pyCFRAM aerosol input.
Requires NASA Earthdata credentials in ~/.netrc:
    machine urs.earthdata.nasa.gov
    login YOUR_USERNAME
    password YOUR_PASSWORD

Usage:
    python scripts/download_merra2_aerosol.py --years 2003 2022 --month 8 \
        --outdir era5_data/merra2

    python scripts/download_merra2_aerosol.py \
        --year-list 1991 1992 1993 2023 --month 4 \
        --outdir era5_data/india_wb23_merra2 --lat 8 37 --lon 63 112

Output:
    era5_data/merra2/M2I3NVAER_{YYYYMMDD}.nc4  (one file per day)
"""

import os
import sys
import argparse
import calendar
import time
from datetime import date


def _merra2_stream(year):
    """Get MERRA-2 file stream number based on year.

    MERRA-2 uses different stream IDs for different periods:
      100: 1980-1991, 200: 1992-2000, 300: 2001-2010, 400: 2011-2020, 2022+
      401: 2021 (reprocessed subset)
    """
    if year <= 1991:
        return 100
    elif year <= 2000:
        return 200
    elif year <= 2010:
        return 300
    elif year == 2021:
        return 401
    else:
        return 400


def build_opendap_url(dt, variables, lat_range, lon_range):
    """Build OPeNDAP subset URL for one day of M2I3NVAER.

    Args:
        dt: date object
        variables: list of variable names
        lat_range: (lat_min, lat_max) in degrees
        lon_range: (lon_min, lon_max) in degrees
    """
    base = ('https://goldsmr5.gesdisc.eosdis.nasa.gov/opendap/MERRA2/'
            'M2I3NVAER.5.12.4')
    stream = _merra2_stream(dt.year)
    fname = f'MERRA2_{stream}.inst3_3d_aer_Nv.{dt.strftime("%Y%m%d")}.nc4'
    url = f'{base}/{dt.year:04d}/{dt.month:02d}/{fname}'

    # OPeNDAP subsetting: variable list + lat/lon constraints
    # MERRA-2 grid: lat -90 to 90 (0.5°), lon -180 to 180 (0.625°)
    # Indices: lat[i] = -90 + i*0.5, i=0..360  → lat_idx = (lat + 90) / 0.5
    # Indices: lon[j] = -180 + j*0.625, j=0..575 → lon_idx = (lon + 180) / 0.625
    lat_min_idx = int((lat_range[0] + 90) / 0.5)
    lat_max_idx = int((lat_range[1] + 90) / 0.5)
    lon_min_idx = int((lon_range[0] + 180) / 0.625)
    lon_max_idx = int((lon_range[1] + 180) / 0.625)

    # Build constraint expression
    # dims: time[0:7], lev[0:71], lat[idx_min:idx_max], lon[idx_min:idx_max]
    dim_str = f'[0:1:7][0:1:71][{lat_min_idx}:1:{lat_max_idx}][{lon_min_idx}:1:{lon_max_idx}]'
    var_str = ','.join(f'{v}{dim_str}' for v in variables)

    # Also need lat, lon, lev, time as 1D
    coord_str = (f'time[0:1:7],lev[0:1:71],'
                 f'lat[{lat_min_idx}:1:{lat_max_idx}],'
                 f'lon[{lon_min_idx}:1:{lon_max_idx}]')

    return f'{url}.nc4?{var_str},{coord_str}'


def download_day(dt, outdir, variables, lat_range, lon_range):
    """Download one day of MERRA-2 aerosol data."""
    outfile = os.path.join(outdir, f'M2I3NVAER_{dt.strftime("%Y%m%d")}.nc4')
    MIN_SIZE = 10 * 1024 * 1024  # 10 MB — valid MERRA-2 files are ~30-60 MB
    MAX_ATTEMPTS = 4
    if os.path.exists(outfile) and os.path.getsize(outfile) > MIN_SIZE:
        print(f'  {dt}: skip (exists, {os.path.getsize(outfile)/1e6:.1f} MB)')
        return True

    url = build_opendap_url(dt, variables, lat_range, lon_range)

    # Use per-process cookie file to avoid conflicts when running multiple instances
    cookie_file = os.path.expanduser(f'~/.urs_cookies_{os.getpid()}')

    # Use curl with Earthdata auth (.netrc) and -g to disable URL globbing
    cmd = (f'curl -s -g -n -c {cookie_file} -b {cookie_file} -L '
           f'-o "{outfile}" "{url}"')

    cmd2 = (f'wget -q --auth-no-challenge '
            f'--load-cookies {cookie_file} --save-cookies {cookie_file} '
            f'--keep-session-cookies '
            f'-O "{outfile}" "{url}"')

    ret = 1
    for attempt in range(1, MAX_ATTEMPTS + 1):
        print(f'  {dt}: downloading... attempt {attempt}/{MAX_ATTEMPTS}',
              flush=True)
        ret = os.system(cmd)
        if ret != 0:
            ret = os.system(cmd2)

        if os.path.exists(outfile) and os.path.getsize(outfile) >= MIN_SIZE:
            break

        if os.path.exists(outfile):
            os.remove(outfile)
        if attempt < MAX_ATTEMPTS:
            time.sleep(10 * attempt)

    if ret != 0 or not os.path.exists(outfile) or os.path.getsize(outfile) < MIN_SIZE:
        print(f'  {dt}: FAILED after {MAX_ATTEMPTS} attempts (ret={ret})')
        if os.path.exists(outfile):
            os.remove(outfile)
        return False

    fsize = os.path.getsize(outfile) / 1e6
    print(f'  {dt}: done ({fsize:.1f} MB)')
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--years', type=int, nargs=2, default=[2003, 2022],
                        help='Year range (inclusive)')
    parser.add_argument('--year-list', type=int, nargs='+', default=None,
                        help='Explicit years to download; overrides --years')
    parser.add_argument('--month', type=int, default=8, help='Month to download')
    parser.add_argument('--outdir', default='era5_data/merra2',
                        help='Output directory')
    parser.add_argument('--lat', type=float, nargs=2, default=[17, 43],
                        help='Latitude range (slightly larger than ERA5 domain)')
    parser.add_argument('--lon', type=float, nargs=2, default=[94, 126],
                        help='Longitude range')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # MERRA-2 aerosol variables needed for GOCART mapping
    variables = [
        'BCPHILIC', 'BCPHOBIC',     # Black carbon
        'OCPHILIC', 'OCPHOBIC',     # Organic carbon
        'SO4',                       # Sulfate
        'SS001', 'SS002', 'SS003', 'SS004', 'SS005',  # Sea salt (5 bins)
        'DU001', 'DU002', 'DU003', 'DU004', 'DU005',  # Dust (5 bins)
        'DELP',                      # Layer pressure thickness (for vertical coord)
    ]

    if args.year_list:
        years = sorted(set(args.year_list))
        yr_start, yr_end = years[0], years[-1]
    else:
        yr_start, yr_end = args.years
        years = list(range(yr_start, yr_end + 1))
    total = 0
    failed = 0

    print(f'=== MERRA-2 M2I3NVAER Download ===')
    if args.year_list:
        print(f'Years: {years}, Month: {args.month}')
    else:
        print(f'Years: {yr_start}-{yr_end}, Month: {args.month}')
    print(f'Domain: lat [{args.lat[0]}, {args.lat[1]}], lon [{args.lon[0]}, {args.lon[1]}]')
    print(f'Variables: {len(variables)}')
    print(f'Output: {args.outdir}')
    print()

    for year in years:
        ndays = calendar.monthrange(year, args.month)[1]
        print(f'Year {year} ({ndays} days):')
        for day in range(1, ndays + 1):
            dt = date(year, args.month, day)
            total += 1
            ok = download_day(dt, args.outdir, variables,
                              args.lat, args.lon)
            if not ok:
                failed += 1

    print(f'\n=== Done: {total - failed}/{total} files downloaded ===')
    if failed:
        print(f'WARNING: {failed} files failed')
        sys.exit(1)


if __name__ == '__main__':
    main()
