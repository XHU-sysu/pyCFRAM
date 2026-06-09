#!/usr/bin/env python3
"""Download selected MERRA-2 aerosol dates from a plain date list."""

import argparse
from datetime import datetime
import os
import sys

from download_merra2_aerosol import download_day


VARIABLES = [
    'BCPHILIC', 'BCPHOBIC',
    'OCPHILIC', 'OCPHOBIC',
    'SO4',
    'SS001', 'SS002', 'SS003', 'SS004', 'SS005',
    'DU001', 'DU002', 'DU003', 'DU004', 'DU005',
    'DELP',
]


def read_dates(path):
    dates = []
    with open(path, 'r', encoding='utf-8') as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            dates.append(datetime.strptime(line, '%Y-%m-%d').date())
    return dates


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--date-file', required=True,
                    help='Text file with one YYYY-MM-DD date per line')
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--lat', type=float, nargs=2, required=True)
    ap.add_argument('--lon', type=float, nargs=2, required=True)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    dates = read_dates(args.date_file)
    print('=== MERRA-2 Missing-Date Fill ===')
    print(f'Date file: {args.date_file}')
    print(f'Dates: {len(dates)}')
    print(f'Domain: lat [{args.lat[0]}, {args.lat[1]}], '
          f'lon [{args.lon[0]}, {args.lon[1]}]')
    print(f'Output: {args.outdir}')
    print()

    failed = 0
    for dt in dates:
        ok = download_day(dt, args.outdir, VARIABLES, args.lat, args.lon)
        if not ok:
            failed += 1

    print(f'\n=== Done: {len(dates) - failed}/{len(dates)} files downloaded ===')
    if failed:
        print(f'WARNING: {failed} files failed')
        sys.exit(1)


if __name__ == '__main__':
    main()
