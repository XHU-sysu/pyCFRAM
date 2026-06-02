#!/usr/bin/env python3
"""Download ERA5 6-hourly data for India-Bangladesh extreme wet heat case.

Event   : 2023-04-17 ~ 2023-04-20 (extreme wet heat core period)
Clim    : April 1991-2020 (background base state)
Domain  : 10-35°N, 65-110°E

Output layout (GRIB, no CDS-side NetCDF conversion):
  OUTDIR/
    era5_pl_all_{YYYY04}.grib              # all PL variables per year-month
    era5_sl_all_{YYYY04}.grib              # all SL variables per year-month

Strategy:
  Requests are downloaded as GRIB to avoid expensive CDS-side NetCDF conversion.
  1 PL request + 1 SL request = 2 requests/year.

Run on hqlx220 (nohup friendly):
    python3 -u scripts/download_era5_india23.py 2>&1 | tee /tmp/dl_india23.log
    python3 -u scripts/download_era5_india23.py --dry-run
    python3 -u scripts/download_era5_india23.py --start-year 1991 --end-year 2000
"""

import os
import argparse
import cdsapi

# ── Configuration ─────────────────────────────────────────────────────────────

OUTDIR = '/disk/r074/lzhenn/workspace/ust-jumper/pyCFRAM/era5_data/india_wb23'

CDS_URL = 'https://cds.climate.copernicus.eu/api'
CDS_KEY = 'b5e21614-ef11-4117-9f80-ec6294fb65bd'

AREA        = [35, 65, 10, 110]              # N, W, S, E
MONTH       = '04'                           # April
CLIM_YEARS  = list(range(1991, 2021))        # 1991-2020
EVENT_YEAR  = 2023
ALL_DAYS    = [f'{d:02d}' for d in range(1, 31)]   # April = 30 days
TIMES_6H    = ['00:00', '06:00', '12:00', '18:00']

# All 37 ERA5 standard pressure levels (hPa)
PRESSURE_LEVELS = [
    '1','2','3','5','7','10','20','30','50','70',
    '100','125','150','175','200','225','250','300',
    '350','400','450','500','550','600','650','700',
    '750','775','800','825','850','875','900','925',
    '950','975','1000',
]

# PL: short name → CDS long name
PL_VARS = {
    't':    'temperature',
    'q':    'specific_humidity',
    'o3':   'ozone_mass_mixing_ratio',
    'cc':   'fraction_of_cloud_cover',
    'clwc': 'specific_cloud_liquid_water_content',
    'ciwc': 'specific_cloud_ice_water_content',
}

SL_VARS = [
    'skin_temperature',
    'surface_pressure',
    'surface_solar_radiation_downwards',
    'surface_net_solar_radiation',
    'toa_incident_solar_radiation',
    'surface_latent_heat_flux',
    'surface_sensible_heat_flux',
    'maximum_2m_temperature_since_previous_post_processing',
]


# ── Helpers ───────────────────────────────────────────────────────────────────

def make_client():
    return cdsapi.Client(url=CDS_URL, key=CDS_KEY, quiet=False)


def _skip(path, dry_run):
    if os.path.exists(path):
        sz = os.path.getsize(path)
        print(f'  SKIP  {os.path.basename(path)}  ({sz/1e6:.1f} MB, exists)')
        return True
    if dry_run:
        print(f'  DRY   {path}')
        return True
    return False


def _is_queue_limit_error(exc):
    msg = str(exc)
    return (
        'Number queued requests for this dataset is temporarily limited' in msg
        or 'temporarily limited' in msg
    )


def download_pl_all(c, year, dry_run):
    """Download all PL variables in one CDS request as GRIB."""
    yyyymm = f'{year}{MONTH}'
    outfile = os.path.join(OUTDIR, f'era5_pl_all_{yyyymm}.grib')

    if _skip(outfile, dry_run):
        return

    cds_names = list(PL_VARS.values())
    print(f'  → PL all {yyyymm} (GRIB, {len(cds_names)} vars, 1 request) ...',
          flush=True)
    c.retrieve(
        'reanalysis-era5-pressure-levels',
        {
            'product_type'    : 'reanalysis',
            'variable'        : cds_names,
            'pressure_level'  : PRESSURE_LEVELS,
            'year'            : str(year),
            'month'           : MONTH,
            'day'             : ALL_DAYS,
            'time'            : TIMES_6H,
            'area'            : AREA,
            'data_format'     : 'grib',
            'download_format' : 'unarchived',
        },
        outfile,
    )
    print(f'     done  {os.path.getsize(outfile)/1e6:.1f} MB', flush=True)


def download_sl_all(c, year, dry_run):
    """Download all SL variables in one CDS request as GRIB."""
    yyyymm  = f'{year}{MONTH}'
    outfile = os.path.join(OUTDIR, f'era5_sl_all_{yyyymm}.grib')
    if _skip(outfile, dry_run):
        return
    print(f'  → SL all {yyyymm} (GRIB, {len(SL_VARS)} vars, 1 request) ...',
          flush=True)
    c.retrieve(
        'reanalysis-era5-single-levels',
        {
            'product_type'    : 'reanalysis',
            'variable'        : SL_VARS,
            'year'            : str(year),
            'month'           : MONTH,
            'day'             : ALL_DAYS,
            'time'            : TIMES_6H,
            'area'            : AREA,
            'data_format'     : 'grib',
            'download_format' : 'unarchived',
        },
        outfile,
    )
    print(f'     done  {os.path.getsize(outfile)/1e6:.1f} MB', flush=True)


# ── Main ──────────────────────────────────────────────────────────────────────

def main(dry_run=False, start_year=None, end_year=None):
    if start_year is None:
        start_year = min(CLIM_YEARS)
    if end_year is None:
        end_year = EVENT_YEAR

    if not dry_run:
        os.makedirs(OUTDIR, exist_ok=True)
    print(f'Output dir : {OUTDIR}')
    print(f'Dry run    : {dry_run}')
    print(f'Year range : {start_year} – {end_year}')

    all_years = [y for y in CLIM_YEARS + [EVENT_YEAR]
                 if start_year <= y <= end_year]
    print(f'Years      : {all_years}')
    print('Requests/yr: 1 PL + 1 SL = 2 total')

    c = None if dry_run else make_client()

    n_total = len(all_years)
    for i, year in enumerate(all_years):
        label = 'EVENT' if year == EVENT_YEAR else 'clim'
        print(f'\n=== {year} [{label}] ({i+1}/{n_total}) ===', flush=True)

        try:
            download_pl_all(c, year, dry_run)
        except Exception as e:
            print(f'  ERROR PL all {year}: {e}', flush=True)
            if _is_queue_limit_error(e):
                print('  STOP  CDS dataset queue limit reached; retry later.',
                      flush=True)
                raise SystemExit(2)

        try:
            download_sl_all(c, year, dry_run)
        except Exception as e:
            print(f'  ERROR SL all {year}: {e}', flush=True)
            if _is_queue_limit_error(e):
                print('  STOP  CDS dataset queue limit reached; retry later.',
                      flush=True)
                raise SystemExit(2)

    print('\n=== All downloads complete ===')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--dry-run', action='store_true',
                    help='Print paths without downloading')
    ap.add_argument('--start-year', type=int, default=min(CLIM_YEARS),
                    help=f'First year to download (default: {min(CLIM_YEARS)})')
    ap.add_argument('--end-year', type=int, default=EVENT_YEAR,
                    help=f'Last year to download inclusive (default: {EVENT_YEAR})')
    args = ap.parse_args()
    main(dry_run=args.dry_run,
         start_year=args.start_year,
         end_year=args.end_year)
