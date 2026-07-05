#!/usr/bin/env python3
"""Generate the synthetic DAMIP smoke fixtures consumed by
``tests/test_damip_source.py`` (Phase 3, WP-M4.2).

This is a one-off *generator*, not a pytest module (leading underscore so
pytest never collects it). Run it directly to (re)materialize the small
NetCDF fixtures checked into this directory:

    python3 tests/data/damip_smoke/_make_damip_fixtures.py

Two model directories + one standalone O3-climatology file are produced:

- ``ipsl_mini/``: a "full variable" DAMIP hist-aer dataset (ta, hus, ts, ps,
  cl/clw/cli on hybrid ap/b levels, rsdt, rsds, rsus, hfls, hfss, huss, o3)
  -- drives cloud=ACTIVE, o3=MODEL, solar=ACTIVE, flux=ACTIVE, albedo=ACTIVE.
  Also carries one deliberately "mountain" grid column (ps well below the
  full plev range) with genuine CMIP6-style fill-valued cells at the two
  lowest target pressure levels, to exercise the hybrid-projection
  subsurface-zeroing branch AND cmip6_common.fill_subsurface's HOLD/ts
  strategy end-to-end.
- ``cesm2_damip_mini/``: the real CESM2 hist-aer "missing cloud + missing
  rsdt + missing o3" situation (docs/plan_ph3.md §1.4) -- no cl/clw/cli,
  no rsdt, no o3 files at all -- drives cloud=SKIPPED, solar=ANALYTIC,
  o3=CLIMATOLOGY (when paired with ``o3_climatology_mini.nc``) or
  o3=SKIPPED (when o3_climatology_path points nowhere / o3=skip is forced).
- ``o3_climatology_mini.nc``: a tiny CESM-style prescribed-O3 climatology
  (same variable/dim names as the real
  ``raw_data/ozone_1.9x2.5_L26_1850clim_c090420.nc``: ``O3(time,lev,lat,lon)``
  mol/mol, ``lev`` ascending hPa, ``lat`` ascending degrees) for the
  o3=auto/climatology injection tests.

Both model datasets share one 72-month (1850-01..1855-12) synthetic time
series on a 4x4 lat/lon grid, with per-YEAR-constant values (i.e. every
month within a given year has the identical value) so that a day-weighted
annual climatology reduces, by construction, to a simple arithmetic mean
of the (small number of) distinct per-year constants selected -- this
keeps the expected numeric values easy to hand-verify in the test file
without re-deriving cmip6_common's own (separately tested in
tests/test_cmip6_common.py) day-weighting arithmetic.

base_years=[1850, 1851] and warm_years=[1854, 1855] are deliberately BOTH
non-leap-year pairs under the gregorian calendar (1852 -- the leap year in
this range -- is present in the file but excluded from both selected
windows), so each pair gets an exact 50/50 day-weight split and a linear
per-year trend produces an exactly-computable (not just approximately
equal) base/warm delta: e.g. a trend of +T per year gives
warm_mean - base_mean == T * ((4.5) - (0.5)) == 4.0*T exactly.
"""
import os

import cftime
import numpy as np
from netCDF4 import Dataset

HERE = os.path.dirname(os.path.abspath(__file__))

LAT = np.array([-60.0, -20.0, 20.0, 60.0])
LON = np.array([0.0, 90.0, 180.0, 270.0])
NLAT, NLON = len(LAT), len(LON)

# Case target grid (docs/plan_ph3.md §4.1 convention: TOA->sfc, hPa, ascending)
TARGET_PLEV_HPA_TOA2SFC = [100.0, 300.0, 500.0, 700.0, 1000.0]

# Native plev for ipsl_mini / cesm2_damip_mini (sfc->TOA, Pa, descending) --
# identical (in value) to the case target grid reversed x100, i.e. plev==target
# so build_states()'s interp_plev_to_target call is an identity (per its own
# documented guarantee) for the "happy path" fixtures. A dedicated
# differing-grid case is exercised separately in the test file using a case
# grid that is NOT identical to this, to touch the genuine-interpolation path.
NATIVE_PLEV_PA_SFC2TOA = np.array([100000.0, 70000.0, 50000.0, 30000.0, 10000.0])
NLEV_P = len(NATIVE_PLEV_PA_SFC2TOA)

# Mountain column (deliberately low ps -> exercises hybrid subsurface
# zero-fill AND cmip6_common.fill_subsurface's HOLD/ts strategy).
MOUNTAIN_J, MOUNTAIN_I = 0, 0
PS_FLAT_PA = 101300.0
PS_MOUNTAIN_PA = 60000.0   # 600 hPa -- well below the 700/1000 hPa target levels

FILLVAL = 1.0e20

# 6 years total; base_years=[1850,1851] (yo=0,1) and warm_years=[1854,1855]
# (yo=4,5) are the two selected windows (1852/1853, yo=2,3, are unused
# filler years -- present in the file, not selected by either test range).
YEARS = [1850, 1851, 1852, 1853, 1854, 1855]
BASE_YEARS = [1850, 1851]
WARM_YEARS = [1854, 1855]


def _time_axis(calendar, units):
    dates = []
    for y in YEARS:
        for m in range(1, 13):
            dates.append(cftime.datetime(y, m, 15, calendar=calendar))
    return np.asarray(cftime.date2num(dates, units=units, calendar=calendar))


def _year_const_3d(values_per_year, nlev, extra_fill=None):
    """Build a (ntime=48, nlev, nlat, nlon) array where every month within a
    given year carries the SAME value: ``values_per_year[year_offset]``
    (a (nlev, nlat, nlon) array), year_offset = 0..3.

    extra_fill: optional list of (lev_idx, j, i) cells to overwrite with
    FILLVAL for every month (simulates CMIP6 topography masking).
    """
    ntime = 12 * len(YEARS)
    out = np.empty((ntime, nlev, NLAT, NLON), dtype=np.float64)
    for yo in range(len(YEARS)):
        out[yo * 12:(yo + 1) * 12, :, :, :] = values_per_year[yo][None, :, :, :]
    if extra_fill:
        for (k, j, i) in extra_fill:
            out[:, k, j, i] = FILLVAL
    return out


def _year_const_2d(values_per_year):
    ntime = 12 * len(YEARS)
    out = np.empty((ntime, NLAT, NLON), dtype=np.float64)
    for yo in range(len(YEARS)):
        out[yo * 12:(yo + 1) * 12, :, :] = values_per_year[yo][None, :, :]
    return out


def _write_plev_var(path, varname, data_tnlevjnlon, time_vals, time_units, calendar,
                     plev_pa_sfc2toa=NATIVE_PLEV_PA_SFC2TOA):
    with Dataset(path, 'w') as nc:
        nc.createDimension('time', None)
        nc.createDimension('plev', len(plev_pa_sfc2toa))
        nc.createDimension('lat', NLAT)
        nc.createDimension('lon', NLON)
        tv = nc.createVariable('time', 'f8', ('time',))
        tv[:] = time_vals
        tv.units = time_units
        tv.calendar = calendar
        pv = nc.createVariable('plev', 'f8', ('plev',))
        pv[:] = plev_pa_sfc2toa
        pv.units = 'Pa'
        nc.createVariable('lat', 'f8', ('lat',))[:] = LAT
        nc.createVariable('lon', 'f8', ('lon',))[:] = LON
        v = nc.createVariable(varname, 'f8', ('time', 'plev', 'lat', 'lon'))
        v[:] = data_tnlevjnlon


def _write_surf_var(path, varname, data_tjnlon, time_vals, time_units, calendar):
    with Dataset(path, 'w') as nc:
        nc.createDimension('time', None)
        nc.createDimension('lat', NLAT)
        nc.createDimension('lon', NLON)
        tv = nc.createVariable('time', 'f8', ('time',))
        tv[:] = time_vals
        tv.units = time_units
        tv.calendar = calendar
        nc.createVariable('lat', 'f8', ('lat',))[:] = LAT
        nc.createVariable('lon', 'f8', ('lon',))[:] = LON
        v = nc.createVariable(varname, 'f8', ('time', 'lat', 'lon'))
        v[:] = data_tjnlon


def _write_hybrid_cloud_var(path, varname, data_tnlevjnlon, ap, b, time_vals, time_units, calendar):
    with Dataset(path, 'w') as nc:
        nc.createDimension('time', None)
        nc.createDimension('lev', len(ap))
        nc.createDimension('lat', NLAT)
        nc.createDimension('lon', NLON)
        tv = nc.createVariable('time', 'f8', ('time',))
        tv[:] = time_vals
        tv.units = time_units
        tv.calendar = calendar
        lv = nc.createVariable('lev', 'f8', ('lev',))
        lv[:] = np.arange(len(ap), dtype=np.float64) + 1.0
        lv.standard_name = 'atmosphere_hybrid_sigma_pressure_coordinate'
        lv.formula_terms = 'ap: ap b: b ps: ps'
        nc.createVariable('ap', 'f8', ('lev',))[:] = ap
        nc.createVariable('b', 'f8', ('lev',))[:] = b
        nc.createVariable('lat', 'f8', ('lat',))[:] = LAT
        nc.createVariable('lon', 'f8', ('lon',))[:] = LON
        v = nc.createVariable(varname, 'f8', ('time', 'lev', 'lat', 'lon'))
        v.formula_terms = 'ap: ap b: b ps: ps'
        v[:] = data_tnlevjnlon


def _ps_field(mountain=True):
    ps = np.full((NLAT, NLON), PS_FLAT_PA)
    if mountain:
        ps[MOUNTAIN_J, MOUNTAIN_I] = PS_MOUNTAIN_PA
    return ps


def make_ipsl_mini():
    outdir = os.path.join(HERE, 'ipsl_mini')
    os.makedirs(outdir, exist_ok=True)
    calendar = 'gregorian'
    time_units = 'days since 1850-01-01'
    time_vals = _time_axis(calendar, time_units)
    tag = 'Amon_IPSL-CM6A-LR_hist-aer_r1i1p1f1_gr'

    # ta: level0=1000hPa warmest .. level4=10hPa coldest; +0.6K/yr warming trend
    ta_per_year = []
    for yo in range(len(YEARS)):
        prof = np.array([292.0, 275.0, 255.0, 230.0, 205.0])  # per plev level, sfc->TOA
        field = prof[:, None, None] * np.ones((NLEV_P, NLAT, NLON))
        field += yo * 0.6   # uniform warming trend
        ta_per_year.append(field)
    # NATIVE_PLEV_PA_SFC2TOA = [100000(1000hPa), 70000(700hPa), 50000(500hPa),
    # 30000(300hPa), 10000(100hPa)] Pa. Mountain column ps=60000 Pa (600 hPa)
    # -> native levels 0 (1000hPa) and 1 (700hPa) are genuinely below-ground
    # (p > ps); levels 2-4 (500/300/100 hPa) are real atmosphere above the
    # mountain and must NOT be masked.
    ta_data = _year_const_3d(ta_per_year, NLEV_P,
                             extra_fill=[(0, MOUNTAIN_J, MOUNTAIN_I), (1, MOUNTAIN_J, MOUNTAIN_I)])
    _write_plev_var(os.path.join(outdir, 'ta_%s_185001-185512.nc' % tag), 'ta', ta_data,
                    time_vals, time_units, calendar)

    # hus: decreasing with height, time-invariant (no trend)
    hus_prof = np.array([1.0e-2, 5.0e-3, 2.0e-3, 5.0e-4, 1.0e-5])
    hus_field = hus_prof[:, None, None] * np.ones((NLEV_P, NLAT, NLON))
    hus_per_year = [hus_field.copy() for _ in YEARS]
    hus_data = _year_const_3d(hus_per_year, NLEV_P,
                              extra_fill=[(0, MOUNTAIN_J, MOUNTAIN_I), (1, MOUNTAIN_J, MOUNTAIN_I)])
    _write_plev_var(os.path.join(outdir, 'hus_%s_185001-185512.nc' % tag), 'hus', hus_data,
                    time_vals, time_units, calendar)

    # o3: constant across years (hist-aer freezes O3) -- mol/mol
    o3_prof = np.array([2.0e-8, 8.0e-8, 3.0e-7, 1.0e-6, 5.0e-6])  # sfc->TOA increasing with height
    o3_field = o3_prof[:, None, None] * np.ones((NLEV_P, NLAT, NLON))
    o3_per_year = [o3_field.copy() for _ in YEARS]
    o3_data = _year_const_3d(o3_per_year, NLEV_P)
    _write_plev_var(os.path.join(outdir, 'o3_%s_185001-185512.nc' % tag), 'o3', o3_data,
                    time_vals, time_units, calendar)

    # ts: modest warming trend, ps: static topography (mountain column)
    ts_per_year = [np.full((NLAT, NLON), 288.0) + yo * 0.6 for yo in range(len(YEARS))]
    _write_surf_var(os.path.join(outdir, 'ts_%s_185001-185512.nc' % tag), 'ts',
                    _year_const_2d(ts_per_year), time_vals, time_units, calendar)

    ps_field = _ps_field(mountain=True)
    ps_per_year = [ps_field.copy() for _ in YEARS]
    _write_surf_var(os.path.join(outdir, 'ps_%s_185001-185512.nc' % tag), 'ps',
                    _year_const_2d(ps_per_year), time_vals, time_units, calendar)

    # rsds/rsus: time-invariant -> albedo identical base/warm
    rsds_field = np.full((NLAT, NLON), 200.0)
    rsus_field = np.full((NLAT, NLON), 40.0)   # albedo = 0.2
    _write_surf_var(os.path.join(outdir, 'rsds_%s_185001-185512.nc' % tag), 'rsds',
                    _year_const_2d([rsds_field.copy() for _ in YEARS]), time_vals, time_units, calendar)
    _write_surf_var(os.path.join(outdir, 'rsus_%s_185001-185512.nc' % tag), 'rsus',
                    _year_const_2d([rsus_field.copy() for _ in YEARS]), time_vals, time_units, calendar)

    # rsdt: time-invariant (solar constant frozen in hist-aer) -> frc_solar==0
    rsdt_field = np.full((NLAT, NLON), 340.0)
    _write_surf_var(os.path.join(outdir, 'rsdt_%s_185001-185512.nc' % tag), 'rsdt',
                    _year_const_2d([rsdt_field.copy() for _ in YEARS]), time_vals, time_units, calendar)

    # huss: time-invariant, optional surf var
    huss_field = np.full((NLAT, NLON), 8.0e-3)
    _write_surf_var(os.path.join(outdir, 'huss_%s_185001-185512.nc' % tag), 'huss',
                    _year_const_2d([huss_field.copy() for _ in YEARS]), time_vals, time_units, calendar)

    # hfls/hfss: upward-positive, trend so warm > base (more evaporation/sensible loss)
    hfls_per_year = [np.full((NLAT, NLON), 80.0) + yo * 5.0 for yo in range(len(YEARS))]
    hfss_per_year = [np.full((NLAT, NLON), 20.0) + yo * 2.0 for yo in range(len(YEARS))]
    _write_surf_var(os.path.join(outdir, 'hfls_%s_185001-185512.nc' % tag), 'hfls',
                    _year_const_2d(hfls_per_year), time_vals, time_units, calendar)
    _write_surf_var(os.path.join(outdir, 'hfss_%s_185001-185512.nc' % tag), 'hfss',
                    _year_const_2d(hfss_per_year), time_vals, time_units, calendar)

    # cl/clw/cli: hybrid ap/b levels (CMOR mainstream convention), 4 levels,
    # TOA->sfc. Cloud fraction (%) decreasing trend (aerosol reduces cloud).
    ap = np.array([1000.0, 8000.0, 15000.0, 0.0])
    b = np.array([0.0, 0.05, 0.35, 0.98])
    cl_prof_pct = np.array([2.0, 15.0, 40.0, 20.0])   # % by hybrid level, TOA->sfc
    clw_prof = np.array([0.0, 1.0e-6, 5.0e-5, 2.0e-5])
    cli_prof = np.array([2.0e-6, 3.0e-5, 1.0e-5, 0.0])

    cl_per_year = [(cl_prof_pct[:, None, None] * np.ones((4, NLAT, NLON))) - yo * 1.0
                   for yo in range(len(YEARS))]
    clw_per_year = [(clw_prof[:, None, None] * np.ones((4, NLAT, NLON))) * (1.0 - 0.05 * yo)
                    for yo in range(len(YEARS))]
    cli_per_year = [(cli_prof[:, None, None] * np.ones((4, NLAT, NLON))) * (1.0 - 0.05 * yo)
                    for yo in range(len(YEARS))]

    cl_data = _year_const_3d(cl_per_year, 4)
    clw_data = _year_const_3d(clw_per_year, 4)
    cli_data = _year_const_3d(cli_per_year, 4)
    _write_hybrid_cloud_var(os.path.join(outdir, 'cl_%s_185001-185512.nc' % tag), 'cl', cl_data,
                            ap, b, time_vals, time_units, calendar)
    _write_hybrid_cloud_var(os.path.join(outdir, 'clw_%s_185001-185512.nc' % tag), 'clw', clw_data,
                            ap, b, time_vals, time_units, calendar)
    _write_hybrid_cloud_var(os.path.join(outdir, 'cli_%s_185001-185512.nc' % tag), 'cli', cli_data,
                            ap, b, time_vals, time_units, calendar)

    print('Wrote ipsl_mini/ (%d files)' % len(os.listdir(outdir)))


def make_cesm2_damip_mini():
    outdir = os.path.join(HERE, 'cesm2_damip_mini')
    os.makedirs(outdir, exist_ok=True)
    calendar = 'noleap'
    time_units = 'days since 1850-01-01'
    time_vals = _time_axis(calendar, time_units)
    tag = 'Amon_CESM2_hist-aer_r1i1p1f1_gn'

    ta_per_year = []
    for yo in range(len(YEARS)):
        prof = np.array([290.0, 273.0, 253.0, 228.0, 203.0])
        field = prof[:, None, None] * np.ones((NLEV_P, NLAT, NLON))
        field += yo * 0.4
        ta_per_year.append(field)
    ta_data = _year_const_3d(ta_per_year, NLEV_P)   # no mountain column here (already
    _write_plev_var(os.path.join(outdir, 'ta_%s_185001-185512.nc' % tag), 'ta', ta_data,
                    time_vals, time_units, calendar)                        # covered by ipsl_mini)

    hus_prof = np.array([9.0e-3, 4.5e-3, 1.8e-3, 4.0e-4, 8.0e-6])
    hus_field = hus_prof[:, None, None] * np.ones((NLEV_P, NLAT, NLON))
    hus_per_year = [hus_field.copy() for _ in YEARS]
    hus_data = _year_const_3d(hus_per_year, NLEV_P)
    _write_plev_var(os.path.join(outdir, 'hus_%s_185001-185512.nc' % tag), 'hus', hus_data,
                    time_vals, time_units, calendar)

    ts_per_year = [np.full((NLAT, NLON), 285.0) + yo * 0.4 for yo in range(len(YEARS))]
    _write_surf_var(os.path.join(outdir, 'ts_%s_185001-185512.nc' % tag), 'ts',
                    _year_const_2d(ts_per_year), time_vals, time_units, calendar)

    ps_field = _ps_field(mountain=False)
    _write_surf_var(os.path.join(outdir, 'ps_%s_185001-185512.nc' % tag), 'ps',
                    _year_const_2d([ps_field.copy() for _ in YEARS]), time_vals, time_units, calendar)

    rsds_field = np.full((NLAT, NLON), 190.0)
    rsus_field = np.full((NLAT, NLON), 55.0)
    _write_surf_var(os.path.join(outdir, 'rsds_%s_185001-185512.nc' % tag), 'rsds',
                    _year_const_2d([rsds_field.copy() for _ in YEARS]), time_vals, time_units, calendar)
    _write_surf_var(os.path.join(outdir, 'rsus_%s_185001-185512.nc' % tag), 'rsus',
                    _year_const_2d([rsus_field.copy() for _ in YEARS]), time_vals, time_units, calendar)

    huss_field = np.full((NLAT, NLON), 7.0e-3)
    _write_surf_var(os.path.join(outdir, 'huss_%s_185001-185512.nc' % tag), 'huss',
                    _year_const_2d([huss_field.copy() for _ in YEARS]), time_vals, time_units, calendar)

    hfls_per_year = [np.full((NLAT, NLON), 70.0) + yo * 4.0 for yo in range(len(YEARS))]
    hfss_per_year = [np.full((NLAT, NLON), 18.0) + yo * 1.5 for yo in range(len(YEARS))]
    _write_surf_var(os.path.join(outdir, 'hfls_%s_185001-185512.nc' % tag), 'hfls',
                    _year_const_2d(hfls_per_year), time_vals, time_units, calendar)
    _write_surf_var(os.path.join(outdir, 'hfss_%s_185001-185512.nc' % tag), 'hfss',
                    _year_const_2d(hfss_per_year), time_vals, time_units, calendar)

    # Deliberately NO cl/clw/cli, NO rsdt, NO o3 files -- matches the real
    # CESM2 hist-aer availability gap (docs/plan_ph3.md §1.4).
    print('Wrote cesm2_damip_mini/ (%d files, no cl/clw/cli/rsdt/o3 by design)' %
          len(os.listdir(outdir)))


def make_o3_climatology_mini():
    path = os.path.join(HERE, 'o3_climatology_mini.nc')
    nlev, nlat, nlon, ntime = 5, 8, 4, 12
    lev = np.array([10.0, 100.0, 300.0, 700.0, 1000.0])    # ascending hPa (matches real file convention)
    lat = np.linspace(-90.0, 90.0, nlat)                    # ascending
    lon = np.linspace(0.0, 270.0, nlon)
    o3_prof = np.array([5.0e-6, 1.0e-6, 3.0e-7, 8.0e-8, 3.0e-8])  # mol/mol, TOA->sfc... wait
    # NB: lev is ascending pressure (TOA->sfc), so index0(10hPa)=stratosphere
    # (O3 maximum), index-1(1000hPa)=surface (O3 minimum) -- matches o3_prof order above.
    with Dataset(path, 'w') as nc:
        nc.createDimension('time', ntime)
        nc.createDimension('lev', nlev)
        nc.createDimension('lat', nlat)
        nc.createDimension('lon', nlon)
        nc.createVariable('lev', 'f8', ('lev',))[:] = lev
        nc.createVariable('lat', 'f8', ('lat',))[:] = lat
        nc.createVariable('lon', 'f8', ('lon',))[:] = lon
        o3 = np.broadcast_to(o3_prof[None, :, None, None], (ntime, nlev, nlat, nlon)).copy()
        nc.createVariable('O3', 'f8', ('time', 'lev', 'lat', 'lon'))[:] = o3
    print('Wrote o3_climatology_mini.nc')


if __name__ == '__main__':
    make_ipsl_mini()
    make_cesm2_damip_mini()
    make_o3_climatology_mini()
