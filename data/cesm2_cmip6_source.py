"""CMIP6 CESM2 raw output reader + hybrid→plev re-projection.

Loads piControl / abrupt-4xCO2 monthly NetCDF, computes annual climatology
across selected year range, and re-projects cl/clw/cli from 32-layer hybrid
sigma-pressure to user-supplied plev (e.g. CMIP6 plev19).

Convention (input files, CMIP6 standard)
----------------------------------------
- plev (ta/hus): surface→TOA in NetCDF, units Pa
- hybrid (cl/clw/cli): TOA→sfc in NetCDF; pressure of layer k =
                       a[k]*p0 + b[k]*ps(lat,lon)
- lat: -90→90 (S→N), lon: 0→358.75 (W→E)
- time: monthly mean, calendar=noleap, time_origin=0001-01-01

Convention (output, pyCFRAM input)
----------------------------------
NetCDF stores arrays in **surface→TOA** lev order (same as CMIP6 plev). pyCFRAM
applies `[::-1]` internally to flip to TOA→sfc for Fortran processing.

Phase 3 refactor (docs/plan_ph3.md §2.2/§4, WP-M4.1)
-----------------------------------------------------
This module is now a **thin CESM2-specific shim** over the model-agnostic
machinery in `data/cmip6_common.py`: calendar handling (noleap), hybrid
coefficient naming (`a,b,p0`), and the mass-conserving hybrid→plev /
albedo numerical kernels are delegated there. This file only supplies the
CESM2-specific defaults (calendar, time units, filename layout, variable
list) that a future `data/cmip6_damip_source.py` (WP-M4.2) will supply
differently per model.

**Regression-safety contract**: `load_climo_pres`, `hybrid_to_plev_mass_conserving`,
and `compute_albedo` — the three functions `scripts/build_cesm2_official.py`
imports — must produce numerically identical output to the pre-refactor
version of this file. This is enforced by `tests/test_damip_regression.py`
(the "回归金标" gate, docs/plan_ph3.md §8 WP-M4.1) against a synthetic mini
CESM2-style fixture in `tests/data/damip_smoke/cesm2_mini/`. The real
CESM2 raw data only lives on hqlx210 (see persistent_context.md); a
matching full-scale check (rerun `cesm2_4xco2_official --step build`, diff
NC output) must be performed separately by whoever has remote access — see
that test file's module docstring.
"""
import os
import glob
import numpy as np
from netCDF4 import Dataset

from data import cmip6_common as common
from data.cmip6_common import hybrid_to_plev_mass_conserving, compute_albedo  # re-export, moved as-is


# CMIP6 plev (surface→TOA in NetCDF) for ta/hus — fixed
PLEV19_PA = np.array([100000, 92500, 85000, 70000, 60000, 50000, 40000,
                      30000, 25000, 20000, 15000, 10000, 7000, 5000,
                      3000, 2000, 1000, 500, 100], dtype=np.float64)

# pyCFRAM TOA→sfc convention (reverse of CMIP6 NetCDF order, in hPa)
PLEV19_HPA_TOP_DOWN = (PLEV19_PA[::-1] / 100.0).copy()  # [1, 5, ..., 925, 1000]

# CESM2-specific calendar/time defaults (used only as a fallback if a raw
# file is missing the `units`/`calendar` attributes — real CESM2 CMIP6
# output always carries both).
CESM2_CALENDAR = 'noleap'
CESM2_TIME_UNITS = 'days since 0001-01-01'


def list_files(raw_dir, exp_subdir):
    """Discover all CMIP6 monthly NetCDFs in raw_dir/exp_subdir."""
    files = sorted(glob.glob(os.path.join(raw_dir, exp_subdir, '*.nc')))
    var_files = {}
    for f in files:
        # filename: <var>_Amon_CESM2_<exp>_r1i1p1f1_gn_<period>.nc
        var = os.path.basename(f).split('_')[0]
        var_files[var] = f
    return var_files


def years_to_month_indices(time_var, year_start, year_end):
    """Return indices of months whose YEAR is in [year_start, year_end]
    inclusive.

    Delegates to `cmip6_common.decode_time` (calendar-aware via cftime),
    using the time variable's own `units`/`calendar` NetCDF attributes if
    present, else CESM2's noleap/'days since 0001-01-01' defaults (real
    CESM2 CMIP6 raw output always carries both attributes, so the fallback
    only matters for hand-built test fixtures). Replaces the pre-Phase-3
    noleap-only `days/365.0` arithmetic (moved to
    `data/cmip6_common.decode_time` and generalized to any calendar).
    """
    units = getattr(time_var, 'units', CESM2_TIME_UNITS)
    calendar = getattr(time_var, 'calendar', CESM2_CALENDAR)
    years, _months, _day_weights = common.decode_time(
        np.asarray(time_var[:]), units, calendar)
    return common.years_to_month_indices(years, year_start, year_end)


def _day_weights_for(time_var, time_indices):
    """CESM2 shim helper: decode the full time axis, return day_weights
    selected at `time_indices` (aligned 1:1, as required by
    `cmip6_common.annual_climo_from_monthly`)."""
    units = getattr(time_var, 'units', CESM2_TIME_UNITS)
    calendar = getattr(time_var, 'calendar', CESM2_CALENDAR)
    _years, _months, day_weights = common.decode_time(
        np.asarray(time_var[:]), units, calendar)
    return day_weights[time_indices]


def annual_climo_from_monthly(field, time_indices, time_var=None, day_weights=None):
    """CESM2 shim: day-weighted annual climatology, delegating the numerics
    to `cmip6_common.annual_climo_from_monthly`.

    Callers may pass either `time_var` (the netCDF4 time Variable, from
    which day_weights are derived via whichever `calendar` attribute it
    carries — noleap for real CESM2) or precomputed `day_weights` (already
    aligned to `time_indices`) directly.
    """
    if day_weights is None:
        if time_var is None:
            raise ValueError('annual_climo_from_monthly requires time_var or day_weights')
        day_weights = _day_weights_for(time_var, time_indices)
    return common.annual_climo_from_monthly(field, time_indices, day_weights)


def load_climo_pres(raw_dir, exp_subdir, year_start, year_end):
    """Load all pres-level + 2D variables, return (data_dict, lat, lon, plev).

    Returns dict with keys: ta, hus, cl, clw, cli, ts, ps, rsdt, rsds, rsus,
                            hfls, hfss, plus hybrid: a, b, p0
    Each value is 3D (lev, lat, lon) for upper-air, 2D (lat, lon) for surface.
    """
    files = list_files(raw_dir, exp_subdir)
    print('  Loading climo year %d-%d from %s/' % (year_start, year_end, exp_subdir))

    # Time indices come from any variable (all aligned)
    f = Dataset(files['ta'])
    time_var = f.variables['time']
    idx = years_to_month_indices(time_var, year_start, year_end)
    day_weights = _day_weights_for(time_var, idx)
    print('  selected %d months' % len(idx))
    lat = np.array(f.variables['lat'][:])
    lon = np.array(f.variables['lon'][:])
    plev = np.array(f.variables['plev'][:])
    f.close()

    out = {'lat': lat, 'lon': lon, 'plev_pres': plev}

    # Variables on plev grid (ta, hus): shape (nlev_plev=19, nlat, nlon)
    for var in ('ta', 'hus'):
        f = Dataset(files[var])
        out[var] = common.annual_climo_from_monthly(f.variables[var], idx, day_weights)
        f.close()

    # Variables on hybrid grid (cl, clw, cli): shape (nlev_hyb=32, nlat, nlon)
    f = Dataset(files['cl'])
    out['hybrid_a'] = np.array(f.variables['a'][:], dtype=np.float64)
    out['hybrid_b'] = np.array(f.variables['b'][:], dtype=np.float64)
    out['hybrid_p0'] = float(f.variables['p0'][...])
    out['cl'] = common.annual_climo_from_monthly(f.variables['cl'], idx, day_weights) / 100.0  # %→fraction
    f.close()

    f = Dataset(files['clw'])
    out['clw'] = common.annual_climo_from_monthly(f.variables['clw'], idx, day_weights)
    f.close()

    f = Dataset(files['cli'])
    out['cli'] = common.annual_climo_from_monthly(f.variables['cli'], idx, day_weights)
    f.close()

    # 2D surface fields. huss = 2m specific humidity (CMIP6 standard, kg/kg);
    # used by Fu RT for ph(nv1) — apple-to-apple OLD CFRAM raw/CFRAM.zip
    # GW-base.f L322-330 reads huss_base.dat for the surface row of /atmosp/.
    for var in ('ts', 'ps', 'rsdt', 'rsds', 'rsus', 'hfls', 'hfss', 'huss'):
        f = Dataset(files[var])
        out[var] = common.annual_climo_from_monthly(f.variables[var], idx, day_weights)
        f.close()

    return out


def hybrid_to_plev(field_hyb, a, b, p0, ps_2d, plev_target_pa):
    """Re-project field on hybrid sigma-pressure to fixed pressure levels.

    Args:
        field_hyb: shape (nlev_hyb, nlat, nlon), TOA→sfc order
        a, b: shape (nlev_hyb,), hybrid coefficients
        p0: scalar reference pressure (Pa)
        ps_2d: shape (nlat, nlon), surface pressure (Pa)
        plev_target_pa: shape (nlev_target,), target pressure levels (Pa)
                        — order doesn't matter, output matches input order

    Returns:
        field_plev: shape (nlev_target, nlat, nlon)

    Method: log-p linear interpolation per column. Above hybrid TOA: 0.
    Below hybrid bottom: extend bottom value to surface.

    NOTE: superseded by `hybrid_to_plev_mass_conserving` (now in
    `data/cmip6_common.py`) for actual case builds — kept here unmodified,
    untouched by the Phase 3 refactor, as a reference/diagnostic
    implementation (see scripts/diag_cloud_column.py).
    """
    nlev_hyb, nlat, nlon = field_hyb.shape
    nlev_target = len(plev_target_pa)
    out = np.zeros((nlev_target, nlat, nlon), dtype=np.float64)

    # log of target levels (stays constant per column)
    log_pt = np.log(plev_target_pa)

    for j in range(nlat):
        for i in range(nlon):
            # Compute hybrid layer pressures at this column
            p_hyb = a * p0 + b * ps_2d[j, i]   # shape (nlev_hyb,), TOA→sfc
            log_phyb = np.log(p_hyb)
            field_col = field_hyb[:, j, i]      # shape (nlev_hyb,)

            # np.interp requires increasing x (log_phyb is already increasing
            # since p_hyb goes TOA→sfc = small→large, log strictly increasing).
            # For target levels above hybrid top: extrapolate to 0.
            # For target levels below bottom: extrapolate to bottom value.
            interp_vals = np.interp(log_pt, log_phyb, field_col,
                                    left=0.0, right=field_col[-1])
            out[:, j, i] = interp_vals

    return out


def reorder_for_pycfram_input(arr_top_down):
    """If input is TOA→sfc (numpy convention used internally in this module),
    flip to sfc→TOA for pyCFRAM input NetCDF (which expects sfc→TOA, with
    `[::-1]` applied inside run_parallel_python.py to recover TOA→sfc).
    """
    return arr_top_down[::-1]
