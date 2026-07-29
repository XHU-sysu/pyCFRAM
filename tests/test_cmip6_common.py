"""Unit tests for data/cmip6_common.py — the model-agnostic CMIP6 machinery
extracted in Phase 3 WP-M4.1 (docs/plan_ph3.md §10).

Covers, per the plan's explicit test list:
  1. decode_time: noleap, 360_day, gregorian (+ proleptic_gregorian) calendars.
  2. hybrid_to_plev_mass_conserving: column mass conservation, BOTH
     coefficient conventions (a,b,p0 via normalize_hybrid_coeffs, and
     ap,b directly).
  3. interp_plev_to_target: identity when target == source plev.
  4. normalize_grid: idempotency (lon wrap, lat flip).
  5. analytic_solar: two anchors (~340 W/m^2 global mean, ~417 W/m^2 equator).

Plus lighter coverage of the remaining new functions (detect_vertical,
o3_climatology, fill_subsurface, discover_files/discover_variant) so nothing
in the new module is left completely untested.
"""
import os
import sys

import cftime
import numpy as np
import pytest
from netCDF4 import Dataset

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from data import cmip6_common as common


# ---------------------------------------------------------------------------
# 1. decode_time — multiple calendars
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('calendar,year,month,day,expected_days_in_month', [
    ('noleap', 2001, 2, 15, 28.0),      # noleap Feb always 28
    ('noleap', 2001, 4, 15, 30.0),
    ('365_day', 2000, 2, 15, 28.0),     # alias of noleap, still 28 (no leap rule)
    ('360_day', 2001, 2, 15, 30.0),     # 360_day: every month exactly 30 days
    ('360_day', 2001, 7, 15, 30.0),
    ('gregorian', 2000, 2, 15, 29.0),   # 2000 is a leap year under gregorian
    ('gregorian', 2001, 2, 15, 28.0),   # 2001 is not
    ('proleptic_gregorian', 2004, 2, 15, 29.0),  # 2004 leap
    ('proleptic_gregorian', 2003, 2, 15, 28.0),
])
def test_decode_time_calendars(calendar, year, month, day, expected_days_in_month):
    units = 'days since 0001-01-01'
    t = cftime.date2num(cftime.datetime(year, month, day, calendar=calendar),
                        units=units, calendar=calendar)
    years, months, day_weights = common.decode_time(np.array([t]), units, calendar)
    assert years[0] == year
    assert months[0] == month
    assert day_weights[0] == expected_days_in_month


def test_decode_time_sequential_noleap_year_selection():
    """36 consecutive noleap months starting Jan year 1 -> selecting years
    [2,2] must return exactly the 12 indices of the second year, and every
    day_weight must match the standard non-leap month-length table."""
    calendar = 'noleap'
    units = 'days since 0001-01-01'
    dates = []
    y, m = 1, 1
    for _ in range(36):
        dates.append(cftime.datetime(y, m, 15, calendar=calendar))
        m += 1
        if m > 12:
            m = 1
            y += 1
    t = cftime.date2num(dates, units=units, calendar=calendar)
    years, months, day_weights = common.decode_time(np.asarray(t), units, calendar)

    idx = common.years_to_month_indices(years, 2, 2)
    assert list(idx) == list(range(12, 24))

    noleap_table = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    for i in idx:
        assert day_weights[i] == noleap_table[months[i] - 1]


def test_annual_climo_from_monthly_fillvalue_and_weighting():
    """Fill-values (>1e15) are masked to NaN; day-weighted mean matches a
    hand-computed reference; an all-fill cell stays NaN."""
    # 3 selected months with distinct day_weights and one fill-value month
    # for one cell.
    field = np.array([
        [[10.0, 1.0e20]],
        [[20.0, 5.0]],
        [[30.0, 5.0]],
    ])  # shape (3, 1, 2)
    time_indices = np.array([0, 1, 2])
    day_weights = np.array([31.0, 28.0, 31.0])

    out = common.annual_climo_from_monthly(field, time_indices, day_weights)
    assert out.shape == (1, 2)
    # cell 0: (10*31 + 20*28 + 30*31) / (31+28+31)
    expected0 = (10.0 * 31 + 20.0 * 28 + 30.0 * 31) / (31 + 28 + 31)
    assert np.isclose(out[0, 0], expected0)
    # cell 1: month 0 is fill -> excluded; (5*28 + 5*31)/(28+31)
    expected1 = (5.0 * 28 + 5.0 * 31) / (28 + 31)
    assert np.isclose(out[0, 1], expected1)


def test_annual_climo_from_monthly_all_fill_stays_nan():
    field = np.array([[[1.0e20]], [[1.0e20]]])
    out = common.annual_climo_from_monthly(field, np.array([0, 1]), np.array([30.0, 31.0]))
    assert np.isnan(out[0, 0])


# ---------------------------------------------------------------------------
# 2. hybrid_to_plev_mass_conserving — mass conservation, both conventions
# ---------------------------------------------------------------------------

def _column_mass(field_hyb_col, p_edges, ps):
    """Reference column-integrated mass (TOA to actual surface pressure
    `ps`) via simple trapezoidal rule on the hybrid layers themselves, plus
    the same "extend last hybrid value down to ps" convention the function
    under test uses — independent re-implementation for the test, not
    reusing the function's internals."""
    total = 0.0
    for k in range(1, len(p_edges)):
        dp = p_edges[k] - p_edges[k - 1]
        total += 0.5 * (field_hyb_col[k] + field_hyb_col[k - 1]) * dp
    total += field_hyb_col[-1] * (ps - p_edges[-1])
    return total


def test_hybrid_to_plev_mass_conserving_ab_p0_convention():
    a = np.array([0.04, 0.03, 0.015, 0.005, 0.0])
    b = np.array([0.0, 0.02, 0.15, 0.55, 0.95])
    p0 = 100000.0
    ps_2d = np.array([[100000.0]])
    field_hyb = np.array([1.0e-5, 2.0e-5, 3.0e-5, 1.5e-5, 0.5e-5]).reshape(5, 1, 1)

    plev_target = np.array([2000.0, 10000.0, 40000.0, 70000.0, 100000.0])  # TOA->sfc

    out = common.hybrid_to_plev_mass_conserving(field_hyb, a, b, p0, ps_2d, plev_target)
    assert out.shape == (5, 1, 1)

    p_hyb = a * p0 + b * ps_2d[0, 0]
    mass_hyb = _column_mass(field_hyb[:, 0, 0], p_hyb, ps_2d[0, 0])

    # Reconstruct mass on the target grid from the mass-conserving layer
    # means: sum(out[k] * dp[k]) where dp[k] spans [plev_target[k],
    # plev_target[k+1]] (last layer extends to ps).
    dp = np.diff(np.concatenate([plev_target, [ps_2d[0, 0]]]))
    mass_tgt = np.sum(out[:, 0, 0] * dp)
    assert np.isclose(mass_tgt, mass_hyb, rtol=1e-6)


def test_hybrid_to_plev_mass_conserving_ap_b_convention_matches_ab_p0():
    """normalize_hybrid_coeffs('hybrid_ap_b', ...) must feed
    hybrid_to_plev_mass_conserving identically to the equivalent a,b,p0
    call — verifying the function (moved as-is, not rewritten) genuinely
    supports both CMIP6 hybrid-coefficient naming conventions."""
    a = np.array([0.04, 0.03, 0.015, 0.005, 0.0])
    b = np.array([0.0, 0.02, 0.15, 0.55, 0.95])
    p0 = 100000.0
    ps_2d = np.array([[100000.0], [95000.0]])
    field_hyb = np.array([1.0e-5, 2.0e-5, 3.0e-5, 1.5e-5, 0.5e-5])
    field_hyb = np.broadcast_to(field_hyb[:, None, None], (5, 2, 1)).copy()

    plev_target = np.array([2000.0, 10000.0, 40000.0, 70000.0, 100000.0])

    out_ab_p0 = common.hybrid_to_plev_mass_conserving(field_hyb, a, b, p0, ps_2d, plev_target)

    # ap,b convention: ap = a*p0 (physically identical hybrid levels),
    # detected/normalized independently of the ab_p0 path.
    ap = a * p0
    a_eff, b_eff, p0_eff = common.normalize_hybrid_coeffs('hybrid_ap_b', ap=ap, b=b)
    assert p0_eff == 1.0
    out_ap_b = common.hybrid_to_plev_mass_conserving(field_hyb, a_eff, b_eff, p0_eff, ps_2d, plev_target)

    np.testing.assert_allclose(out_ab_p0, out_ap_b)


def test_normalize_hybrid_coeffs_ab_p0_roundtrip():
    a = np.array([0.1, 0.2])
    b = np.array([0.3, 0.4])
    p0 = 100000.0
    a_eff, b_eff, p0_eff = common.normalize_hybrid_coeffs('hybrid_ab_p0', a=a, b=b, p0=p0)
    np.testing.assert_array_equal(a_eff, a)
    np.testing.assert_array_equal(b_eff, b)
    assert p0_eff == p0


def test_normalize_hybrid_coeffs_unknown_scheme_raises():
    with pytest.raises(ValueError):
        common.normalize_hybrid_coeffs('bogus_scheme')


def test_detect_vertical_ab_p0():
    info = common.detect_vertical(formula_terms='p0: p0 a: a b: b ps: ps')
    assert info['scheme'] == 'hybrid_ab_p0'
    assert info['a'] == 'a' and info['b'] == 'b' and info['p0'] == 'p0'


def test_detect_vertical_ap_b():
    info = common.detect_vertical(formula_terms='ap: ap b: b ps: ps')
    assert info['scheme'] == 'hybrid_ap_b'
    assert info['ap'] == 'ap' and info['b'] == 'b'


def test_detect_vertical_override_bypasses_probing():
    override = {'scheme': 'hybrid_ap_b', 'ap': 'my_ap', 'b': 'my_b'}
    assert common.detect_vertical(override=override) == override


def test_detect_vertical_unrecognized_raises():
    with pytest.raises(ValueError):
        common.detect_vertical(formula_terms='something: weird')
    with pytest.raises(ValueError):
        common.detect_vertical()


def test_compute_albedo_polar_night_safe():
    rsus = np.array([50.0, 0.0])
    rsds = np.array([100.0, 0.0])
    alb = common.compute_albedo(rsus, rsds)
    assert np.isclose(alb[0], 0.5)
    assert alb[1] == 0.0  # rsds~=0 -> albedo undefined -> 0, not NaN/inf


# ---------------------------------------------------------------------------
# 3. interp_plev_to_target — identity when target == source
# ---------------------------------------------------------------------------

def test_interp_plev_to_target_identity_same_order():
    plev = np.array([100000.0, 70000.0, 50000.0, 30000.0, 10000.0, 3000.0])
    field = np.random.RandomState(0).uniform(200, 300, size=(6, 3, 4))
    out = common.interp_plev_to_target(field, plev, plev)
    np.testing.assert_array_equal(out, field)


def test_interp_plev_to_target_identity_reordered_target():
    """Identity must hold even if plev_target is given in a different
    (but value-identical) order than plev_src — the function sorts
    internally."""
    plev_src = np.array([3000.0, 10000.0, 30000.0, 50000.0, 70000.0, 100000.0])
    field = np.random.RandomState(1).uniform(200, 300, size=(6, 2, 2))
    plev_target = plev_src[::-1].copy()  # same values, reversed order

    out = common.interp_plev_to_target(field, plev_src, plev_target)
    # out is ordered according to plev_target; reversing it back should
    # equal the original field ordered by plev_src.
    np.testing.assert_allclose(out[::-1], field)


def test_interp_plev_to_target_monotone_between_knots():
    """Basic correctness (not just the identity case): interpolated value
    at a level between two source levels should lie between their values
    for a monotone field (log-p linear interpolation preserves monotonicity
    for a linear-in-log-p field by construction)."""
    plev_src = np.array([100000.0, 10000.0, 1000.0])
    field = np.array([300.0, 220.0, 200.0]).reshape(3, 1, 1)
    plev_target = np.array([50000.0])
    out = common.interp_plev_to_target(field, plev_src, plev_target)
    assert 220.0 < out[0, 0, 0] < 300.0


# ---------------------------------------------------------------------------
# 4. normalize_grid — idempotency
# ---------------------------------------------------------------------------

def test_normalize_grid_lon_wrap_and_lat_flip():
    lat = np.array([60.0, -60.0, 0.0])        # not S->N ascending
    lon = np.array([-170.0, 200.0, 10.0])     # mixed -180..180 / 0..360, no aliasing ties
    data = np.zeros((3, 3))
    for j in range(3):
        for i in range(3):
            data[j, i] = j * 10 + i          # tag by original (lat_idx, lon_idx)

    lat_out, lon_out, data_out = common.normalize_grid(lat, lon, data)

    assert np.all(np.diff(lat_out) > 0)
    assert np.all(np.diff(lon_out) > 0)
    assert np.all((lon_out >= 0) & (lon_out < 360))
    # Data was permuted consistently: verify by re-deriving expected order.
    order_lat = np.argsort(lat)
    order_lon = np.argsort(lon % 360.0)
    expected = data[np.ix_(order_lat, order_lon)]
    np.testing.assert_array_equal(data_out, expected)


def test_normalize_grid_idempotent():
    lat = np.array([60.0, -60.0, 0.0, 30.0])
    lon = np.array([-90.0, 270.1, 10.0, 359.0])
    data = np.arange(16, dtype=np.float64).reshape(4, 4)

    lat1, lon1, data1 = common.normalize_grid(lat, lon, data)
    lat2, lon2, data2 = common.normalize_grid(lat1, lon1, data1)

    np.testing.assert_array_equal(lat1, lat2)
    np.testing.assert_array_equal(lon1, lon2)
    np.testing.assert_array_equal(data1, data2)


def test_normalize_grid_no_data_arg():
    lat_out, lon_out = common.normalize_grid(np.array([10.0, -10.0]), np.array([200.0, -20.0]))
    assert lat_out[0] < lat_out[1]
    assert np.all((lon_out >= 0) & (lon_out < 360))


# ---------------------------------------------------------------------------
# 4a2. fill_nan_hold_toward_surface — CMIP6 below-ground-plev masking guard
# (real bug: MRI-ESM2-0's o3 has ~23% NaN at 1000 hPa, tapering to 0% by
# 600 hPa -- below-ground terrain masking -- which bilinear regridding
# would otherwise spread into valid neighboring cells)
# ---------------------------------------------------------------------------

def test_fill_nan_hold_toward_surface_fills_surface_end_run():
    # (nlev=4, nlat=1, nlon=2), sfc->TOA. Column 0 has NaN at levels 0-1
    # (below-ground at that terrain point); column 1 has no missing data.
    field = np.array([
        [[np.nan, 10.0]],
        [[np.nan, 11.0]],
        [[5.0, 12.0]],
        [[6.0, 13.0]],
    ])
    out = common.fill_nan_hold_toward_surface(field)
    assert not np.any(np.isnan(out))
    # Column 0's masked levels held at the shallowest valid level's value (5.0).
    np.testing.assert_array_equal(out[:, 0, 0], [5.0, 5.0, 5.0, 6.0])
    # Column 1 (no NaN) passes through unchanged.
    np.testing.assert_array_equal(out[:, 0, 1], [10.0, 11.0, 12.0, 13.0])


def test_fill_nan_hold_toward_surface_leaves_all_nan_column_untouched():
    field = np.full((3, 1, 1), np.nan)
    out = common.fill_nan_hold_toward_surface(field)
    assert np.all(np.isnan(out))  # not silently zero-filled or otherwise papered over


def test_fill_nan_hold_toward_surface_no_nan_is_noop():
    field = np.arange(24.0).reshape(4, 2, 3)
    out = common.fill_nan_hold_toward_surface(field)
    np.testing.assert_array_equal(out, field)


# ---------------------------------------------------------------------------
# 4b. regrid_horizontal_bilinear — cross-model horizontal grid mismatch
# (real bug: MRI-ESM2-0's o3 is published on a 64x128 grid while every
# other Amon variable in the same model is 160x320)
# ---------------------------------------------------------------------------

def test_regrid_horizontal_bilinear_identity_when_grids_match():
    """Regridding onto the exact same grid must reproduce the input closely
    (bilinear identity, modulo tiny floating-point interpolation noise)."""
    lat_src = np.linspace(-90, 90, 8)
    lon_src = np.linspace(0, 360, 16, endpoint=False)
    LA, LO = np.meshgrid(lat_src, lon_src, indexing='ij')
    field = np.sin(np.radians(LA)) * np.cos(np.radians(LO))
    out = common.regrid_horizontal_bilinear(lat_src, lon_src, field, lat_src, lon_src)
    np.testing.assert_allclose(out, field, atol=1e-8)


def test_regrid_horizontal_bilinear_coarse_to_fine_shape_and_range():
    """Coarse (64x128-like) -> fine (160x320-like) grid: output shape must
    match the target, and values must stay within the source field's range
    (bilinear doesn't overshoot for a smooth field)."""
    lat_src = np.linspace(-90, 90, 8)
    lon_src = np.linspace(0, 360, 16, endpoint=False)
    LA, LO = np.meshgrid(lat_src, lon_src, indexing='ij')
    field = 10.0 + 5.0 * np.sin(np.radians(LA))  # smooth, no lon dependence for a simple bound check

    lat_tgt = np.linspace(-90, 90, 20)
    lon_tgt = np.linspace(0, 360, 40, endpoint=False)
    out = common.regrid_horizontal_bilinear(lat_src, lon_src, field, lat_tgt, lon_tgt)

    assert out.shape == (20, 40)
    assert out.min() >= field.min() - 1e-6
    assert out.max() <= field.max() + 1e-6


def test_regrid_horizontal_bilinear_handles_leading_level_dimension():
    """A (nlev, nlat, nlon) field regrids level-by-level, preserving nlev."""
    lat_src = np.linspace(-90, 90, 8)
    lon_src = np.linspace(0, 360, 16, endpoint=False)
    field = np.random.RandomState(0).rand(5, 8, 16)  # (nlev=5, nlat, nlon)

    lat_tgt = np.linspace(-90, 90, 20)
    lon_tgt = np.linspace(0, 360, 40, endpoint=False)
    out = common.regrid_horizontal_bilinear(lat_src, lon_src, field, lat_tgt, lon_tgt)
    assert out.shape == (5, 20, 40)
    assert np.all(np.isfinite(out))


def test_regrid_horizontal_bilinear_longitude_periodic_wrap():
    """A field with a feature straddling the lon=0/360 seam must interpolate
    smoothly across it (periodic padding), not show a discontinuity."""
    lat_src = np.linspace(-90, 90, 8)
    lon_src = np.linspace(0, 360, 16, endpoint=False)
    LA, LO = np.meshgrid(lat_src, lon_src, indexing='ij')
    field = np.cos(np.radians(LO))  # smooth periodic function of lon only

    # Target point right at the seam (lon=359) should be close to cos(359deg),
    # not jump to some unrelated value from bad edge handling.
    lat_tgt = np.array([0.0])
    lon_tgt = np.array([359.0])
    out = common.regrid_horizontal_bilinear(lat_src, lon_src, field, lat_tgt, lon_tgt)
    expected = np.cos(np.radians(359.0))
    assert abs(out[0, 0] - expected) < 0.05


# ---------------------------------------------------------------------------
# 5. analytic_solar — two anchors
# ---------------------------------------------------------------------------

def test_analytic_solar_equatorial_anchor():
    q_eq = common.analytic_solar(0.0)
    assert 380.0 < q_eq < 450.0  # plan anchor ~417 W/m^2, generous band


def test_analytic_solar_global_area_weighted_mean_anchor():
    """Area-weighted (cos-latitude) global mean should be close to S0/4."""
    lat = np.linspace(-89.5, 89.5, 180)
    q = common.analytic_solar(lat)
    weights = np.cos(np.deg2rad(lat))
    global_mean = np.sum(q * weights) / np.sum(weights)
    assert abs(global_mean - 1361.0 / 4.0) < 15.0  # within ~15 W/m^2 of 340.25


def test_analytic_solar_scalar_vs_array_shape():
    assert isinstance(common.analytic_solar(30.0), float)
    arr = common.analytic_solar(np.array([0.0, 30.0, 60.0]))
    assert arr.shape == (3,)


def test_analytic_solar_symmetric_about_equator():
    """No solstice/hemisphere asymmetry baked in (annual mean) -> +lat and
    -lat give identical insolation."""
    q_pos = common.analytic_solar(45.0)
    q_neg = common.analytic_solar(-45.0)
    assert np.isclose(q_pos, q_neg)


# ---------------------------------------------------------------------------
# o3_climatology
# ---------------------------------------------------------------------------

@pytest.fixture
def synthetic_o3_nc(tmp_path):
    path = str(tmp_path / 'ozone_synthetic.nc')
    nlev, nlat, nlon, ntime = 4, 5, 6, 2
    lev = np.array([1000.0, 100.0, 10.0, 1.0])   # ascending pressure (hPa), as in real file
    lat = np.linspace(-90, 90, nlat)
    with Dataset(path, 'w') as nc:
        nc.createDimension('time', ntime)
        nc.createDimension('lev', nlev)
        nc.createDimension('lat', nlat)
        nc.createDimension('lon', nlon)
        nc.createVariable('lev', 'f8', ('lev',))[:] = lev
        nc.createVariable('lat', 'f8', ('lat',))[:] = lat
        o3 = np.empty((ntime, nlev, nlat, nlon))
        for k in range(nlev):
            o3[:, k, :, :] = 1e-6 * (k + 1)   # constant per level, mol/mol
        nc.createVariable('O3', 'f8', ('time', 'lev', 'lat', 'lon'))[:] = o3
    return path, lev, lat


def test_o3_climatology_unit_conversion_and_shape(synthetic_o3_nc):
    path, lev_src, lat_src = synthetic_o3_nc
    lev_tgt = np.array([900.0, 50.0, 2.0])
    lat_tgt = np.linspace(-80, 80, 4)
    lon_tgt = np.array([0.0, 120.0, 240.0])

    o3_3d = common.o3_climatology(path, lev_tgt, lat_tgt, lon_tgt)
    assert o3_3d.shape == (3, 4, 3)
    # mol/mol -> kg/kg: should be scaled up by 48/29 relative to source vmr
    # order of magnitude check (source ~1e-6 mol/mol level-1 -> ~1.65e-6 kg/kg)
    assert 1e-7 < o3_3d.min() < 1e-4
    # broadcast over lon: identical across the last axis
    np.testing.assert_array_equal(o3_3d[:, :, 0], o3_3d[:, :, 1])


# ---------------------------------------------------------------------------
# fill_subsurface
# ---------------------------------------------------------------------------

def test_fill_subsurface_hold_and_zero_strategies():
    lev_hpa = np.array([10.0, 500.0, 1000.0])   # TOA->sfc
    ps_hpa = np.array([[950.0]])                # column: last layer (1000) is subsurface
    ts = np.array([[280.0]])

    fields = {
        'ta_lay': np.array([220.0, 260.0, 999.0]).reshape(3, 1, 1),
        'q': np.array([1e-6, 5e-3, 0.0]).reshape(3, 1, 1),
        'cliq': np.array([0.0, 1e-5, 1e-5]).reshape(3, 1, 1),
        'untouched_var': np.array([1.0, 2.0, 3.0]).reshape(3, 1, 1),
    }

    out, summary = common.fill_subsurface(fields, lev_hpa, ps_hpa, ts)

    assert summary['n_subsurface'] == 1
    # ta_lay subsurface -> ts
    assert out['ta_lay'][2, 0, 0] == 280.0
    assert out['ta_lay'][1, 0, 0] == 260.0  # untouched real layer
    # q HOLD: subsurface copies lowest-real-layer (index 1) value
    assert out['q'][2, 0, 0] == 5e-3
    # cliq ZERO strategy
    assert out['cliq'][2, 0, 0] == 0.0
    # pass-through
    np.testing.assert_array_equal(out['untouched_var'], fields['untouched_var'])


def test_fill_subsurface_clips_negative_cliq():
    lev_hpa = np.array([500.0, 1000.0])
    ps_hpa = np.array([[1000.0]])   # no subsurface layers
    ts = np.array([[280.0]])
    fields = {
        'ta_lay': np.array([260.0, 270.0]).reshape(2, 1, 1),
        'cliq': np.array([-1e-6, 1e-5]).reshape(2, 1, 1),
    }
    out, summary = common.fill_subsurface(fields, lev_hpa, ps_hpa, ts)
    assert summary['n_subsurface'] == 0
    assert out['cliq'][0, 0, 0] == 0.0  # negative clipped even though not subsurface
    assert summary['cliq_neg_clipped'] == 1


def test_fill_subsurface_nan_ta_treated_as_subsurface():
    lev_hpa = np.array([500.0, 1000.0])
    ps_hpa = np.array([[1000.0]])   # ps says both layers are "real"
    ts = np.array([[280.0]])
    fields = {'ta_lay': np.array([260.0, np.nan]).reshape(2, 1, 1),
              'q': np.array([1e-3, 1e-3]).reshape(2, 1, 1)}
    out, summary = common.fill_subsurface(fields, lev_hpa, ps_hpa, ts)
    assert summary['n_subsurface'] == 1
    assert summary['n_nan_only'] == 1
    assert out['ta_lay'][1, 0, 0] == 280.0


# ---------------------------------------------------------------------------
# discover_files / discover_variant
# ---------------------------------------------------------------------------

def test_discover_files_and_variant(tmp_path):
    raw_dir = str(tmp_path)
    names = [
        'ta_Amon_IPSL-CM6A-LR_hist-aer_r1i1p1f1_gr_185001-202012.nc',
        'hus_Amon_IPSL-CM6A-LR_hist-aer_r1i1p1f1_gr_185001-202012.nc',
        'ta_Amon_IPSL-CM6A-LR_hist-aer_r2i1p1f1_gr_185001-202012.nc',
    ]
    for n in names:
        open(os.path.join(raw_dir, n), 'w').close()

    files = common.discover_files(raw_dir, ['ta', 'hus', 'cl'],
                                   model='IPSL-CM6A-LR', experiment='hist-aer')
    assert set(files.keys()) == {'ta', 'hus'}   # 'cl' absent -> not in dict
    assert len(files['ta']) == 2
    assert len(files['hus']) == 1

    variants = common.discover_variant(raw_dir, 'IPSL-CM6A-LR', 'hist-aer', var='ta')
    assert variants == ['r1i1p1f1', 'r2i1p1f1']
