"""Model-agnostic CMIP6 machinery shared by every per-model data source.

Extracted from `data/cesm2_cmip6_source.py` (Phase 3, WP-M4.1, see
`docs/plan_ph3.md` §2.2/§4) so that a DAMIP multi-model source
(`data/cmip6_damip_source.py`, WP-M4.2) can reuse the same numerics without
re-deriving CESM2's hardcodes for every new model. Nothing in this module
assumes a particular calendar, hybrid-coefficient naming, plev list, or
variable-availability pattern — those are exactly the axes of cross-model
heterogeneity documented in `docs/plan_ph3.md` §3.

`data/cesm2_cmip6_source.py` is now a thin CESM2-specific shim over this
module (noleap calendar, `a,b,p0` hybrid coefficients) — see that file's
module docstring for the regression-safety contract.
"""
import glob
import os
import re

import cftime
import numpy as np
from netCDF4 import Dataset


# mol/mol -> kg/kg mass-mixing-ratio conversion (M_O3 / M_air = 48/29).
# Lifted from scripts/inject_cesm_o3.py:40 (VMR_TO_MMR).
VMR_TO_MMR = 48.0 / 29.0

# Default TOA insolation constant, W/m^2 (docs/plan_ph3.md §6.2).
DEFAULT_SOLAR_CONSTANT = 1361.0

# CMIP6 fill-value threshold (typical fill ~1e20; anything with |x|>1e15 is
# unambiguously a fill/missing marker, never a physical field value).
FILLVALUE_ABS_THRESHOLD = 1e15


# ---------------------------------------------------------------------------
# Calendar-aware time decoding (replaces the noleap-only days/365.0 arithmetic
# formerly baked into cesm2_cmip6_source.years_to_month_indices()).
# ---------------------------------------------------------------------------

def _days_in_month(year, month, calendar):
    """Days in (year, month) under `calendar`, via cftime date arithmetic.
    Correct for noleap/365_day (Feb always 28), 360_day (every month 30),
    gregorian/proleptic_gregorian (real leap years)."""
    start = cftime.datetime(year, month, 1, calendar=calendar)
    if month == 12:
        end = cftime.datetime(year + 1, 1, 1, calendar=calendar)
    else:
        end = cftime.datetime(year, month + 1, 1, calendar=calendar)
    return float((end - start).days)


def decode_time(time_values, units, calendar):
    """Decode a CF time coordinate into per-record (year, month, day_weight).

    Calendar-aware via `cftime`, correctly handling noleap/365_day, 360_day,
    gregorian, proleptic_gregorian, and any other calendar cftime recognizes
    (docs/plan_ph3.md §3 "日历" row). This is the generalization of the
    CESM2-only `days/365.0` arithmetic in the pre-Phase-3
    `cesm2_cmip6_source.years_to_month_indices()`, which is silently wrong
    for every calendar except noleap/365_day.

    Parameters
    ----------
    time_values : array-like
        Raw numeric values of the NetCDF `time` variable.
    units : str
        CF `units` attribute, e.g. ``'days since 0001-01-01'``.
    calendar : str
        CF `calendar` attribute, e.g. ``'noleap'``, ``'360_day'``,
        ``'gregorian'``, ``'proleptic_gregorian'``, ``'365_day'``.

    Returns
    -------
    years : (n,) int64 ndarray
    months : (n,) int64 ndarray, 1-12
    day_weights : (n,) float64 ndarray — days in that record's calendar
        month (``360_day`` => always 30.0; noleap Feb => always 28.0; etc).
    """
    time_values = np.asarray(time_values, dtype=np.float64)
    dates = np.atleast_1d(cftime.num2date(time_values, units=units, calendar=calendar))
    years = np.array([d.year for d in dates], dtype=np.int64)
    months = np.array([d.month for d in dates], dtype=np.int64)
    day_weights = np.array(
        [_days_in_month(int(y), int(m), calendar) for y, m in zip(years, months)],
        dtype=np.float64)
    return years, months, day_weights


def years_to_month_indices(years, year_start, year_end):
    """Indices where ``year_start <= year <= year_end`` (inclusive), given a
    pre-decoded `years` array (see `decode_time`)."""
    years = np.asarray(years)
    mask = (years >= year_start) & (years <= year_end)
    return np.where(mask)[0]


def annual_climo_from_monthly(field, time_indices, day_weights):
    """Day-weighted annual climatology across `time_indices`, robust to
    CMIP6 fill-values (~1e20) that appear at plev cells below local surface
    pressure in some months but not others (seasonal ps variation).

    Fill-values are converted to NaN first, then a per-cell day-weighted mean
    is taken over only the valid selected months; cells with no valid data
    in ANY selected month remain NaN (caller/downstream `fill_subsurface`
    handles that).

    Parameters
    ----------
    field : netCDF4 Variable or ndarray, shape (ntime, ...).
    time_indices : (n,) int array — which records to include.
    day_weights : (n,) float array, ALREADY selected/aligned to
        `time_indices` (i.e. ``day_weights[i]`` weights ``time_indices[i]``,
        not the full unselected time axis).

    Returns
    -------
    ndarray, shape ``field.shape[1:]``, NaN where no valid data.
    """
    n = len(time_indices)
    if n == 0:
        raise ValueError('No time records selected')
    if len(day_weights) != n:
        raise ValueError('day_weights length (%d) must match time_indices (%d)' %
                          (len(day_weights), n))

    sub = np.asarray(field[time_indices], dtype=np.float64)
    sub = np.where(np.abs(sub) > FILLVALUE_ABS_THRESHOLD, np.nan, sub)

    w_b = np.asarray(day_weights, dtype=np.float64).reshape((n,) + (1,) * (sub.ndim - 1))

    valid = ~np.isnan(sub)
    num = np.nansum(np.where(valid, w_b * sub, 0.0), axis=0)
    den = np.sum(np.where(valid, w_b, 0.0), axis=0)
    with np.errstate(divide='ignore', invalid='ignore'):
        out = np.where(den > 0, num / den, np.nan)
    return out


# ---------------------------------------------------------------------------
# Hybrid sigma-pressure vertical coordinate: naming-convention detection +
# the (unmodified) mass-conserving hybrid->plev projection.
# ---------------------------------------------------------------------------

_FORMULA_TERM_RE = re.compile(r'(\w+):\s*(\S+)')


def _parse_formula_terms(formula_terms_str):
    return dict(_FORMULA_TERM_RE.findall(formula_terms_str))


def detect_vertical(formula_terms=None, override=None):
    """Determine which hybrid-sigma-pressure coefficient-naming convention a
    model uses, from the vertical coordinate's `formula_terms` CF attribute.

    - CESM2-style (``p = a*p0 + b*ps``): `formula_terms` contains ``a`` and
      ``p0`` -> ``scheme='hybrid_ab_p0'``.
    - CMOR-mainstream (``p = ap + b*ps``, e.g. IPSL/MRI/CNRM/HadGEM):
      `formula_terms` contains ``ap`` (no ``p0``) -> ``scheme='hybrid_ap_b'``.

    Parameters
    ----------
    formula_terms : str, optional
        The `formula_terms` attribute string of the vertical-coordinate
        (typically ``cl``) variable, e.g. ``"a: a b: b ps: ps p0: p0"`` or
        ``"ap: ap b: b ps: ps"``.
    override : dict, optional
        Explicit scheme dict (e.g. from a model's `configs/damip_models.d/
        <model>.yaml` ``vertical:`` block) that bypasses probing entirely and
        is returned verbatim.

    Returns
    -------
    dict with key ``'scheme'`` (``'hybrid_ab_p0'`` | ``'hybrid_ap_b'``) plus
    the source variable name for each coefficient term, e.g.
    ``{'scheme': 'hybrid_ab_p0', 'a': 'a', 'b': 'b', 'p0': 'p0', 'ps': 'ps'}``.
    """
    if override is not None:
        return dict(override)
    if not formula_terms:
        raise ValueError('detect_vertical: no formula_terms given and no override')
    terms = _parse_formula_terms(formula_terms)
    if 'p0' in terms and 'a' in terms:
        return {'scheme': 'hybrid_ab_p0', 'a': terms['a'], 'b': terms.get('b'),
                'p0': terms['p0'], 'ps': terms.get('ps')}
    if 'ap' in terms:
        return {'scheme': 'hybrid_ap_b', 'ap': terms['ap'], 'b': terms.get('b'),
                'ps': terms.get('ps')}
    raise ValueError('detect_vertical: unrecognized formula_terms %r' % formula_terms)


def normalize_hybrid_coeffs(scheme, a=None, b=None, p0=None, ap=None):
    """Normalize hybrid coefficients from either naming convention into
    ``(a_eff, b, p0_eff)`` such that layer pressure = ``a_eff*p0_eff + b*ps``
    — i.e. ready to feed directly into `hybrid_to_plev_mass_conserving`,
    which is NOT rewritten (docs/plan_ph3.md §4: "MOVE AS-IS"). For the CMOR
    mainstream ``ap,b`` convention this is exact since
    ``p = ap + b*ps = ap*1 + b*ps``.
    """
    if scheme == 'hybrid_ab_p0':
        if a is None or b is None or p0 is None:
            raise ValueError('hybrid_ab_p0 scheme requires a, b, p0')
        return np.asarray(a, dtype=np.float64), np.asarray(b, dtype=np.float64), float(p0)
    if scheme == 'hybrid_ap_b':
        if ap is None or b is None:
            raise ValueError('hybrid_ap_b scheme requires ap, b')
        return np.asarray(ap, dtype=np.float64), np.asarray(b, dtype=np.float64), 1.0
    raise ValueError('normalize_hybrid_coeffs: unknown scheme %r' % scheme)


def hybrid_to_plev_mass_conserving(field_hyb, a, b, p0, ps_2d, plev_target_pa):
    """Mass-conserving re-projection of mixing ratio from hybrid to plev.

    Builds cumulative mass M(p) = int_0^p f * dp on the hybrid grid
    (trapezoidal), interpolates M at target plev boundaries, and differences
    to extract layer-mean mixing ratios. By construction,
    sum(target_field * dp_target) over a column equals
    sum(hybrid_field * dp_hybrid) over the column — total mass conserved.

    Moved as-is from ``data/cesm2_cmip6_source.py`` (was L196) — this
    numerical kernel is model-agnostic (docs/plan_ph3.md §2.2/§4). For
    ``ap,b``-convention models, call `normalize_hybrid_coeffs('hybrid_ap_b',
    ap=..., b=...)` first to get ``(a_eff, b, p0_eff=1.0)`` and pass those
    here — the formula ``a*p0 + b*ps`` already covers both conventions.

    Args:
        field_hyb: (nlev_hyb, nlat, nlon) mixing ratio on hybrid (TOA->sfc)
        a, b: hybrid coefficients (TOA->sfc)
        p0: reference pressure (Pa)
        ps_2d: (nlat, nlon) surface pressure (Pa)
        plev_target_pa: target plev (Pa), TOA->sfc ordering (ascending)

    Returns:
        field_plev: (nlev_target, nlat, nlon) layer-mean mixing ratio,
                    TOA->sfc ordering. The value at index k represents the
                    mass-mean mixing ratio for the layer between
                    plev_target[k] (top, smaller p) and plev_target[k+1]
                    (bottom, larger p). For the LAST index, the layer is
                    between plev_target[-2] and ps (capped).
    """
    nlev_hyb = len(a)
    nlev_t = len(plev_target_pa)
    nlat, nlon = ps_2d.shape

    out = np.zeros((nlev_t, nlat, nlon), dtype=np.float64)

    for j in range(nlat):
        for i in range(nlon):
            ps = ps_2d[j, i]
            p_h = a * p0 + b * ps              # (nlev_hyb,) TOA->sfc, ascending
            f_h = field_hyb[:, j, i]
            f_h = np.where(np.isnan(f_h), 0.0, f_h)

            # Cumulative mass from TOA via trapezoidal. Anchors: M(p=0)=0,
            # M(p=ps) = M(p_h[-1]) + f_h[-1] * (ps - p_h[-1]).
            M_hyb = np.zeros(nlev_hyb)
            for kh in range(1, nlev_hyb):
                dp = p_h[kh] - p_h[kh-1]
                f_avg = 0.5 * (f_h[kh] + f_h[kh-1])
                M_hyb[kh] = M_hyb[kh-1] + f_avg * dp
            M_at_ps = M_hyb[-1] + f_h[-1] * (ps - p_h[-1])

            p_anchor = np.concatenate([[0.0], p_h, [ps]])
            M_anchor = np.concatenate([[0.0], M_hyb, [M_at_ps]])

            # Target layer k is between plev_target[k] (top) and plev_target[k+1]
            # (bottom). Last index: layer between plev_target[-1] and ps.
            for kt in range(nlev_t):
                p_top = plev_target_pa[kt]
                if kt < nlev_t - 1:
                    p_bot_nominal = plev_target_pa[kt + 1]
                else:
                    p_bot_nominal = ps    # last index: integrate to surface
                p_bot = min(p_bot_nominal, ps)
                if p_bot <= p_top:
                    out[kt, j, i] = 0.0   # subsurface or zero-thickness
                    continue
                M_top = np.interp(p_top, p_anchor, M_anchor)
                M_bot = np.interp(p_bot, p_anchor, M_anchor)
                out[kt, j, i] = (M_bot - M_top) / (p_bot - p_top)

    return out


def compute_albedo(rsus, rsds):
    """Surface albedo = rsus / rsds, clipped [0, 1]. Where rsds~=0 (polar
    night), set albedo=0 (no SW input -> albedo undefined, doesn't matter).

    Moved as-is from ``data/cesm2_cmip6_source.py`` (was L269).
    """
    with np.errstate(divide='ignore', invalid='ignore'):
        alb = np.where(rsds > 1.0, rsus / rsds, 0.0)
    return np.clip(alb, 0.0, 1.0)


# ---------------------------------------------------------------------------
# Arbitrary-plev re-interpolation (for ta/hus/o3 when a model's published
# plev list differs from the case's target grid).
# ---------------------------------------------------------------------------

def interp_plev_to_target(field_src, plev_src_pa, plev_target_pa):
    """Log-pressure linear interpolation of a pressure-level field onto an
    arbitrary target plev grid, along axis 0.

    Identity guarantee: when ``plev_target_pa == plev_src_pa`` (same values,
    any order), the output equals `field_src` exactly — `np.interp`
    evaluated exactly at a source abscissa returns that abscissa's ordinate
    with no additional rounding. Required by docs/plan_ph3.md §10 item 3.

    Parameters
    ----------
    field_src : ndarray, shape (nlev_src, ...) — any number of trailing dims.
    plev_src_pa, plev_target_pa : 1D arrays, Pa. Need not be sorted;
        internally sorted ascending before interpolation.

    Returns
    -------
    ndarray, shape ``(nlev_target,) + field_src.shape[1:]``.
    """
    field_src = np.asarray(field_src, dtype=np.float64)
    plev_src_pa = np.asarray(plev_src_pa, dtype=np.float64)
    plev_target_pa = np.asarray(plev_target_pa, dtype=np.float64)

    nlev_src = field_src.shape[0]
    trailing_shape = field_src.shape[1:]
    ncol = int(np.prod(trailing_shape)) if trailing_shape else 1
    flat = field_src.reshape(nlev_src, ncol)

    order = np.argsort(plev_src_pa)
    log_src = np.log(plev_src_pa[order])
    flat_sorted = flat[order, :]

    log_tgt = np.log(plev_target_pa)
    nlev_tgt = len(plev_target_pa)
    out = np.empty((nlev_tgt, ncol), dtype=np.float64)
    for c in range(ncol):
        out[:, c] = np.interp(log_tgt, log_src, flat_sorted[:, c])

    return out.reshape((nlev_tgt,) + trailing_shape)


# ---------------------------------------------------------------------------
# Grid normalization (lon 0-360 ascending, lat S->N ascending).
# ---------------------------------------------------------------------------

def normalize_grid(lat, lon, data=None):
    """Normalize `lon` to ascending [0, 360) and `lat` to ascending S->N.

    Idempotent: applying this twice in a row is a no-op on its own output
    (docs/plan_ph3.md §10 item 4).

    Parameters
    ----------
    lat, lon : 1D arrays (degrees).
    data : ndarray, optional. If given, its LAST TWO axes must be
        ``(..., lat, lon)``; they are reordered to match ``lat_out/lon_out``.

    Returns
    -------
    ``(lat_out, lon_out)`` if `data` is None, else
    ``(lat_out, lon_out, data_out)``.
    """
    lat = np.asarray(lat, dtype=np.float64)
    lon = np.asarray(lon, dtype=np.float64) % 360.0

    order_lat = np.argsort(lat, kind='stable')
    order_lon = np.argsort(lon, kind='stable')
    lat_out = lat[order_lat]
    lon_out = lon[order_lon]

    if data is None:
        return lat_out, lon_out

    data = np.asarray(data)
    data_out = np.take(data, order_lat, axis=-2)
    data_out = np.take(data_out, order_lon, axis=-1)
    return lat_out, lon_out, data_out


def regrid_horizontal_bilinear(lat_src, lon_src, field_src, lat_tgt, lon_tgt):
    """Bilinear regrid of a (..., nlat_src, nlon_src) field to (nlat_tgt,
    nlon_tgt), periodic in longitude.

    Some models publish a variable on a coarser native horizontal grid than
    the model's other Amon variables (confirmed against real downloaded
    data: MRI-ESM2-0's `o3` is on a 64x128 grid while ta/hus/ts/ps/cl are all
    on 160x320) -- CMIP6 has no requirement that every variable in a
    table_id share one grid. Any field whose own (lat, lon) doesn't match
    the case's reference grid (established from `ta`) must be regridded
    onto it before combining with other fields; otherwise a "shape
    mismatch" IndexError (or a wrong, unaligned broadcast if shapes happen
    to be compatible) awaits downstream in fill_subsurface/state assembly.

    Same bilinear-periodic-in-longitude approach as
    core/kernels.py's `_bilinear_regrid_2d` (independently implemented here,
    not imported -- data/ and core/ are separate layers and core/ must stay
    untouched per Phase 3's zero-core-changes acceptance criterion). Handles
    an arbitrary number of leading dimensions (e.g. (nlev, nlat, nlon)) by
    looping the 2D interpolation over them.

    Parameters
    ----------
    lat_src, lon_src : 1D arrays (degrees) -- field_src's native grid.
    field_src : (..., nlat_src, nlon_src) ndarray.
    lat_tgt, lon_tgt : 1D arrays (degrees) -- the reference grid to regrid onto.

    Returns
    -------
    (..., nlat_tgt, nlon_tgt) ndarray.
    """
    from scipy.interpolate import RegularGridInterpolator

    lat_src = np.asarray(lat_src, dtype=np.float64)
    lon_src = np.asarray(lon_src, dtype=np.float64)
    lat_tgt = np.asarray(lat_tgt, dtype=np.float64)
    lon_tgt = np.asarray(lon_tgt, dtype=np.float64)
    field_src = np.asarray(field_src, dtype=np.float64)

    if lat_src[0] > lat_src[-1]:
        lat_src = lat_src[::-1]
        field_src = field_src[..., ::-1, :]

    dlon = lon_src[1] - lon_src[0]
    lon_pad = np.concatenate(([lon_src[0] - dlon], lon_src, [lon_src[-1] + dlon]))
    lat_clip = np.clip(lat_tgt, lat_src.min(), lat_src.max())
    lon_wrapped = np.mod(lon_tgt - lon_src[0], 360.0) + lon_src[0]
    LA, LO = np.meshgrid(lat_clip, lon_wrapped, indexing='ij')
    pts = np.stack([LA.ravel(), LO.ravel()], axis=-1)

    leading_shape = field_src.shape[:-2]
    field_flat = field_src.reshape(-1, field_src.shape[-2], field_src.shape[-1])
    out_flat = np.empty((field_flat.shape[0], lat_tgt.size, lon_tgt.size))
    for i in range(field_flat.shape[0]):
        data_pad = np.concatenate(
            [field_flat[i, :, -1:], field_flat[i], field_flat[i, :, :1]], axis=1)
        interp = RegularGridInterpolator(
            (lat_src, lon_pad), data_pad, method='linear',
            bounds_error=False, fill_value=None)
        out_flat[i] = interp(pts).reshape(lat_tgt.size, lon_tgt.size)

    return out_flat.reshape(leading_shape + (lat_tgt.size, lon_tgt.size))


# ---------------------------------------------------------------------------
# Analytic (rsdt-free) TOA insolation fallback.
# ---------------------------------------------------------------------------

def analytic_solar(lat_deg, s0=DEFAULT_SOLAR_CONSTANT, ndays=365):
    """Annual-mean daily-mean TOA insolation as a function of latitude only
    (no rsdt data needed) — standard astronomical daily-mean-insolation
    formula. Used as the DAMIP ``solar=ANALYTIC`` fallback for models
    lacking ``rsdt`` (docs/plan_ph3.md §6.2: 5/13 hist-aer models). Since
    hist-aer freezes the solar constant, base and perturbed states get the
    SAME field from this function, so ``frc_solar == 0`` by construction.

    Method: for each day of year, solar declination
    ``delta = 23.45 deg * sin(2*pi*(284+day)/365)``; hour angle
    ``H0 = arccos(clip(-tan(phi)*tan(delta), -1, 1))`` (clipped for polar
    day/night); daily-mean flux
    ``Qday = (S0/pi) * (H0*sin(phi)*sin(delta) + cos(phi)*cos(delta)*sin(H0))``;
    annual mean = mean over `ndays`.

    Anchors (approximate, per docs/plan_ph3.md §6.2): area-weighted global
    mean ~= S0/4 ~= 340 W/m^2; equatorial annual mean ~= 417 W/m^2.

    Parameters
    ----------
    lat_deg : scalar or array-like, degrees.
    s0 : solar constant, W/m^2 (default 1361).
    ndays : int, number of days averaged over one year (default 365).

    Returns
    -------
    Same shape as `lat_deg` (python float if scalar input): W/m^2.
    """
    scalar_input = np.ndim(lat_deg) == 0
    lat_arr = np.atleast_1d(np.asarray(lat_deg, dtype=np.float64))
    phi = np.deg2rad(lat_arr)

    day = np.arange(1, ndays + 1, dtype=np.float64)
    delta = np.deg2rad(23.45 * np.sin(2.0 * np.pi * (284.0 + day) / 365.0))  # (ndays,)

    tan_phi = np.tan(phi)[:, None]        # (nlat, 1)
    tan_delta = np.tan(delta)[None, :]    # (1, ndays)
    x = np.clip(-tan_phi * tan_delta, -1.0, 1.0)
    H0 = np.arccos(x)                     # (nlat, ndays) hour angle, radians

    sin_phi = np.sin(phi)[:, None]
    cos_phi = np.cos(phi)[:, None]
    sin_delta = np.sin(delta)[None, :]
    cos_delta = np.cos(delta)[None, :]

    q_day = (s0 / np.pi) * (H0 * sin_phi * sin_delta + cos_phi * cos_delta * np.sin(H0))
    q_ann = q_day.mean(axis=1)

    if scalar_input:
        return float(q_ann[0])
    return q_ann.reshape(np.shape(lat_deg))


# ---------------------------------------------------------------------------
# O3 climatology injection (for models publishing no o3 at all).
# ---------------------------------------------------------------------------

def _interp_lev_lat(field_zm, lev_src_hpa, lat_src, lev_tgt_hpa, lat_tgt):
    """1D linear interpolation: first in log(p) at each source latitude,
    then in latitude at each target level. `field_zm`: (nlev_src, nlat_src)."""
    logp_src = np.log(lev_src_hpa)
    logp_tgt = np.log(lev_tgt_hpa)
    nlat_src = len(lat_src)

    o3_v = np.empty((len(lev_tgt_hpa), nlat_src))
    for j in range(nlat_src):
        o3_v[:, j] = np.interp(logp_tgt, logp_src, field_zm[:, j],
                               left=field_zm[0, j], right=field_zm[-1, j])

    out = np.empty((len(lev_tgt_hpa), len(lat_tgt)))
    for k in range(len(lev_tgt_hpa)):
        out[k, :] = np.interp(lat_tgt, lat_src, o3_v[k, :])
    return out


def o3_climatology(o3_source_path, lev_target_hpa, lat_target, lon_target,
                    var='O3'):
    """Load a CESM-style prescribed-O3 climatology NetCDF (mol/mol,
    ``(time, lev, lat, lon)``), reduce to an annual+zonal-mean profile,
    convert mol/mol -> kg/kg, and interpolate onto an arbitrary target grid.

    Lifted from ``scripts/inject_cesm_o3.py`` (``load_cesm_o3`` +
    ``interp_to_pycfram_grid``), generalized: the source path and target
    grid are now arguments instead of a hardcoded candidate-path list and
    the cesm2_4xco2_official case grid. Feeding the SAME output array to
    both base and perturbed states gives ``frc_o3 == 0`` by construction —
    the DAMIP ``o3=auto``/``o3=climatology`` branch (docs/plan_ph3.md §6.1).

    Parameters
    ----------
    o3_source_path : str — path to the O3 climatology NetCDF.
    lev_target_hpa : (nlev,) array, hPa.
    lat_target : (nlat,) array, degrees.
    lon_target : (nlon,) array, degrees (broadcast only; O3 is zonal-mean).
    var : str — name of the O3 variable in the source file (default 'O3').

    Returns
    -------
    o3_3d : (nlev, nlat, nlon) ndarray, kg/kg.
    """
    nc = Dataset(o3_source_path)
    o3_vmr = np.array(nc.variables[var][:])       # (time, lev, lat, lon) mol/mol
    lev_src = np.array(nc.variables['lev'][:])     # ascending pressure
    lat_src = np.array(nc.variables['lat'][:])
    nc.close()

    o3_zm_vmr = np.nanmean(o3_vmr, axis=(0, 3))    # (nlev_src, nlat_src)
    o3_zm_mmr = o3_zm_vmr * VMR_TO_MMR

    o3_zm = _interp_lev_lat(o3_zm_mmr, lev_src, lat_src, lev_target_hpa, lat_target)
    o3_3d = np.broadcast_to(
        o3_zm[:, :, None], (len(lev_target_hpa), len(lat_target), len(lon_target))
    ).copy()
    return o3_3d


# ---------------------------------------------------------------------------
# Subsurface (below-ground) layer filling.
# ---------------------------------------------------------------------------

DEFAULT_HOLD_VARS = ('q', 'o3')
DEFAULT_ZERO_VARS = ('cliq', 'cice', 'camt',
                      'bc', 'ocphi', 'ocpho', 'sulf', 'ss', 'dust')
DEFAULT_CLIP_NEG_VARS = ('cliq', 'cice')


def fill_subsurface(fields, lev_hpa, ps_hpa, ts,
                     hold_vars=DEFAULT_HOLD_VARS,
                     zero_vars=DEFAULT_ZERO_VARS,
                     clip_neg_vars=DEFAULT_CLIP_NEG_VARS):
    """Fill below-surface ("subsurface") pressure-level cells in-memory,
    BEFORE any NetCDF is written (`build_case_input.validate_states` rejects
    non-finite values at write time — docs/plan_ph3.md §2.1 — so this must
    happen inside a source's ``build_states()``, not as a build-after hook).

    Lifted from ``scripts/mask_subsurface_layers.py``'s ``mask_one()``,
    generalized to operate on in-memory arrays rather than a NetCDF file
    opened in ``r+`` mode.

    Strategy (unchanged from the original tool)
    --------------------------------------------
    - ``ta_lay``: subsurface layers set to `ts` (zero-thickness layer at
      ground temperature).
    - HOLD vars (default ``q``, ``o3``): copy the lowest-real-layer value
      down through the subsurface column. RRTMG cannot tolerate an
      exact-zero H2O/O3 layer, so HOLD avoids a hard zero-discontinuity
      while the column stays radiatively shadowed by the surface below it.
    - ZERO vars (default clouds + aerosol species): set to 0 in the
      subsurface (RRTMG tolerates exact-zero cloud/aerosol).
    - `clip_neg_vars`: additionally clip globally-negative values (numerical
      advection ghosts) to 0.

    A cell is also treated as subsurface if ``ta_lay`` is NaN or a raw
    fill-value there, regardless of the ``lev > ps`` test (CMIP6 plev
    products may mask extra cells beyond a simple ps comparison).

    Parameters
    ----------
    fields : dict[str, ndarray]
        Must include ``'ta_lay'`` (nlev, nlat, nlon); any of hold_vars /
        zero_vars present will be masked, absent ones silently skipped.
        Any other keys are passed through unchanged.
    lev_hpa : (nlev,) array, hPa — same vertical order as ``fields['ta_lay']``.
    ps_hpa : (nlat, nlon) array, hPa — surface pressure.
    ts : (nlat, nlon) array, K — surface skin temperature.

    Returns
    -------
    out_fields : dict[str, ndarray] — copies with subsurface filling applied
        (plus untouched pass-through of any other keys in `fields`).
    summary : dict — diagnostic counts/ranges, same spirit as
        ``mask_subsurface_layers.mask_one()``'s printed summary.
    """
    lev_hpa = np.asarray(lev_hpa, dtype=np.float64)
    ps_hpa = np.asarray(ps_hpa, dtype=np.float64)
    ts = np.asarray(ts, dtype=np.float64)
    ta_raw = np.asarray(fields['ta_lay'], dtype=np.float64)

    submask_ps = lev_hpa[:, None, None] > ps_hpa[None, :, :]
    submask_nan = np.isnan(ta_raw) | (np.abs(ta_raw) > FILLVALUE_ABS_THRESHOLD)
    submask_3d = submask_ps | submask_nan
    n_sub = int(submask_3d.sum())

    summary = {
        'n_subsurface': n_sub,
        'n_total': int(submask_3d.size),
        'n_nan_only': int((submask_nan & ~submask_ps).sum()),
    }

    out = {}

    ts_b = np.broadcast_to(ts[None, :, :], ta_raw.shape)
    ta_new = np.where(submask_3d, ts_b, ta_raw)
    out['ta_lay'] = ta_new
    if n_sub > 0:
        summary['ta_lay_old_range'] = (float(np.nanmin(ta_raw[submask_3d])),
                                        float(np.nanmax(ta_raw[submask_3d])))
    else:
        summary['ta_lay_old_range'] = (float('nan'), float('nan'))

    # "Lowest real layer" = largest lev among non-subsurface layers per column.
    real_mask_3d = ~submask_3d
    lev_3d = lev_hpa[:, None, None] * np.ones_like(submask_3d, dtype=np.float64)
    lev_masked = np.where(real_mask_3d, lev_3d, -np.inf)
    lowest_real_k = lev_masked.argmax(axis=0)

    for v in hold_vars:
        if v in fields:
            arr = np.asarray(fields[v], dtype=np.float64)
            held = np.take_along_axis(arr, lowest_real_k[None, :, :], axis=0)[0]
            held_b = np.broadcast_to(held[None, :, :], arr.shape)
            out[v] = np.where(submask_3d, held_b, arr)
            summary[v + '_hold_mean'] = float(held.mean())

    for v in zero_vars:
        if v in fields:
            arr = np.asarray(fields[v], dtype=np.float64)
            old_max = float(arr[submask_3d].max()) if n_sub > 0 else 0.0
            arr_new = np.where(submask_3d, 0.0, arr)
            if v in clip_neg_vars:
                neg_count = int((arr_new < 0).sum())
                if neg_count > 0:
                    arr_new = np.where(arr_new < 0, 0.0, arr_new)
                summary[v + '_neg_clipped'] = neg_count
            out[v] = arr_new
            summary[v + '_sub_max'] = old_max

    for k, v in fields.items():
        if k not in out:
            out[k] = v

    return out, summary


# ---------------------------------------------------------------------------
# Generic CMIP6 filename discovery.
# ---------------------------------------------------------------------------

def discover_files(raw_dir, variables, model=None, experiment=None,
                    variant=None, grid=None, table_id='Amon'):
    """Glob-based CMIP6 filename discovery.

    Pattern: ``<var>_<table_id>_<model>_<experiment>_<variant>_<grid>_
    <time-range>.nc`` (standard CMOR output filename convention). Any of
    model/experiment/variant/grid left as None is wildcarded.

    Returns
    -------
    dict[str, list[str]]: variable -> sorted list of matching file paths,
    present only for variables that had >= 1 match. A variable simply being
    absent from the returned dict is NOT an error here — "not published for
    this model/experiment" is DAMIP's main path, not an edge case
    (docs/plan_ph3.md §3) — callers decide what to do (skip/fallback).
    """
    out = {}
    for var in variables:
        pattern = '%s_%s_%s_%s_%s_%s_*.nc' % (
            var, table_id, model or '*', experiment or '*',
            variant or '*', grid or '*')
        matches = sorted(glob.glob(os.path.join(raw_dir, pattern)))
        if matches:
            out[var] = matches
    return out


def discover_variant(raw_dir, model, experiment, var='ta', table_id='Amon',
                      grid=None):
    """List variant_labels (e.g. ``r1i1p1f1``) available for
    ``(model, experiment)`` by globbing on a single reference variable.

    Returns sorted list of distinct variant tokens found — the
    dictionary-order first entry is the conventional glob-fallback default
    (docs/plan_ph3.md §3 "variant_label" row).
    """
    pattern = '%s_%s_%s_%s_*_%s_*.nc' % (var, table_id, model, experiment, grid or '*')
    matches = sorted(glob.glob(os.path.join(raw_dir, pattern)))
    variants = set()
    for m in matches:
        parts = os.path.basename(m).split('_')
        # <var>_<table>_<model>_<experiment>_<variant>_<grid>_<trange>.nc
        if len(parts) >= 5:
            variants.add(parts[4])
    return sorted(variants)
