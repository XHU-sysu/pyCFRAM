"""CMIP6 DAMIP single-forcing data source (Phase 3, WP-M4.2).

Implements ``@register_source('cmip6_damip')`` per the ``data/source_base.py``
plugin contract used by ``data/era5_source.py`` and (as of WP-M4.1)
``data/cesm2_cmip6_source.py``: ``build_states()`` returns
``(base_state, perturbed_state, nonrad_forcing)`` and is consumed by the
generic writer ``scripts/build_case_input.py``.

This module is model-agnostic: all cross-model heterogeneity (calendar,
hybrid-coefficient naming, published plev list, variant/grid labels,
variable availability) is either auto-detected at runtime or supplied via
``configs/damip_models.d/<model>.yaml`` -- see ``docs/plan_ph3.md`` §3/§4.
The actual numerics (calendar-aware time decoding, mass-conserving
hybrid->plev projection, log-p interpolation, grid normalization, analytic
solar fallback, O3 climatology injection, subsurface filling) all live in
``data/cmip6_common.py`` (WP-M4.1) and are reused here unchanged.

The ten-step ``build_states()`` flow below follows ``docs/plan_ph3.md`` §5.1
step-by-step (each step is commented with its plan reference). §6 is the
"missing-variable decision tree": cloud (ACTIVE/SKIPPED), O3
(MODEL/CLIMATOLOGY/SKIPPED), solar (ACTIVE/ANALYTIC), aerosol (SKIPPED),
flux (ACTIVE/UNAVAILABLE), albedo (ACTIVE/FALLBACK). Note: §6.1's prose uses
per-family tags (e.g. "o3=MODEL") that are more specific than the six-tag
roll-up vocabulary given in §6.3 (ACTIVE/SKIPPED/CLIMATOLOGY/ANALYTIC/
FALLBACK/UNAVAILABLE) -- this module's ``self.provenance`` records BOTH: the
precise per-family tag from §6.1 as ``status``, plus that same value already
falls inside the §6.3 set for every family except O3 (which uses "MODEL"
where a strict six-tag reading would say "ACTIVE" -- §8 WP-M4.4's own worked
example explicitly writes "o3=MODEL", so this module follows the plan's more
specific, worked-example wording rather than the terser summary list). This
is a judgement call flagged explicitly in this WP's hand-off report so the
downstream summary-writer WP (M4.4, not part of this WP) can normalize it if
it wants strict conformance to the six-tag list instead.

All missing-data handling happens *inside* ``build_states()`` -- there is no
build-after hook -- because ``scripts/build_case_input.py``'s
``validate_states()`` rejects any non-finite value at write time
(docs/plan_ph3.md §2.1).

Known integration gap (not fixed here, out of this WP's file allowlist)
-------------------------------------------------------------------------
``scripts/build_case_input.py`` triggers source registration via a hardcoded
``import data.era5_source`` (line ~239) and ``run_case.py`` (~line 58)
whitelists ``{era5_daily, era5_date_range, era5_merra2, None}`` for the
build-dispatch that calls ``build_case_input.py``. Neither file imports this
module or lists ``'cmip6_damip'`` in that whitelist yet, so
``run_case.py <damip_case> --step build`` will currently print "No source
block in case.yaml" and do nothing, even though
``get_source(cfg)`` would resolve correctly once this module has been
imported by *something*. This is exactly the gap ``docs/plan_ph3.md`` §2.3
assigns to WP-M4.3 ("分发泛化 + 写盘器两处小改"), a separate/later work
package that owns ``run_case.py`` and ``scripts/build_case_input.py`` in this
phase's file-ownership map (§9.3) -- this module intentionally does not
touch either file. Until WP-M4.3 lands, exercise this source directly via
``from data.cmip6_damip_source import CMIP6DamipSource`` (see
``tests/test_damip_source.py``) or via a one-line
``import data.cmip6_damip_source`` before calling
``data.source_base.get_source(cfg)``.
"""
import glob
import json
import os

import numpy as np
import yaml
from netCDF4 import Dataset

from . import cmip6_common as common
from .source_base import DataSource, register_source


# ---------------------------------------------------------------------------
# Variable groups (CMIP6 Amon short names -> pyCFRAM names handled specially
# below; everything else keeps its CMIP6 name as the dict key until final
# assembly).
# ---------------------------------------------------------------------------

REQUIRED_VARS = ('ta', 'hus', 'ts', 'ps', 'rsds', 'rsus', 'hfls', 'hfss')
CLOUD_VARS = ('cl', 'clw', 'cli')
OPTIONAL_SURF_VARS = ('rsdt', 'huss')
ALL_DISCOVER_VARS = (
    'ta', 'hus', 'ts', 'ps', 'cl', 'clw', 'cli',
    'rsdt', 'rsds', 'rsus', 'hfls', 'hfss', 'huss', 'o3',
)

VALID_O3_MODES = ('auto', 'use_model', 'climatology', 'skip')

# O3 fallback constant (kg/kg) when o3=skip or the climatology injection
# source is unavailable. NOT exact zero: RRTMG cannot tolerate an
# exact-zero H2O/O3 layer (docs/plan_ph3.md §6.1, same rationale as
# cmip6_common.fill_subsurface's HOLD strategy / mask_subsurface_layers.py).
O3_SKIP_CONSTANT_KGKG = 1e-12

# Candidate paths for the CESM 1850 O3 climatology injection source, mirrors
# scripts/inject_cesm_o3.py's CESM_O3_CANDIDATES (that script's hardcoded
# list is CESM2-4xCO2-case-specific; this is the model-agnostic DAMIP
# equivalent). A case.yaml may override via `source.o3_climatology_path`
# (an additive convenience key, not in the docs/plan_ph3.md §4.1 template,
# used mainly so tests can point at a small synthetic fixture instead of
# the real ~17MB raw_data file).
_O3_CLIMATOLOGY_CANDIDATES = (
    'raw_data/ozone_1.9x2.5_L26_1850clim_c090420.nc',
)


# ---------------------------------------------------------------------------
# configs/damip_models.d/*.yaml + configs/damip_experiments.yaml loaders
# ---------------------------------------------------------------------------

def load_model_config(model_name):
    """Load ``configs/damip_models.d/<model_name>.yaml``.

    Returns ``{}`` (with a printed warning) if no such file exists --
    graceful runtime-detection fallback, NOT a hard failure. This is the M5
    "zero-yaml new model" path (docs/plan_ph3.md §4.2): a model with no
    models.d entry still works via variant/grid auto-discovery
    (``discover_variant``/``_discover_grid_label``) and vertical-coordinate
    auto-probing (``detect_vertical`` with no override).
    """
    from core.config import PROJECT_ROOT
    path = os.path.join(PROJECT_ROOT, 'configs', 'damip_models.d', '%s.yaml' % model_name)
    if not os.path.exists(path):
        print("WARNING: cmip6_damip: no configs/damip_models.d/%s.yaml -- falling back "
              "to runtime detection (auto variant/grid discovery, vertical scheme "
              "auto-probe). See docs/plan_ph3.md §4.2." % model_name)
        return {}
    with open(path) as f:
        return yaml.safe_load(f) or {}


def list_available_models():
    """List model names with a ``configs/damip_models.d/<model>.yaml`` file."""
    from core.config import PROJECT_ROOT
    paths = sorted(glob.glob(os.path.join(PROJECT_ROOT, 'configs', 'damip_models.d', '*.yaml')))
    return [os.path.splitext(os.path.basename(p))[0] for p in paths]


def load_experiment_config(experiment_name):
    """Load the ``experiment_name`` entry from ``configs/damip_experiments.yaml``.

    Returns ``{}`` (with a printed warning) if the file or the entry is
    missing -- used only for the §5.2 single-forcing sanity-check warnings,
    never for CO2/O3 values themselves (those come from case.yaml's
    ``source:`` block).
    """
    from core.config import PROJECT_ROOT
    path = os.path.join(PROJECT_ROOT, 'configs', 'damip_experiments.yaml')
    if not os.path.exists(path):
        print("WARNING: cmip6_damip: configs/damip_experiments.yaml not found; "
              "single-forcing sanity checks (§5.2) will be skipped.")
        return {}
    with open(path) as f:
        all_exp = yaml.safe_load(f) or {}
    if experiment_name not in all_exp:
        print("WARNING: cmip6_damip: experiment '%s' not in configs/damip_experiments.yaml "
              "-- proceeding without experiment-level sanity checks." % experiment_name)
        return {}
    return all_exp[experiment_name]


def _resolve_raw_dir(raw_dir):
    if not os.path.isabs(raw_dir):
        from core.config import PROJECT_ROOT
        raw_dir = os.path.join(PROJECT_ROOT, raw_dir)
    return raw_dir


def _discover_grid_label(raw_dir, model, experiment, variant, var='ta', table_id='Amon'):
    """Like ``cmip6_common.discover_variant`` but for the grid_label token
    (CMOR filename position 5: ``<var>_<table>_<model>_<experiment>_
    <variant>_<grid>_<trange>.nc``), for the case where grid_label is not
    given in case.yaml or models.d. Kept local to this module (rather than
    added to cmip6_common.py) since it is DAMIP-source-specific plumbing,
    not reusable numerics.
    """
    pattern = '%s_%s_%s_%s_%s_*_*.nc' % (var, table_id, model, experiment, variant)
    matches = sorted(glob.glob(os.path.join(raw_dir, pattern)))
    grids = set()
    for m in matches:
        parts = os.path.basename(m).split('_')
        if len(parts) >= 6:
            grids.add(parts[5])
    return sorted(grids)


def _resolve_o3_climatology_path(source_cfg):
    from core.config import PROJECT_ROOT
    override = source_cfg.get('o3_climatology_path')
    candidates = []
    if override:
        candidates.append(override if os.path.isabs(override) else os.path.join(PROJECT_ROOT, override))
    candidates.extend(os.path.join(PROJECT_ROOT, p) for p in _O3_CLIMATOLOGY_CANDIDATES)
    for c in candidates:
        if os.path.exists(c):
            return c
    return None


# ---------------------------------------------------------------------------
# Time-aware multi-file climatology loading (docs/plan_ph3.md §5.1 steps 1-3)
# ---------------------------------------------------------------------------

def _check_month_completeness(varname, label, y0, y1, years_sel, months_sel):
    """Fail-fast month/year completeness check (docs/plan_ph3.md §5.1 step 2:
    "所选年段的月数必须 == 年数×12，缺月即 fail-fast 报哪个文件缺哪段")."""
    expected_years = list(range(int(y0), int(y1) + 1))
    got_years = set(int(y) for y in years_sel)
    missing_years = [y for y in expected_years if y not in got_years]
    if missing_years:
        raise ValueError(
            "cmip6_damip: %s (%s %d-%d) missing year(s) %s -- check raw_dir coverage"
            % (varname, label, y0, y1, missing_years))
    for y in expected_years:
        months_y = sorted(int(m) for m, yy in zip(months_sel, years_sel) if int(yy) == y)
        expected_months = list(range(1, 13))
        if months_y != expected_months:
            missing_months = [m for m in expected_months if m not in months_y]
            raise ValueError(
                "cmip6_damip: %s (%s year %d) missing month(s) %s -- check raw_dir coverage"
                % (varname, label, y, missing_months))


def _climo_pair_for_variable(paths, varname, base_years, warm_years):
    """Day-weighted annual climatology for both base and warm periods from
    one or more (time-chunked) CMIP6 files, in a single I/O pass per file.

    Returns
    -------
    (climo_base, climo_warm, lat, lon, plev)
        `plev` is None if the file has no 'plev' dimension (e.g. hybrid-level
        cl/clw/cli, or 2D surface variables).
    """
    base_datas, base_dws, base_yrs, base_mos = [], [], [], []
    warm_datas, warm_dws, warm_yrs, warm_mos = [], [], [], []
    lat = lon = plev = None

    for p in paths:
        nc = Dataset(p)
        time_var = nc.variables['time']
        units = getattr(time_var, 'units', 'days since 1850-01-01')
        calendar = getattr(time_var, 'calendar', 'standard')
        years, months, day_weights = common.decode_time(np.asarray(time_var[:]), units, calendar)
        idx_base = common.years_to_month_indices(years, base_years[0], base_years[1])
        idx_warm = common.years_to_month_indices(years, warm_years[0], warm_years[1])

        if lat is None:
            lat = np.array(nc.variables['lat'][:], dtype=np.float64)
            lon = np.array(nc.variables['lon'][:], dtype=np.float64)
            if 'plev' in nc.variables:
                plev = np.array(nc.variables['plev'][:], dtype=np.float64)

        if len(idx_base):
            base_datas.append(np.asarray(nc.variables[varname][idx_base], dtype=np.float64))
            base_dws.append(day_weights[idx_base])
            base_yrs.append(years[idx_base])
            base_mos.append(months[idx_base])
        if len(idx_warm):
            warm_datas.append(np.asarray(nc.variables[varname][idx_warm], dtype=np.float64))
            warm_dws.append(day_weights[idx_warm])
            warm_yrs.append(years[idx_warm])
            warm_mos.append(months[idx_warm])
        nc.close()

    def _finalize(datas, dws, yrs, mos, label, y0, y1):
        if not datas:
            raise ValueError(
                "cmip6_damip: %s: no data found for %s years %d-%d in %s"
                % (varname, label, y0, y1, paths))
        data_all = np.concatenate(datas, axis=0)
        dw_all = np.concatenate(dws)
        yr_all = np.concatenate(yrs)
        mo_all = np.concatenate(mos)
        _check_month_completeness(varname, label, y0, y1, yr_all, mo_all)
        idx = np.arange(len(yr_all))
        return common.annual_climo_from_monthly(data_all, idx, dw_all)

    climo_base = _finalize(base_datas, base_dws, base_yrs, base_mos, 'base', base_years[0], base_years[1])
    climo_warm = _finalize(warm_datas, warm_dws, warm_yrs, warm_mos, 'warm', warm_years[0], warm_years[1])
    return climo_base, climo_warm, lat, lon, plev


def _maybe_interp(field, plev_src_pa, plev_target_pa):
    """Skip the (identity) interpolation call when source == target plev
    exactly -- cheap correctness-preserving optimization for the common M4
    case (target grid == native CMIP6 plev19)."""
    if plev_src_pa is None:
        raise ValueError('cmip6_damip: no plev coordinate found for a plev-based variable')
    plev_src_pa = np.asarray(plev_src_pa, dtype=np.float64)
    plev_target_pa = np.asarray(plev_target_pa, dtype=np.float64)
    if plev_src_pa.shape == plev_target_pa.shape and np.array_equal(plev_src_pa, plev_target_pa):
        return np.asarray(field, dtype=np.float64)
    return common.interp_plev_to_target(field, plev_src_pa, plev_target_pa)


def _read_hybrid_coeffs(nc, vinfo, nlev_data=None):
    """Read + normalize hybrid-sigma-pressure coefficients per the naming
    convention `vinfo['scheme']` (from `common.detect_vertical`).

    Some models (confirmed against real downloaded data: IPSL-CM6A-LR)
    publish their ap/b (or a/b) coefficients at layer INTERFACES
    (nlev_data+1 values, on a dim like `klevp1`) rather than at layer
    MIDPOINTS (nlev_data values, matching the field data's actual level
    dimension, e.g. `presnivs`). `hybrid_to_plev_mass_conserving` (moved
    as-is from cesm2_cmip6_source.py, docs/plan_ph3.md §4 -- not rewritten,
    since it's already validated bit-exact against real CESM2 data) requires
    a/b the same length as the field. When `nlev_data` is given and the raw
    coefficient array is exactly one longer, adjacent interface pairs are
    averaged down to midpoint values -- exact, not an approximation: since
    layer pressure p = ap*p0 + b*ps is linear in ap/b for fixed p0/ps,
    averaging interface ap (or a) and b pairwise gives exactly the average
    of the two interface pressures, which is the standard definition of a
    layer's midpoint pressure.

    Returns
    -------
    (a_eff, b_eff, p0_eff, flipped) -- `hybrid_to_plev_mass_conserving`
    requires TOA->sfc order (b ascending from ~0 at TOA to ~1 at surface).
    CESM2's native hyam/hybm are already stored this way (the function's
    original, only-tested convention), but other models (confirmed:
    IPSL-CM6A-LR) store levels surface->TOA instead (b descending from ~1 to
    ~0). Direction is detected from b and the coefficient arrays are
    reversed if needed; `flipped` tells the caller to reverse the field data
    (cl/clw/cli) the same way before calling that function -- reversing only
    one side silently scrambles the profile (this crashed with IndexError in
    an earlier attempt at this fix for a different reason, then produced a
    plausible-looking but physically wrong near-zero cloud field once that
    was fixed, until this order mismatch was found by sanity-checking
    against the real raw data's actual cloud fraction magnitude).
    """
    def _to_midlayer(raw):
        raw = np.asarray(raw, dtype=np.float64)
        if nlev_data is not None and raw.shape[0] == nlev_data + 1:
            return 0.5 * (raw[:-1] + raw[1:])
        return raw

    if vinfo['scheme'] == 'hybrid_ab_p0':
        a_raw = _to_midlayer(nc.variables[vinfo['a']][:])
        b_raw = _to_midlayer(nc.variables[vinfo['b']][:])
        p0_raw = float(np.asarray(nc.variables[vinfo['p0']][:]).reshape(-1)[0])
        a_eff, b_eff, p0_eff = common.normalize_hybrid_coeffs('hybrid_ab_p0', a=a_raw, b=b_raw, p0=p0_raw)
    elif vinfo['scheme'] == 'hybrid_ap_b':
        ap_raw = _to_midlayer(nc.variables[vinfo['ap']][:])
        b_raw = _to_midlayer(nc.variables[vinfo['b']][:])
        a_eff, b_eff, p0_eff = common.normalize_hybrid_coeffs('hybrid_ap_b', ap=ap_raw, b=b_raw)
    else:
        raise ValueError('cmip6_damip: unrecognized hybrid scheme %r' % vinfo['scheme'])

    flipped = bool(b_eff[0] > b_eff[-1])
    if flipped:
        a_eff = a_eff[::-1].copy()
        b_eff = b_eff[::-1].copy()
    return a_eff, b_eff, p0_eff, flipped


@register_source('cmip6_damip')
class CMIP6DamipSource(DataSource):
    """Build CFRAM states from a CMIP6 DAMIP single-forcing experiment.

    Expects a case.yaml ``source:`` block per docs/plan_ph3.md §4.1:

        source:
          type: cmip6_damip
          model: IPSL-CM6A-LR
          experiment: hist-aer
          variant: r1i1p1f1        # optional
          grid_label: gr           # optional
          raw_dir: raw_data/cmip6_damip
          base_years: [1850, 1859]
          warm_years: [2011, 2020]
          co2: {source: constant, base_ppmv: 284.7, perturbed_ppmv: 284.7}
          o3: auto                 # auto|use_model|climatology|skip
          aerosol: {source: zero}
    """

    def build_states(self):
        scfg = self.source_cfg
        model = scfg['model']
        experiment = scfg['experiment']
        raw_dir = _resolve_raw_dir(scfg['raw_dir'])
        base_years = [int(scfg['base_years'][0]), int(scfg['base_years'][1])]
        warm_years = [int(scfg['warm_years'][0]), int(scfg['warm_years'][1])]

        model_cfg = load_model_config(model)
        exp_cfg = load_experiment_config(experiment)

        # ── variant / grid_label resolution (case.yaml > models.d > glob) ──
        variant = scfg.get('variant') or model_cfg.get('default_variant')
        grid_label_hint = scfg.get('grid_label') or model_cfg.get('default_grid')
        if not variant:
            variants = common.discover_variant(raw_dir, model, experiment, grid=grid_label_hint)
            if not variants:
                raise FileNotFoundError(
                    "cmip6_damip: no variant found for %s/%s under %s (tried table_id=Amon, var=ta)"
                    % (model, experiment, raw_dir))
            variant = variants[0]
            print("WARNING: cmip6_damip: no variant configured for %s -> auto-discovered '%s' "
                  "(available: %s)" % (model, variant, variants))

        grid_label = grid_label_hint
        if not grid_label:
            grids = _discover_grid_label(raw_dir, model, experiment, variant)
            if grids:
                grid_label = grids[0]
                if len(grids) > 1:
                    print("WARNING: cmip6_damip: multiple grid_labels found for %s/%s/%s: %s "
                          "-> using '%s'" % (model, experiment, variant, grids, grid_label))

        print("cmip6_damip: model=%s experiment=%s variant=%s grid=%s raw_dir=%s"
              % (model, experiment, variant, grid_label, raw_dir))

        # docs/plan_ph3.md §5.1 step 1: discover_files/discover_variant
        files = common.discover_files(raw_dir, ALL_DISCOVER_VARS, model=model,
                                       experiment=experiment, variant=variant, grid=grid_label)

        missing_required = [v for v in REQUIRED_VARS if v not in files]
        if missing_required:
            raise FileNotFoundError(
                "cmip6_damip: required variable(s) %s not found for %s/%s/%s/%s under %s"
                % (missing_required, model, experiment, variant, grid_label, raw_dir))

        # Target plev: case.yaml grid.pressure_levels is TOA->surface hPa,
        # ascending (pyCFRAM internal convention, e.g. [1, 5, ..., 1000]).
        # NetCDF write convention (docs/input_spec.md) is surface->TOA
        # (descending) -- matches native CMIP6 plev order directly.
        target_plev_toa2sfc_hpa = np.asarray(self.cfg['grid']['pressure_levels'], dtype=np.float64)
        target_lev_hpa = target_plev_toa2sfc_hpa[::-1].copy()          # sfc->TOA, hPa
        target_plev_pa_sfc2toa = target_lev_hpa * 100.0                # sfc->TOA, Pa
        target_plev_pa_toa2sfc = target_plev_pa_sfc2toa[::-1].copy()   # TOA->sfc, Pa (hybrid proj needs this order)

        provenance = {
            'source_type': 'cmip6_damip',
            'model': model, 'experiment': experiment,
            'variant': variant, 'grid_label': grid_label,
            'raw_dir': raw_dir,
            'base_years': base_years, 'warm_years': warm_years,
            'target_plev_hpa_toa2sfc': target_plev_toa2sfc_hpa.tolist(),
            'processes': {},
            'sanity_checks': {},
        }

        base_state = {}
        warm_state = {}

        # ── docs/plan_ph3.md §5.1 steps 2-3, 5: ta/hus (plev, interp if needed) ──
        climo_ta_b, climo_ta_w, lat, lon, plev_ta = _climo_pair_for_variable(
            files['ta'], 'ta', base_years, warm_years)
        climo_hus_b, climo_hus_w, _, _, plev_hus = _climo_pair_for_variable(
            files['hus'], 'hus', base_years, warm_years)

        base_state['ta_lay'] = _maybe_interp(climo_ta_b, plev_ta, target_plev_pa_sfc2toa)
        warm_state['ta_lay'] = _maybe_interp(climo_ta_w, plev_ta, target_plev_pa_sfc2toa)
        base_state['q'] = _maybe_interp(climo_hus_b, plev_hus if plev_hus is not None else plev_ta,
                                         target_plev_pa_sfc2toa)
        warm_state['q'] = _maybe_interp(climo_hus_w, plev_hus if plev_hus is not None else plev_ta,
                                         target_plev_pa_sfc2toa)

        nlat, nlon = len(lat), len(lon)
        nlev_t = len(target_lev_hpa)
        shape3d = (nlev_t, nlat, nlon)
        shape2d = (nlat, nlon)

        # ── ts / ps (required 2D) ──
        climo_ts_b, climo_ts_w, _, _, _ = _climo_pair_for_variable(files['ts'], 'ts', base_years, warm_years)
        climo_ps_b, climo_ps_w, _, _, _ = _climo_pair_for_variable(files['ps'], 'ps', base_years, warm_years)
        base_state['ts'] = climo_ts_b
        warm_state['ts'] = climo_ts_w
        base_state['ps'] = climo_ps_b
        warm_state['ps'] = climo_ps_w

        # ── albedo (rsus/rsds -- docs/plan_ph3.md §5.1 step 6) ──
        climo_rsds_b, climo_rsds_w, _, _, _ = _climo_pair_for_variable(files['rsds'], 'rsds', base_years, warm_years)
        climo_rsus_b, climo_rsus_w, _, _, _ = _climo_pair_for_variable(files['rsus'], 'rsus', base_years, warm_years)
        base_state['albedo'] = common.compute_albedo(climo_rsus_b, climo_rsds_b)
        warm_state['albedo'] = common.compute_albedo(climo_rsus_w, climo_rsds_w)
        provenance['processes']['albedo'] = {
            'status': 'ACTIVE', 'reason': 'rsus/rsds both present -> compute_albedo'}

        # ── solar (rsdt if present, else analytic_solar -- §5.1 step 6, §6.2) ──
        if 'rsdt' in files:
            climo_rsdt_b, climo_rsdt_w, _, _, _ = _climo_pair_for_variable(
                files['rsdt'], 'rsdt', base_years, warm_years)
            base_state['solar'] = climo_rsdt_b
            warm_state['solar'] = climo_rsdt_w
            solar_status = 'ACTIVE'
            solar_reason = 'model rsdt present'
        else:
            solar_1d = common.analytic_solar(lat)
            solar_2d = np.broadcast_to(solar_1d[:, None], shape2d).copy()
            base_state['solar'] = solar_2d
            warm_state['solar'] = solar_2d.copy()
            solar_status = 'ANALYTIC'
            solar_reason = ('rsdt not published for this model; used analytic_solar(lat) fallback '
                             '-- hist-aer freezes the solar constant so frc_solar==0 by construction')
        provenance['processes']['solar'] = {'status': solar_status, 'reason': solar_reason}

        # ── huss (optional 2D; state key added even though the current
        #    write_surf_nc (pre-WP-M4.3) does not yet emit it -- see module
        #    docstring's "Known integration gap") ──
        if 'huss' in files:
            climo_huss_b, climo_huss_w, _, _, _ = _climo_pair_for_variable(
                files['huss'], 'huss', base_years, warm_years)
            base_state['huss'] = climo_huss_b
            warm_state['huss'] = climo_huss_w

        # ── flux (hfls/hfss -> nonrad forcing, §5.1 step 9, §2.5 sign law) ──
        climo_hfls_b, climo_hfls_w, _, _, _ = _climo_pair_for_variable(files['hfls'], 'hfls', base_years, warm_years)
        climo_hfss_b, climo_hfss_w, _, _, _ = _climo_pair_for_variable(files['hfss'], 'hfss', base_years, warm_years)
        # CMIP6 hfls/hfss are positive-UP; CFRAM forcing convention needs
        # frc = Delta(downward surface energy flux) = -(warm - base).
        nonrad = {
            'lhflx': -(climo_hfls_w - climo_hfls_b),
            'shflx': -(climo_hfss_w - climo_hfss_b),
        }
        provenance['processes']['flux'] = {
            'status': 'ACTIVE',
            'reason': 'hfls/hfss both present; frc_lhflx=-(warm.hfls-base.hfls), '
                      'frc_shflx=-(warm.hfss-base.hfss)'}

        # ── cloud (cl/clw/cli -- §5.1 step 4, §6.1 decision tree) ──
        cloud_ok = all(v in files for v in CLOUD_VARS)
        cloud_reason = None
        if cloud_ok:
            climo_cl_b, climo_cl_w, _, _, _ = _climo_pair_for_variable(files['cl'], 'cl', base_years, warm_years)
            climo_clw_b, climo_clw_w, _, _, _ = _climo_pair_for_variable(files['clw'], 'clw', base_years, warm_years)
            climo_cli_b, climo_cli_w, _, _, _ = _climo_pair_for_variable(files['cli'], 'cli', base_years, warm_years)

            cl_b_frac = climo_cl_b / 100.0   # %->fraction
            cl_w_frac = climo_cl_w / 100.0

            if (np.all(np.isnan(cl_b_frac)) or np.all(np.isnan(cl_w_frac)) or
                    np.all(np.isnan(climo_clw_b)) or np.all(np.isnan(climo_clw_w)) or
                    np.all(np.isnan(climo_cli_b)) or np.all(np.isnan(climo_cli_w))):
                cloud_ok = False
                cloud_reason = 'cl/clw/cli files present but data is all-NaN/fill'

        if cloud_ok:
            with Dataset(files['cl'][0]) as nc0:
                cl_var = nc0.variables['cl']
                formula_terms = getattr(cl_var, 'formula_terms', None)
                if formula_terms is None and 'lev' in nc0.variables:
                    formula_terms = getattr(nc0.variables['lev'], 'formula_terms', None)

                vertical_cfg = model_cfg.get('vertical') or {'scheme': 'auto'}
                override = None if vertical_cfg.get('scheme', 'auto') == 'auto' else vertical_cfg
                is_hybrid = (override is not None) or (formula_terms is not None)

                if is_hybrid:
                    vinfo = common.detect_vertical(formula_terms=formula_terms, override=override)
                    a_eff, b_eff, p0_eff, coeffs_flipped = _read_hybrid_coeffs(
                        nc0, vinfo, nlev_data=cl_b_frac.shape[0])

                    def _to_target(field_frac, ps_2d_pa):
                        # hybrid_to_plev_mass_conserving wants TOA->sfc order
                        # for both the field and a_eff/b_eff. _read_hybrid_coeffs
                        # already normalized the coefficients to that order and
                        # reports whether a reversal was needed -- the field
                        # (native file order, matching the coefficients' native
                        # order) must be reversed the same way, or the profile
                        # is silently scrambled (see _read_hybrid_coeffs docstring).
                        field_toa2sfc = field_frac[::-1] if coeffs_flipped else field_frac
                        proj_toa2sfc = common.hybrid_to_plev_mass_conserving(
                            field_toa2sfc, a_eff, b_eff, p0_eff, ps_2d_pa, target_plev_pa_toa2sfc)
                        return proj_toa2sfc[::-1]

                    base_state['camt'] = np.clip(_to_target(cl_b_frac, climo_ps_b), 0.0, 1.0)
                    warm_state['camt'] = np.clip(_to_target(cl_w_frac, climo_ps_w), 0.0, 1.0)
                    base_state['cliq'] = _to_target(climo_clw_b, climo_ps_b)
                    warm_state['cliq'] = _to_target(climo_clw_w, climo_ps_w)
                    base_state['cice'] = _to_target(climo_cli_b, climo_ps_b)
                    warm_state['cice'] = _to_target(climo_cli_w, climo_ps_w)
                else:
                    plev_cl = np.array(nc0.variables['plev'][:], dtype=np.float64) if 'plev' in nc0.variables else plev_ta
                    base_state['camt'] = np.clip(_maybe_interp(cl_b_frac, plev_cl, target_plev_pa_sfc2toa), 0.0, 1.0)
                    warm_state['camt'] = np.clip(_maybe_interp(cl_w_frac, plev_cl, target_plev_pa_sfc2toa), 0.0, 1.0)
                    base_state['cliq'] = _maybe_interp(climo_clw_b, plev_cl, target_plev_pa_sfc2toa)
                    warm_state['cliq'] = _maybe_interp(climo_clw_w, plev_cl, target_plev_pa_sfc2toa)
                    base_state['cice'] = _maybe_interp(climo_cli_b, plev_cl, target_plev_pa_sfc2toa)
                    warm_state['cice'] = _maybe_interp(climo_cli_w, plev_cl, target_plev_pa_sfc2toa)
            cloud_status = 'ACTIVE'
            cloud_reason = cloud_reason or 'cl/clw/cli all present and not all-NaN'
        else:
            base_state['camt'] = np.zeros(shape3d)
            warm_state['camt'] = np.zeros(shape3d)
            base_state['cliq'] = np.zeros(shape3d)
            warm_state['cliq'] = np.zeros(shape3d)
            base_state['cice'] = np.zeros(shape3d)
            warm_state['cice'] = np.zeros(shape3d)
            cloud_status = 'SKIPPED'
            missing_cloud = [v for v in CLOUD_VARS if v not in files]
            cloud_reason = cloud_reason or ('cl/clw/cli not all published (missing: %s)' % missing_cloud)
        provenance['processes']['cloud'] = {'status': cloud_status, 'reason': cloud_reason}

        # ── O3 (§5.1 step 7, §6.1 decision tree) ──
        o3_mode = scfg.get('o3', 'auto')
        if o3_mode not in VALID_O3_MODES:
            raise ValueError("cmip6_damip: source.o3 must be one of %s, got %r" % (VALID_O3_MODES, o3_mode))
        o3_available = 'o3' in files

        if o3_mode == 'use_model' or (o3_mode == 'auto' and o3_available):
            if not o3_available:
                raise FileNotFoundError(
                    "cmip6_damip: o3='use_model' requested but %s/%s/%s has no o3 data"
                    % (model, experiment, variant))
            climo_o3_b, climo_o3_w, lat_o3, lon_o3, plev_o3 = _climo_pair_for_variable(
                files['o3'], 'o3', base_years, warm_years)
            # Some models publish o3 on a coarser native horizontal grid than
            # the model's other Amon variables (confirmed: MRI-ESM2-0's o3 is
            # 64x128 while ta/hus/ts/ps/cl are all 160x320 -- CMIP6 has no
            # requirement that every variable in a table_id share one grid).
            # Regrid onto the reference (ta) grid before anything downstream
            # assumes every field shares one lat/lon.
            if lat_o3.shape != lat.shape or lon_o3.shape != lon.shape or \
                    not np.array_equal(lat_o3, lat) or not np.array_equal(lon_o3, lon):
                climo_o3_b = common.regrid_horizontal_bilinear(lat_o3, lon_o3, climo_o3_b, lat, lon)
                climo_o3_w = common.regrid_horizontal_bilinear(lat_o3, lon_o3, climo_o3_w, lat, lon)
            plev_o3_use = plev_o3 if plev_o3 is not None else plev_ta
            o3_b_mmr = _maybe_interp(climo_o3_b, plev_o3_use, target_plev_pa_sfc2toa) * common.VMR_TO_MMR
            o3_w_mmr = _maybe_interp(climo_o3_w, plev_o3_use, target_plev_pa_sfc2toa) * common.VMR_TO_MMR
            base_state['o3'] = o3_b_mmr
            warm_state['o3'] = o3_w_mmr
            o3_status = 'MODEL'
            o3_reason = 'model o3 (mol/mol) x 48/29 -> kg/kg'
        elif o3_mode in ('auto', 'climatology'):
            o3_path = _resolve_o3_climatology_path(scfg)
            if o3_path is None:
                base_state['o3'] = np.full(shape3d, O3_SKIP_CONSTANT_KGKG)
                warm_state['o3'] = np.full(shape3d, O3_SKIP_CONSTANT_KGKG)
                o3_status = 'SKIPPED'
                o3_reason = ('o3=%s but no model o3 and no climatology injection source found -> '
                             'constant %.0e kg/kg (RRTMG cannot tolerate exact zero)'
                             % (o3_mode, O3_SKIP_CONSTANT_KGKG))
            else:
                o3_3d = common.o3_climatology(o3_path, target_lev_hpa, lat, lon)
                base_state['o3'] = o3_3d
                warm_state['o3'] = o3_3d.copy()
                o3_status = 'CLIMATOLOGY'
                o3_reason = ('no model o3 (or o3=climatology forced); injected CESM 1850 climatology '
                             'from %s (same field both states -> frc_o3==0)' % o3_path)
        else:  # o3_mode == 'skip'
            base_state['o3'] = np.full(shape3d, O3_SKIP_CONSTANT_KGKG)
            warm_state['o3'] = np.full(shape3d, O3_SKIP_CONSTANT_KGKG)
            o3_status = 'SKIPPED'
            o3_reason = 'o3=skip requested -> constant %.0e kg/kg' % O3_SKIP_CONSTANT_KGKG
        provenance['processes']['o3'] = {'status': o3_status, 'reason': o3_reason}

        # ── CO2 (§5.1 step 7) ──
        co2_base = self.get_co2('base')
        co2_pert = self.get_co2('perturbed')
        base_state['co2'] = np.full(shape3d, co2_base, dtype=np.float64)
        warm_state['co2'] = np.full(shape3d, co2_pert, dtype=np.float64)

        # ── Aerosol (§6.1: hist-aer Amon has no 3D aerosol mmr -> zero) ──
        base_state.update(self.get_aerosol(shape3d))
        warm_state.update(self.get_aerosol(shape3d))
        aer_src = scfg.get('aerosol', {}).get('source', 'zero')
        provenance['processes']['aerosol'] = {
            'status': 'SKIPPED' if aer_src == 'zero' else 'ACTIVE',
            'reason': ('hist-aer Amon has no 3D aerosol mixing ratio; source=zero -> '
                       'runner auto-detects skip_aerosol' if aer_src == 'zero'
                       else 'aerosol.source=%s' % aer_src)}

        # ── lat/lon coordinates for both states (pre-normalize) ──
        base_state['lat'] = lat
        base_state['lon'] = lon
        base_state['lev'] = target_lev_hpa
        warm_state['lat'] = lat
        warm_state['lon'] = lon
        warm_state['lev'] = target_lev_hpa

        # ── docs/plan_ph3.md §5.1 step 8: fill_subsurface (no NaN may
        #    survive past this point -- validate_states rejects them) ──
        for state in (base_state, warm_state):
            ps_hpa = state['ps'] / 100.0
            fields_in = {
                'ta_lay': state['ta_lay'], 'q': state['q'], 'o3': state['o3'],
                'camt': state['camt'], 'cliq': state['cliq'], 'cice': state['cice'],
                'bc': state['bc'], 'ocphi': state['ocphi'], 'ocpho': state['ocpho'],
                'sulf': state['sulf'], 'ss': state['ss'], 'dust': state['dust'],
            }
            out, _summary = common.fill_subsurface(fields_in, target_lev_hpa, ps_hpa, state['ts'])
            state.update(out)

        # ── docs/plan_ph3.md §5.1 step 6/normalize: lon 0-360 ascending,
        #    lat S->N ascending, applied consistently to every field ──
        for state in (base_state, warm_state):
            for key in list(state.keys()):
                if key in ('lat', 'lon', 'lev'):
                    continue
                arr = state[key]
                if isinstance(arr, np.ndarray) and arr.ndim >= 2:
                    _, _, state[key] = common.normalize_grid(lat, lon, arr)
        lat_n, lon_n = common.normalize_grid(lat, lon)
        base_state['lat'] = lat_n
        base_state['lon'] = lon_n
        warm_state['lat'] = lat_n
        warm_state['lon'] = lon_n

        # ── docs/plan_ph3.md §5.2: single-forcing consistency sanity
        #    checks (warnings only, never a hard gate) ──
        sanity = provenance['sanity_checks']
        co2_equal = bool(np.isclose(co2_base, co2_pert))
        sanity['co2_equal'] = co2_equal
        exp_co2_kind = exp_cfg.get('co2') if exp_cfg else None
        if exp_co2_kind == 'fixed_1850' and not co2_equal:
            print("WARNING: cmip6_damip: experiment '%s' expects co2 fixed at 1850 but "
                  "base_ppmv(%.2f) != perturbed_ppmv(%.2f)"
                  % (experiment, co2_base * 1e6, co2_pert * 1e6))
        elif exp_co2_kind == 'time_varying' and co2_equal:
            print("WARNING: cmip6_damip: experiment '%s' expects time-varying CO2 but "
                  "base_ppmv == perturbed_ppmv" % experiment)

        if o3_status == 'MODEL':
            denom = max(float(np.nanmean(np.abs(base_state['o3']))), 1e-30)
            o3_rel_diff = float(np.nanmean(np.abs(warm_state['o3'] - base_state['o3']))) / denom
            sanity['o3_rel_diff'] = o3_rel_diff
            exp_o3_default = exp_cfg.get('o3_default') if exp_cfg else None
            if exp_o3_default not in (None, 'use_model') and o3_rel_diff > 0.01:
                print("WARNING: cmip6_damip: o3 domain-mean relative diff %.3f%% exceeds 1%% for "
                      "an experiment where O3 is expected frozen" % (o3_rel_diff * 100.0))

        if solar_status == 'ACTIVE':
            solar_diff = float(np.nanmax(np.abs(warm_state['solar'] - base_state['solar'])))
            sanity['solar_max_abs_diff_wm2'] = solar_diff
            if solar_diff > 0.5:
                print("WARNING: cmip6_damip: solar (rsdt) max abs diff %.3f W/m2 exceeds 0.5 W/m2 "
                      "-- expected near-zero if the solar constant is frozen" % solar_diff)

        self.provenance = provenance
        self._write_provenance()

        return base_state, warm_state, nonrad

    def _write_provenance(self):
        """Write provenance.json next to the case's input NCs.

        scripts/build_case_input.py has no provenance-dump hook (a later
        WP owns that file), so this writes directly -- see module
        docstring's "Known integration gap".
        """
        case_dir = self.cfg.get('_case_dir')
        if not case_dir:
            return
        input_dir = os.path.join(case_dir, 'input')
        os.makedirs(input_dir, exist_ok=True)
        path = os.path.join(input_dir, 'provenance.json')
        with open(path, 'w') as f:
            json.dump(self.provenance, f, indent=2)
        print("  Wrote %s" % path)
