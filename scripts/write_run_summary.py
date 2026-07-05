#!/usr/bin/env python3
"""Write cases/<case>/output/<case>.summary.txt after a CFRAM run.

docs/plan_ph3.md §6.3 (WP-M4.4). Standalone post-processing script -- it is
NOT imported by, and does not couple to, scripts/run_parallel_python.py
(that file is already-accepted Phase 2 code; per the plan we read its
*output* NetCDF only, never its internals). run_case.py calls this script
via the same subprocess.call `run_step()` helper used for build/plot, at
the end of the `run` step only.

Reads:
    cases/<case>/input/provenance.json   (written by build_states(); DAMIP
        source plugin only -- may be absent for non-DAMIP cases (ERA5,
        cesm2_cmip6-collab, idealized climlab cases, ...). This script
        degrades gracefully rather than failing the run step for those.)
    cases/<case>/output/cfram_result.nc  (written by the CFRAM run step;
        has dT_observed + whichever dT_<process> variables the run
        actually produced -- the set varies per case/config, see
        scripts/run_parallel_python.py's `output_terms` filter.)
    cases/<case>/input/*.nc              (for the md5 footer)
    cases/<case>/case.yaml               (radiation.scheme, via core.config)
    configs/damip_experiments.yaml       (optional; only used to judge
        whether a given sanity_checks value is expected to be near-zero
        for this experiment_id -- e.g. hist-nat *should* show a solar
        diff since solar is that experiment's own varying forcing).

Writes:
    cases/<case>/output/<case>.summary.txt

Usage:
    python3 scripts/write_run_summary.py --case damip_ipsl_histaer
"""
import os
import sys
import json
import glob
import hashlib
import argparse

import numpy as np
from netCDF4 import Dataset

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.config import load_case, PROJECT_ROOT

FILL_ABS_THRESHOLD = 900.0   # matches scripts/compute_lr_attribution.py / plot_closure_profile.py
O3_REL_DIFF_OK = 0.01        # docs/plan_ph3.md §5.2: "两态均相对差 < 1%"
SOLAR_DIFF_OK_WM2 = 0.5      # docs/plan_ph3.md §5.2: "|solar_warm - solar_base| 域均 < 0.5 W/m^2"

M3_DOC_REL = 'docs/m3_methodology_comparison.md'

# ---------------------------------------------------------------------------
# Additivity: which dT_<term> variables to sum, and which to skip because
# they are the *same* physical quantity expressed a second (redundant) way.
#
# All of the exclusions below are read straight off the variable list in
# cfram_result.nc -- this module never imports run_parallel_python.py, it
# just knows the *naming convention* that module documents in its own
# comments (cited below) and in session_log.md.
# ---------------------------------------------------------------------------

# Raw full-state radiative-response diagnostics used internally to *derive*
# atmdyn/sfcdyn/ocndyn (run_parallel_python.py ~L423-467: "F_dyn = -F_rad",
# "rad_1d_full call uses full warm state"); not independent physical
# processes in their own right. Always excluded when present.
_DIAGNOSTIC_ONLY_TERMS = frozenset({'warm', 'full'})

# name -> set of finer-grained dT_<x> names that reconstruct it (near-)exactly.
# If ALL of the components are present as separate variables, the aggregate
# name is dropped from the additivity sum in favour of the components.
_EXACT_AGGREGATES = {
    # OLD CFRAM "dt_dyn": full-column drdt^-1 . frc_full = atmdyn + sfcdyn
    # (run_parallel_python.py comment: "dry ≡ OLD dt_dyn ... = atmdyn+sfcdyn").
    'dry': frozenset({'atmdyn', 'sfcdyn'}),
    # OLD CFRAM "dt_sfc_dyn": ocndyn + lhflx + shflx, verified exact to
    # machine precision after the surface-only-forcing fix
    # (session_log.md 2026-06-01: "residual ... max 3.6e-14 K").
    'sfcdyn': frozenset({'ocndyn', 'lhflx', 'shflx'}),
    # Exact LW/SW split of bulk cloud forcing, same RRTMG call
    # (session_log.md 2026-04-17: "dT_cloud == dT_cloud_lw + dT_cloud_sw,
    # rel=1.5e-15").
    'cloud': frozenset({'cloud_lw', 'cloud_sw'}),
}

# Per-species aerosol breakdown (two possible vocabularies: legacy path-B
# 5-term OC-merged, or Fortran 6-term OC-split) sums to ~bulk dT_aerosol but
# with a real (~2.5-4%) non-linear residual (session_log.md 2026-04-18:
# "dT_aerosol=+1.22 K, 6-species sum=+1.19 K, 2.5% nonlinear residual") --
# NOT an exact identity like the ones above, but still "the same quantity
# two ways": if the bulk term is present we prefer it and drop the species
# breakdown from the additivity sum (species stay in the raw NC either way).
_AEROSOL_SPECIES_NAMES = frozenset(
    {'bc', 'oc', 'ocphi', 'ocpho', 'sulf', 'ss', 'seas', 'dust'})


def select_additive_terms(all_dT_names):
    """Pick a non-double-counting subset of `dT_<term>` names to sum.

    Returns
    -------
    terms : sorted list[str] -- the process names to sum (without the
        'dT_' prefix), e.g. ['atmdyn', 'ocndyn', 'q', ...].
    excluded : dict[str, str] -- name (or comma-joined group) -> reason it
        was left out of the sum.
    notes : list[str] -- cases this function could not confidently resolve
        (per task instructions: flag ambiguity rather than silently guess).
    """
    names = set(all_dT_names) - {'observed'}
    excluded = {}
    notes = []

    for term in sorted(_DIAGNOSTIC_ONLY_TERMS):
        if term in names:
            names.discard(term)
            excluded[term] = (
                "diagnostic full-state radiative-response term used "
                "internally to derive atmdyn/sfcdyn/ocndyn "
                "(scripts/run_parallel_python.py), not an independent "
                "physical process")

    if 'dry' in names:
        names.discard('dry')
        excluded['dry'] = (
            "= dT_atmdyn + dT_sfcdyn (OLD CFRAM 'dt_dyn', full-column "
            "aggregate) -- components counted separately")

    if 'sfcdyn' in names:
        components = _EXACT_AGGREGATES['sfcdyn']
        if components <= names:
            names.discard('sfcdyn')
            excluded['sfcdyn'] = (
                "= dT_ocndyn + dT_lhflx + dT_shflx (exact to machine "
                "precision, session_log.md 2026-06-01) -- components "
                "counted separately")
        else:
            notes.append(
                "dT_sfcdyn is present but not all of "
                "dT_ocndyn/dT_lhflx/dT_shflx are -- kept dT_sfcdyn in the "
                "sum as-is; cannot confirm it is a pure duplicate for this "
                "case's variable set")

    if 'cloud' in names and _EXACT_AGGREGATES['cloud'] <= names:
        names.discard('cloud_lw')
        names.discard('cloud_sw')
        excluded['cloud_lw+cloud_sw'] = (
            "= dT_cloud split into LW/SW components (exact, "
            "session_log.md 2026-04-17) -- bulk dT_cloud counted instead")

    present_species = names & _AEROSOL_SPECIES_NAMES
    if present_species:
        if 'aerosol' in names:
            for s in present_species:
                names.discard(s)
            excluded[",".join(sorted(present_species))] = (
                "per-species aerosol breakdown of dT_aerosol (~2.5-4% "
                "non-linear residual vs bulk, session_log.md 2026-04-18) "
                "-- bulk dT_aerosol counted instead")
        else:
            notes.append(
                "per-species aerosol terms present (%s) without a bulk "
                "dT_aerosol -- summing the species directly; this is only "
                "an approximate (not exact) reconstruction of the aerosol "
                "process, per session_log.md 2026-04-18"
                % ", ".join(sorted(present_species)))

    return sorted(names), excluded, notes


def _mask_fill(arr):
    arr = np.asarray(arr, dtype=np.float64)
    return np.where(np.abs(arr) > FILL_ABS_THRESHOLD, np.nan, arr)


def compute_additivity(nc_path):
    """dT_observed[sfc] - Sum(dT_X[sfc]) domain mean/max, plus term bookkeeping.

    Surface-level convention: the *last* index along the `lev` dimension.
    Confirmed against two other Phase-2 modules that already index the
    surface this way:
      - scripts/plot_closure_profile.py: `nlev_atm = len(lev) - 1  # last
        index is surface` and reads `dT_obs[-1]` as the surface value.
      - scripts/compute_lr_attribution.py: `skin_idx = np.argmax(lev >
        1005.0) if np.any(lev > 1005.0) else lev.size - 1` -- i.e. it looks
        for an explicit "skin" marker level above 1005 hPa but falls back
        to the last index, and in the one real fixture inspected
        (tests/data/smoke/cfram_result_mini.nc) the skin marker (1013 hPa)
        *is* the last element, so both conventions agree in practice.
    """
    with Dataset(nc_path) as d:
        lev = np.array(d.variables['lev'][:], dtype=np.float64)
        nlat = d.dimensions['lat'].size
        nlon = d.dimensions['lon'].size
        nlev = d.dimensions['lev'].size
        sfc_idx = -1  # last index along lev == surface (see docstring)

        all_dT_names = sorted(
            v[3:] for v in d.variables if v.startswith('dT_') and v != 'dT_observed')

        if 'dT_observed' not in d.variables:
            return {
                'ok': False,
                'reason': "dT_observed not found in %s" % nc_path,
                'nlat': nlat, 'nlon': nlon, 'nlev': nlev,
            }

        obs_sfc = _mask_fill(d.variables['dT_observed'][sfc_idx, :, :])

        terms, excluded, notes = select_additive_terms(all_dT_names)

        total = np.zeros_like(obs_sfc)
        used_terms = []
        for t in terms:
            vname = 'dT_' + t
            arr_sfc = _mask_fill(d.variables[vname][sfc_idx, :, :])
            if np.all(np.isnan(arr_sfc)):
                # Present in the file (created) but never actually computed
                # for this case/config (e.g. output_terms filtering, or an
                # optional forcing that wasn't supplied) -- not a real
                # contributor, don't count it and don't claim we did.
                continue
            total = total + np.nan_to_num(arr_sfc, nan=0.0)
            used_terms.append(t)

        residual = obs_sfc - total
        return {
            'ok': True,
            'nlat': nlat, 'nlon': nlon, 'nlev': nlev,
            'all_dT_names': all_dT_names,
            'used_terms': used_terms,
            'excluded': excluded,
            'notes': notes,
            'resid_mean': float(np.nanmean(residual)),
            'resid_max_abs': float(np.nanmax(np.abs(residual))),
        }


def md5_of_file(path, chunk_size=1 << 20):
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(chunk_size), b''):
            h.update(chunk)
    return h.hexdigest()


def load_damip_experiments_cfg():
    path = os.path.join(PROJECT_ROOT, 'configs', 'damip_experiments.yaml')
    if not os.path.exists(path):
        return {}
    import yaml
    with open(path) as f:
        return yaml.safe_load(f) or {}


def judge_sanity_checks(sanity_checks, experiment, damip_experiments_cfg):
    """One-line PASS/WARNING judgment per docs/plan_ph3.md §5.2.

    Where possible, the expectation is looked up from
    configs/damip_experiments.yaml (e.g. hist-nat legitimately *should* show
    a large solar diff, since solar is that experiment's own varying
    forcing) rather than assuming every experiment is CO2/O3/solar-frozen.
    Never raises/exits on an unexpected value -- worst case is a WARNING
    line, per the task's "don't fail/crash" instruction.
    """
    exp_cfg = damip_experiments_cfg.get(experiment) if experiment else None
    varying = set(exp_cfg.get('varying', [])) if exp_cfg else None

    rows = []

    if 'co2_equal' in sanity_checks:
        co2_equal = sanity_checks['co2_equal']
        co2_mode = exp_cfg.get('co2') if exp_cfg else None
        if co2_mode == 'fixed_1850':
            expect = True
        elif co2_mode == 'time_varying':
            expect = False
        else:
            expect = None
        if expect is None:
            judgment = ("(no CO2 expectation on file for experiment=%r; "
                        "raw value=%r)" % (experiment, co2_equal))
        elif co2_equal == expect:
            judgment = "PASS (matches expected co2_equal=%r for experiment=%r)" % (
                expect, experiment)
        else:
            judgment = "WARNING: expected co2_equal=%r for experiment=%r, got %r" % (
                expect, experiment, co2_equal)
        rows.append(('co2_equal', co2_equal, judgment))

    if 'o3_rel_diff' in sanity_checks:
        o3_rel_diff = sanity_checks['o3_rel_diff']
        o3_is_varying = bool(varying) and 'o3' in varying
        if o3_is_varying:
            judgment = ("(o3 is the varying forcing for experiment=%r -- a "
                        "large diff is expected here, not checked)" % experiment)
        elif o3_rel_diff is not None and o3_rel_diff < O3_REL_DIFF_OK:
            judgment = "PASS (< %.0f%% frozen-O3 threshold)" % (O3_REL_DIFF_OK * 100)
        else:
            judgment = ("WARNING: o3_rel_diff=%r exceeds the %.0f%% frozen-O3 "
                        "threshold" % (o3_rel_diff, O3_REL_DIFF_OK * 100))
        rows.append(('o3_rel_diff', o3_rel_diff, judgment))

    if 'solar_max_abs_diff_wm2' in sanity_checks:
        solar_diff = sanity_checks['solar_max_abs_diff_wm2']
        solar_is_varying = bool(varying) and 'solar' in varying
        if solar_is_varying:
            judgment = ("(solar is the varying forcing for experiment=%r -- a "
                        "large diff is expected here, not checked)" % experiment)
        elif solar_diff is not None and solar_diff < SOLAR_DIFF_OK_WM2:
            judgment = "PASS (< %.1f W/m^2 frozen-solar threshold)" % SOLAR_DIFF_OK_WM2
        else:
            judgment = ("WARNING: solar_max_abs_diff_wm2=%r exceeds the %.1f "
                        "W/m^2 frozen-solar threshold" % (solar_diff, SOLAR_DIFF_OK_WM2))
        rows.append(('solar_max_abs_diff_wm2', solar_diff, judgment))

    # Any other sanity_checks keys we don't have a specific rule for yet --
    # report the raw value with no judgment rather than silently dropping it.
    for k, v in sanity_checks.items():
        if k not in ('co2_equal', 'o3_rel_diff', 'solar_max_abs_diff_wm2'):
            rows.append((k, v, "(no judgment rule for this key)"))

    return rows


def build_summary_text(case, cfg, provenance, additivity, sanity_rows, checksums):
    lines = []
    W = 78
    lines.append("=" * W)
    lines.append("pyCFRAM run summary: %s" % case)
    lines.append("=" * W)

    # ---- Header --------------------------------------------------------
    prov = provenance or {}
    src_cfg = cfg.get('source', {}) or {}
    lines.append("Case name         : %s" % cfg.get('case_name', case))
    lines.append("Description       : %s" % cfg.get('description', ''))
    lines.append("Model             : %s" % prov.get('model', src_cfg.get('model', 'N/A')))
    lines.append("Experiment        : %s" % prov.get('experiment', src_cfg.get('experiment', 'N/A')))
    lines.append("Variant           : %s" % prov.get('variant', src_cfg.get('variant', 'N/A')))
    lines.append("Grid label        : %s" % prov.get('grid_label', src_cfg.get('grid_label', 'N/A')))
    lines.append("Base years        : %s" % prov.get('base_years', src_cfg.get('base_years', 'N/A')))
    lines.append("Warm years        : %s" % prov.get('warm_years', src_cfg.get('warm_years', 'N/A')))
    lines.append("Target plev (hPa, TOA->sfc): %s"
                 % prov.get('target_plev_hpa_toa2sfc', cfg.get('grid', {}).get('pressure_levels', 'N/A')))
    if additivity.get('ok'):
        lines.append("Grid dimensions   : lat=%d lon=%d lev=%d"
                     % (additivity['nlat'], additivity['nlon'], additivity['nlev']))
    else:
        lines.append("Grid dimensions   : N/A (%s)" % additivity.get('reason', 'cfram_result.nc unreadable'))
    lines.append("Radiation engine  : %s" % cfg.get('radiation', {}).get('scheme', 'N/A'))
    if provenance is None:
        lines.append("")
        lines.append("NOTE: cases/%s/input/provenance.json not found -- this is expected "
                      "for non-DAMIP-sourced cases (ERA5, cesm2_cmip6-collab, idealized "
                      "climlab cases, ...) or a DAMIP case built before WP-M4.2/M4.4. "
                      "The process-activity table below is omitted." % case)

    # ---- Process activity table ----------------------------------------
    lines.append("")
    lines.append("-" * W)
    lines.append("Process activity (cases/%s/input/provenance.json)" % case)
    lines.append("-" * W)
    if provenance is not None and 'processes' in provenance:
        for proc, info in provenance['processes'].items():
            status = info.get('status', 'UNKNOWN')
            reason = info.get('reason', '')
            # Status vocabulary in practice includes ACTIVE/SKIPPED/
            # CLIMATOLOGY/ANALYTIC/FALLBACK/UNAVAILABLE (docs/plan_ph3.md
            # §6.3) plus MODEL (o3: "used the model's own o3, not
            # injected", data/cmip6_damip_source.py). All formatted the
            # same way -- we don't validate against the vocabulary, just
            # report exactly what build_states() recorded.
            lines.append("  %-10s %-13s %s" % (proc, status, reason))
    else:
        lines.append("  (not available -- see NOTE above)")

    # ---- Additivity residual --------------------------------------------
    lines.append("")
    lines.append("-" * W)
    lines.append("Additivity residual: dT_observed[sfc] - Sum(dT_X[sfc])")
    lines.append("-" * W)
    if additivity.get('ok'):
        lines.append("Surface-level convention: last index of the `lev` dimension "
                     "(cross-checked against scripts/plot_closure_profile.py's "
                     "`dT_obs[-1]` / `nlev_atm = len(lev)-1` and "
                     "scripts/compute_lr_attribution.py's `skin_idx` fallback "
                     "`lev.size - 1`).")
        lines.append("")
        lines.append("Terms summed (N=%d): %s"
                     % (len(additivity['used_terms']), ", ".join(additivity['used_terms']) or "(none)"))
        if additivity['excluded']:
            lines.append("Excluded as duplicate/derived aggregates (would double-count):")
            for name, reason in sorted(additivity['excluded'].items()):
                lines.append("  dT_%-20s %s" % (name, reason))
        if additivity['notes']:
            lines.append("Ambiguity notes (not silently resolved):")
            for note in additivity['notes']:
                lines.append("  - %s" % note)
        lines.append("")
        lines.append("Domain-mean residual  : %+.6f K" % additivity['resid_mean'])
        lines.append("Domain-max |residual| : %.6f K" % additivity['resid_max_abs'])
        lines.append("")
        lines.append(
            "This residual is EXPECTED CFRAM first-order-expansion "
            "non-linearity, NOT a bug: CFRAM's decomposition linearizes the "
            "Planck feedback matrix once per process, but the true radiative "
            "response to the *sum* of all processes' temperature changes is "
            "not exactly the sum of the individual processes' responses. "
            "See %s (§5 'M3 results', §7 'Known limitations') -- e.g. a ~22%% "
            "pointwise |residual| with < 0.5%% domain-mean bias was documented "
            "there for a comparable (cesm2_4xco2_official) full-decomposition "
            "case; this run's residual above is case-specific." % M3_DOC_REL)
    else:
        lines.append("N/A: %s" % additivity.get('reason'))

    # ---- Single-forcing consistency -------------------------------------
    lines.append("")
    lines.append("-" * W)
    lines.append("Single-forcing consistency (provenance.json sanity_checks)")
    lines.append("-" * W)
    if sanity_rows:
        for key, value, judgment in sanity_rows:
            lines.append("  %-24s = %-10s %s" % (key, value, judgment))
    else:
        lines.append("  (not available -- no provenance.json / no sanity_checks block)")

    # ---- Footer -----------------------------------------------------------
    lines.append("")
    lines.append("-" * W)
    lines.append("Reproduction")
    lines.append("-" * W)
    lines.append("python3 run_case.py %s --step build && python3 run_case.py %s --step run"
                 % (case, case))
    lines.append("")
    lines.append("Input file checksums (md5):")
    if checksums:
        for name, digest in checksums:
            lines.append("  %-24s %s" % (name, digest))
    else:
        lines.append("  (no *.nc files found under cases/%s/input/)" % case)
    lines.append("")
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', required=True, help='Case name (directory under cases/)')
    args = parser.parse_args()

    cfg = load_case(args.case)
    case_dir = cfg['_case_dir']
    output_dir = cfg['_output_dir']
    input_dir = os.path.join(case_dir, 'input')

    prov_path = os.path.join(input_dir, 'provenance.json')
    provenance = None
    if os.path.exists(prov_path):
        with open(prov_path) as f:
            provenance = json.load(f)

    result_nc = os.path.join(output_dir, 'cfram_result.nc')
    if not os.path.exists(result_nc):
        sys.exit("ERROR: %s not found -- run `python3 run_case.py %s --step run` first."
                  % (result_nc, args.case))
    additivity = compute_additivity(result_nc)

    experiment = (provenance or {}).get('experiment')
    damip_experiments_cfg = load_damip_experiments_cfg()
    sanity_rows = []
    if provenance is not None and 'sanity_checks' in provenance:
        sanity_rows = judge_sanity_checks(
            provenance['sanity_checks'], experiment, damip_experiments_cfg)

    checksums = []
    for nc_path in sorted(glob.glob(os.path.join(input_dir, '*.nc'))):
        checksums.append((os.path.basename(nc_path), md5_of_file(nc_path)))

    text = build_summary_text(args.case, cfg, provenance, additivity, sanity_rows, checksums)

    out_path = os.path.join(output_dir, '%s.summary.txt' % args.case)
    os.makedirs(output_dir, exist_ok=True)
    with open(out_path, 'w') as f:
        f.write(text)
    print("Wrote %s" % out_path)


if __name__ == '__main__':
    main()
