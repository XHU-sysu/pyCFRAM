"""Tests for scripts/write_run_summary.py (docs/plan_ph3.md §6.3, WP-M4.4).

The real cfram_result.nc for the M4 DAMIP cases does not exist yet (the
Fortran engine hasn't been run on damip_ipsl_histaer in this worktree --
that's WP-M4.5, pending remote data + a remote CFRAM run). So these tests
exercise the writer against the existing Phase-2 smoke fixture
tests/data/smoke/cfram_result_mini.nc (dT_observed + dT_ocndyn/dT_lhflx/
dT_shflx/dT_sfcdyn/dT_q/dT_atmdyn, 8x8 lat/lon, 20 lev), paired with a
synthetic provenance.json modeled on the real one written by
data/cmip6_damip_source.py's build_states() (captured by running
`python3 run_case.py damip_ipsl_histaer --step build` against the
committed tests/data/damip_smoke/ipsl_mini/ fixture).
"""
import glob
import json
import os
import shutil
import subprocess
import sys

import numpy as np
import pytest
from netCDF4 import Dataset

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_ROOT)

from scripts import write_run_summary as wrs  # noqa: E402

SMOKE_DATA = os.path.join(PROJECT_ROOT, 'tests', 'data', 'smoke')
SMOKE_RESULT_NC = os.path.join(SMOKE_DATA, 'cfram_result_mini.nc')
SMOKE_SURF_NC = os.path.join(SMOKE_DATA, 'input', 'perturbed_surf_mini.nc')

CASE_NAME = '_write_run_summary_smoke'
CASE_DIR = os.path.join(PROJECT_ROOT, 'cases', CASE_NAME)

CASE_YAML = """
case_name: WriteRunSummarySmoke
description: "write_run_summary.py smoke test (synthetic fixture)"
grid:
  pressure_levels: [100, 300, 500, 700, 1000]
input:
  base_pres: input/base_pres.nc
radiation:
  scheme: rrtmg
"""


def _default_provenance(**overrides):
    """Modeled on the real provenance.json written for damip_ipsl_histaer
    (data/cmip6_damip_source.py build_states(), verified 2026-07-05 by
    running `run_case.py damip_ipsl_histaer --step build` against
    tests/data/damip_smoke/ipsl_mini/)."""
    prov = {
        "source_type": "cmip6_damip",
        "model": "IPSL-CM6A-LR",
        "experiment": "hist-aer",
        "variant": "r1i1p1f1",
        "grid_label": "gr",
        "raw_dir": "tests/data/damip_smoke/ipsl_mini",
        "base_years": [1850, 1851],
        "warm_years": [1854, 1855],
        "target_plev_hpa_toa2sfc": [100.0, 300.0, 500.0, 700.0, 1000.0],
        "processes": {
            "albedo": {"status": "ACTIVE", "reason": "rsus/rsds both present -> compute_albedo"},
            "solar": {"status": "ACTIVE", "reason": "model rsdt present"},
            "flux": {"status": "ACTIVE", "reason": "hfls/hfss both present"},
            "cloud": {"status": "ACTIVE", "reason": "cl/clw/cli all present and not all-NaN"},
            "o3": {"status": "MODEL", "reason": "model o3 (mol/mol) x 48/29 -> kg/kg"},
            "aerosol": {"status": "SKIPPED", "reason": "hist-aer Amon has no 3D aerosol mixing ratio"},
        },
        "sanity_checks": {
            "co2_equal": True,
            "o3_rel_diff": 0.0,
            "solar_max_abs_diff_wm2": 0.0,
        },
    }
    prov.update(overrides)
    return prov


@pytest.fixture
def smoke_case():
    """A real cases/<name>/ directory (write_run_summary.py loads
    case.yaml via core.config.load_case, which hardcodes the cases/ path,
    so -- like tests/test_e2e_smoke.py -- we stage a throwaway case dir
    rather than mocking the loader)."""
    os.makedirs(os.path.join(CASE_DIR, 'input'), exist_ok=True)
    os.makedirs(os.path.join(CASE_DIR, 'output'), exist_ok=True)
    with open(os.path.join(CASE_DIR, 'case.yaml'), 'w') as f:
        f.write(CASE_YAML)
    with open(os.path.join(CASE_DIR, 'input', 'provenance.json'), 'w') as f:
        json.dump(_default_provenance(), f, indent=2)
    # Two small input NCs so the md5-checksum footer has something to list.
    shutil.copy(SMOKE_SURF_NC, os.path.join(CASE_DIR, 'input', 'base_pres.nc'))
    shutil.copy(SMOKE_SURF_NC, os.path.join(CASE_DIR, 'input', 'perturbed_pres.nc'))
    shutil.copy(SMOKE_RESULT_NC, os.path.join(CASE_DIR, 'output', 'cfram_result.nc'))

    yield CASE_DIR

    shutil.rmtree(CASE_DIR, ignore_errors=True)


def _run_summary_script(case_name=CASE_NAME):
    result = subprocess.run(
        [sys.executable, os.path.join(PROJECT_ROOT, 'scripts', 'write_run_summary.py'),
         '--case', case_name],
        cwd=PROJECT_ROOT, capture_output=True, text=True)
    return result


# ---------------------------------------------------------------------------
# 1. Header formatting
# ---------------------------------------------------------------------------

def test_header_fields_present(smoke_case):
    result = _run_summary_script()
    assert result.returncode == 0, "stdout:\n%s\nstderr:\n%s" % (result.stdout, result.stderr)

    summary_path = os.path.join(smoke_case, 'output', '%s.summary.txt' % CASE_NAME)
    assert os.path.exists(summary_path)
    text = open(summary_path).read()

    assert 'pyCFRAM run summary: %s' % CASE_NAME in text
    assert 'Model             : IPSL-CM6A-LR' in text
    assert 'Experiment        : hist-aer' in text
    assert 'Variant           : r1i1p1f1' in text
    assert 'Grid label        : gr' in text
    assert 'Base years        : [1850, 1851]' in text
    assert 'Warm years        : [1854, 1855]' in text
    assert '[100.0, 300.0, 500.0, 700.0, 1000.0]' in text
    # Grid dims come from the NC itself (8x8 lat/lon, 20 lev), not provenance.
    assert 'Grid dimensions   : lat=8 lon=8 lev=20' in text
    assert 'Radiation engine  : rrtmg' in text


def test_header_degrades_gracefully_without_provenance(smoke_case):
    """Non-DAMIP cases (ERA5, cesm2_cmip6-collab, ...) have no
    provenance.json -- the writer must not crash the `run` step for them."""
    os.remove(os.path.join(smoke_case, 'input', 'provenance.json'))
    result = _run_summary_script()
    assert result.returncode == 0, "stdout:\n%s\nstderr:\n%s" % (result.stdout, result.stderr)

    text = open(os.path.join(smoke_case, 'output', '%s.summary.txt' % CASE_NAME)).read()
    assert 'provenance.json not found' in text
    assert 'Model             : N/A' in text
    assert '(not available -- see NOTE above)' in text


# ---------------------------------------------------------------------------
# 2. Process-activity table: six-tag vocabulary + the MODEL o3 special case
# ---------------------------------------------------------------------------

def test_process_table_six_tag_vocabulary_and_model_status(smoke_case):
    prov = _default_provenance(processes={
        "albedo":  {"status": "FALLBACK",     "reason": "rsus/rsds missing -> albedo=0.15 constant"},
        "solar":   {"status": "ANALYTIC",     "reason": "no rsdt -> analytic_solar(lat)"},
        "flux":    {"status": "UNAVAILABLE",  "reason": "hfls/hfss missing -> no sfcdyn family"},
        "cloud":   {"status": "CLIMATOLOGY",  "reason": "synthetic test value"},
        "o3":      {"status": "MODEL",        "reason": "model o3 (mol/mol) x 48/29 -> kg/kg"},
        "aerosol": {"status": "SKIPPED",      "reason": "hist-aer Amon has no 3D aerosol"},
    })
    with open(os.path.join(smoke_case, 'input', 'provenance.json'), 'w') as f:
        json.dump(prov, f)

    result = _run_summary_script()
    assert result.returncode == 0, "stdout:\n%s\nstderr:\n%s" % (result.stdout, result.stderr)

    text = open(os.path.join(smoke_case, 'output', '%s.summary.txt' % CASE_NAME)).read()
    for status in ('ACTIVE', 'SKIPPED', 'CLIMATOLOGY', 'ANALYTIC', 'FALLBACK', 'UNAVAILABLE'):
        # ACTIVE isn't in this particular provenance -- check the other five
        # plus MODEL are all rendered without the process misclassifying or
        # crashing on an unrecognized status string.
        if status == 'ACTIVE':
            continue
        assert status in text, "status %r missing from summary:\n%s" % (status, text)
    assert 'o3         MODEL' in text
    assert 'model o3 (mol/mol) x 48/29 -> kg/kg' in text


# ---------------------------------------------------------------------------
# 3. Additivity residual, cross-checked with an independent numpy calc
# ---------------------------------------------------------------------------

def _mask_fill(arr):
    arr = np.asarray(arr, dtype=np.float64)
    return np.where(np.abs(arr) > 900.0, np.nan, arr)


def test_additivity_residual_matches_independent_numpy_calc():
    """core computation, called directly (no subprocess) so we can compare
    against a hand-rolled numpy reference computed independently in this
    test, per the task's validation requirement."""
    result = wrs.compute_additivity(SMOKE_RESULT_NC)
    assert result['ok']

    # Independent reference: sfc = last lev index (confirmed against
    # scripts/plot_closure_profile.py's `dT_obs[-1]` /
    # `nlev_atm = len(lev) - 1` convention).
    with Dataset(SMOKE_RESULT_NC) as d:
        obs = _mask_fill(d.variables['dT_observed'][-1])
        q = _mask_fill(d.variables['dT_q'][-1])
        atmdyn = _mask_fill(d.variables['dT_atmdyn'][-1])
        ocndyn = _mask_fill(d.variables['dT_ocndyn'][-1])
        lhflx = _mask_fill(d.variables['dT_lhflx'][-1])
        shflx = _mask_fill(d.variables['dT_shflx'][-1])
        sfcdyn = _mask_fill(d.variables['dT_sfcdyn'][-1])

    # Sanity-check the premise this fixture is built to exercise: sfcdyn
    # really is (to machine precision) ocndyn+lhflx+shflx here, so summing
    # both sfcdyn AND its components would double count.
    np.testing.assert_allclose(sfcdyn, ocndyn + lhflx + shflx, atol=1e-9)

    expected_total = q + atmdyn + ocndyn + lhflx + shflx
    expected_residual = obs - expected_total

    assert result['used_terms'] == sorted(['q', 'atmdyn', 'ocndyn', 'lhflx', 'shflx'])
    assert 'sfcdyn' in result['excluded']
    assert result['resid_mean'] == pytest.approx(float(np.nanmean(expected_residual)), abs=1e-9)
    assert result['resid_max_abs'] == pytest.approx(
        float(np.nanmax(np.abs(expected_residual))), abs=1e-9)


def test_additivity_residual_section_in_summary_text(smoke_case):
    result = _run_summary_script()
    assert result.returncode == 0
    text = open(os.path.join(smoke_case, 'output', '%s.summary.txt' % CASE_NAME)).read()

    expected = wrs.compute_additivity(SMOKE_RESULT_NC)
    assert ("Domain-mean residual  : %+.6f K" % expected['resid_mean']) in text
    assert ("Domain-max |residual| : %.6f K" % expected['resid_max_abs']) in text
    assert 'docs/m3_methodology_comparison.md' in text
    assert 'not a bug' in text.lower() or 'NOT a bug' in text


# ---------------------------------------------------------------------------
# select_additive_terms unit coverage (double-counting rules in isolation)
# ---------------------------------------------------------------------------

def test_select_additive_terms_excludes_sfcdyn_when_components_present():
    names = ['ocndyn', 'lhflx', 'shflx', 'sfcdyn', 'q', 'atmdyn']
    terms, excluded, notes = wrs.select_additive_terms(names)
    assert 'sfcdyn' not in terms
    assert set(terms) == {'ocndyn', 'lhflx', 'shflx', 'q', 'atmdyn'}
    assert 'sfcdyn' in excluded
    assert not notes


def test_select_additive_terms_keeps_sfcdyn_when_components_incomplete():
    # Only ocndyn present, no lhflx/shflx -- can't confirm sfcdyn is a pure
    # duplicate, so it must be kept (with an ambiguity note), not dropped.
    names = ['ocndyn', 'sfcdyn', 'q']
    terms, excluded, notes = wrs.select_additive_terms(names)
    assert 'sfcdyn' in terms
    assert 'sfcdyn' not in excluded
    assert any('sfcdyn' in n for n in notes)


def test_select_additive_terms_excludes_dry_warm_full():
    names = ['q', 'atmdyn', 'sfcdyn', 'ocndyn', 'lhflx', 'shflx', 'dry', 'warm', 'full']
    terms, excluded, notes = wrs.select_additive_terms(names)
    assert 'dry' not in terms and 'dry' in excluded
    assert 'warm' not in terms and 'warm' in excluded
    assert 'full' not in terms and 'full' in excluded


def test_select_additive_terms_cloud_lw_sw_vs_bulk():
    names = ['cloud', 'cloud_lw', 'cloud_sw', 'q']
    terms, excluded, notes = wrs.select_additive_terms(names)
    assert 'cloud' in terms
    assert 'cloud_lw' not in terms and 'cloud_sw' not in terms
    assert 'cloud_lw+cloud_sw' in excluded


def test_select_additive_terms_aerosol_species_vs_bulk():
    names = ['aerosol', 'bc', 'ocphi', 'ocpho', 'sulf', 'ss', 'dust', 'q']
    terms, excluded, notes = wrs.select_additive_terms(names)
    assert 'aerosol' in terms
    for s in ('bc', 'ocphi', 'ocpho', 'sulf', 'ss', 'dust'):
        assert s not in terms
    assert not notes


def test_select_additive_terms_aerosol_species_without_bulk_notes_ambiguity():
    names = ['bc', 'sulf', 'q']  # no bulk 'aerosol' present
    terms, excluded, notes = wrs.select_additive_terms(names)
    assert {'bc', 'sulf'} <= set(terms)
    assert any('aerosol' in n for n in notes)


# ---------------------------------------------------------------------------
# 4. Single-forcing sanity_checks pass/warn judgment
# ---------------------------------------------------------------------------

_HIST_AER_CFG = {"hist-aer": {"varying": ["aerosol"], "co2": "fixed_1850"}}
_HIST_NAT_CFG = {"hist-nat": {"varying": ["solar", "volcanic"], "co2": "fixed_1850"}}


def test_sanity_checks_pass_case():
    sanity = {"co2_equal": True, "o3_rel_diff": 0.0, "solar_max_abs_diff_wm2": 0.0}
    rows = wrs.judge_sanity_checks(sanity, "hist-aer", _HIST_AER_CFG)
    judgments = {k: j for k, v, j in rows}
    assert judgments['co2_equal'].startswith('PASS')
    assert judgments['o3_rel_diff'].startswith('PASS')
    assert judgments['solar_max_abs_diff_wm2'].startswith('PASS')


def test_sanity_checks_warn_case():
    """A hist-aer run where CO2 accidentally differs between base/warm (a
    real misconfiguration this check exists to catch) must WARN, not
    silently pass and must not raise."""
    sanity = {"co2_equal": False, "o3_rel_diff": 0.05, "solar_max_abs_diff_wm2": 5.0}
    rows = wrs.judge_sanity_checks(sanity, "hist-aer", _HIST_AER_CFG)
    judgments = {k: j for k, v, j in rows}
    assert judgments['co2_equal'].startswith('WARNING')
    assert judgments['o3_rel_diff'].startswith('WARNING')
    assert judgments['solar_max_abs_diff_wm2'].startswith('WARNING')


def test_sanity_checks_does_not_warn_on_the_experiments_own_varying_forcing():
    """hist-nat's whole point is a varying solar forcing -- a large solar
    diff there is correct, not a bug, and must not be flagged WARNING."""
    sanity = {"solar_max_abs_diff_wm2": 12.3}
    rows = wrs.judge_sanity_checks(sanity, "hist-nat", _HIST_NAT_CFG)
    judgments = {k: j for k, v, j in rows}
    assert 'WARNING' not in judgments['solar_max_abs_diff_wm2']


def test_sanity_checks_unknown_experiment_reports_without_crashing():
    sanity = {"co2_equal": True}
    rows = wrs.judge_sanity_checks(sanity, "some-unlisted-experiment", {})
    judgments = {k: j for k, v, j in rows}
    assert 'no CO2 expectation' in judgments['co2_equal']


def test_sanity_checks_section_in_summary_text_warn_case(smoke_case):
    prov = _default_provenance(sanity_checks={
        "co2_equal": False, "o3_rel_diff": 0.05, "solar_max_abs_diff_wm2": 5.0})
    with open(os.path.join(smoke_case, 'input', 'provenance.json'), 'w') as f:
        json.dump(prov, f)

    result = _run_summary_script()
    assert result.returncode == 0
    text = open(os.path.join(smoke_case, 'output', '%s.summary.txt' % CASE_NAME)).read()
    assert 'WARNING' in text


# ---------------------------------------------------------------------------
# 5. md5 footer + reproduction command
# ---------------------------------------------------------------------------

def test_footer_has_reproduction_command_and_checksums(smoke_case):
    result = _run_summary_script()
    assert result.returncode == 0
    text = open(os.path.join(smoke_case, 'output', '%s.summary.txt' % CASE_NAME)).read()

    assert ("python3 run_case.py %s --step build && python3 run_case.py %s --step run"
            % (CASE_NAME, CASE_NAME)) in text

    expected_files = sorted(
        os.path.basename(p) for p in glob.glob(os.path.join(smoke_case, 'input', '*.nc')))
    assert expected_files == ['base_pres.nc', 'perturbed_pres.nc']
    for fname in expected_files:
        digest = wrs.md5_of_file(os.path.join(smoke_case, 'input', fname))
        assert digest in text


def test_missing_cfram_result_exits_nonzero(smoke_case):
    os.remove(os.path.join(smoke_case, 'output', 'cfram_result.nc'))
    result = _run_summary_script()
    assert result.returncode != 0
    assert 'cfram_result.nc' in (result.stderr + result.stdout)
