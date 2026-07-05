"""CESM2 4xCO2 regression gold — WP-M4.1, docs/plan_ph3.md §8/§10 "门①".

This is the highest-priority gate in Phase 3: WP-M4.1 extracted the
model-agnostic CMIP6 machinery out of `data/cesm2_cmip6_source.py` into
`data/cmip6_common.py`. This test is the safety net proving that refactor
did not change CESM2's numerical output.

Two-layer design (docs/plan_ph3.md §8 WP-M4.1)
------------------------------------------------
1. **Local layer (this file)**: a synthetic, hand-built CESM2-style mini
   dataset (`tests/data/damip_smoke/cesm2_mini/`, ~4x4 lat/lon x a handful
   of levels x 3 years monthly, noleap calendar, `a/b/p0` hybrid
   coefficients, deliberately including CMIP6 fill-value (~1e20) cells to
   exercise the fill->NaN masking path). The golden values
   (`golden_output.npz`) were captured by running the PRE-REFACTOR
   `data/cesm2_cmip6_source.py` (`load_climo_pres` +
   `hybrid_to_plev_mass_conserving` + `compute_albedo`, i.e. exactly the
   pipeline `scripts/build_cesm2_official.py` runs) against this fixture,
   BEFORE any of the WP-M4.1 refactor edits were made. This test re-runs
   the SAME pipeline through the now-refactored (cmip6_common-backed) code
   and asserts bit-exact equality.
2. **Remote layer (NOT covered by this file — must be run separately by an
   operator with hqlx210 access)**: rerun `run_case.py cesm2_4xco2_official
   --step build` against the REAL CESM2 CMIP6 raw data (6.9+4.1+21 GB, only
   on hqlx210, see `.claude/persistent_context.md`) both before and after
   this refactor, and diff the resulting `base_pres.nc` / `perturbed_pres.nc`
   / etc. (md5 or numeric diff). **This has NOT been done as part of this
   change** — the local-layer pass here is necessary but not sufficient
   evidence that the real CESM2 pipeline is unaffected. See the final
   report of WP-M4.1 for explicit hand-off language.

If this test ever needs regenerating (e.g. the fixture changes), rerun the
pre-refactor `load_climo_pres`/`hybrid_to_plev_mass_conserving`/
`compute_albedo` pipeline against `tests/data/damip_smoke/cesm2_mini/` and
overwrite `golden_output.npz` — but only from a git ref that predates the
`cmip6_common` refactor, otherwise the "golden" values are no longer
independent of the code under test.
"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from data.cesm2_cmip6_source import (
    load_climo_pres, hybrid_to_plev_mass_conserving, compute_albedo,
)

FIXTURE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            'data', 'damip_smoke', 'cesm2_mini')
GOLDEN_PATH = os.path.join(FIXTURE_DIR, 'golden_output.npz')

# Matches the mini fixture's target plev grid (sfc->TOA, Pa) — see
# tests/data/damip_smoke/cesm2_mini generation notes below.
PLEV_MINI_PA_SFC2TOA = np.array([100000.0, 70000.0, 50000.0, 30000.0, 10000.0, 3000.0])
YEAR_START, YEAR_END = 1, 2


@pytest.fixture(scope='module')
def golden():
    if not os.path.exists(GOLDEN_PATH):
        pytest.skip('golden_output.npz missing — regenerate per this file\'s module docstring')
    return dict(np.load(GOLDEN_PATH))


@pytest.fixture(scope='module')
def refactored_output():
    """Run the (post-refactor) CESM2 pipeline exactly as
    scripts/build_cesm2_official.py does: load_climo_pres for both
    piControl and abrupt-4XCO2, mass-conserving hybrid->plev projection of
    cl/clw/cli, and compute_albedo."""
    base = load_climo_pres(FIXTURE_DIR, 'piControl', YEAR_START, YEAR_END)
    warm = load_climo_pres(FIXTURE_DIR, 'abrupt-4XCO2', YEAR_START, YEAR_END)

    a, b, p0 = base['hybrid_a'], base['hybrid_b'], base['hybrid_p0']
    for state in (base, warm):
        for var in ('cl', 'clw', 'cli'):
            field_hyb = state[var]
            ps_2d = state['ps']
            field_plev_top_down = hybrid_to_plev_mass_conserving(
                field_hyb, a, b, p0, ps_2d, PLEV_MINI_PA_SFC2TOA[::-1])
            state[var + '_plev'] = field_plev_top_down[::-1]

    base['albedo'] = compute_albedo(base['rsus'], base['rsds'])
    warm['albedo'] = compute_albedo(warm['rsus'], warm['rsds'])

    out = {}
    for prefix, state in (('base', base), ('warm', warm)):
        for k, v in state.items():
            out['%s__%s' % (prefix, k)] = np.asarray(v)
    return out


def test_golden_keys_match(golden, refactored_output):
    """No variable silently dropped or renamed by the refactor."""
    assert set(refactored_output.keys()) == set(golden.keys())


def test_golden_values_bit_exact(golden, refactored_output):
    """The refactor (cesm2_cmip6_source.py -> thin shim over
    cmip6_common.py) must not change any numeric output for the CESM2
    path — bit-exact, not just close, since the underlying arithmetic is
    literally relocated, not rewritten."""
    mismatches = []
    for key in sorted(golden.keys()):
        g = golden[key]
        n = refactored_output[key]
        if g.shape != n.shape:
            mismatches.append('%s: shape %s != %s' % (key, g.shape, n.shape))
            continue
        if not np.array_equal(g, n, equal_nan=True):
            if np.issubdtype(g.dtype, np.floating):
                max_diff = float(np.nanmax(np.abs(g - n)))
            else:
                max_diff = None
            mismatches.append('%s: max_abs_diff=%s' % (key, max_diff))
    assert not mismatches, 'Regression in cmip6_common refactor:\n' + '\n'.join(mismatches)


def test_fillvalue_masking_preserved(golden):
    """Sanity check on the fixture itself: the deliberate CMIP6 fill-value
    (~1e20) injections in the mini fixture must have produced the two
    documented behaviors (see data/cmip6_common.py:annual_climo_from_monthly
    docstring) in the GOLDEN (pre-refactor) output — partial-month fill
    averages away to a finite value, all-selected-months fill stays NaN.
    This anchors that the fixture is actually exercising the fill-masking
    code path, not just incidentally passing."""
    # base__ta: cell (k=0, j=0, i=0) had fill in 8/24 selected months -> finite
    assert np.isfinite(golden['base__ta'][0, 0, 0])
    # base__ta: cell (k=-1, j=-1, i=-1) had fill in ALL 24 selected months -> NaN
    assert np.isnan(golden['base__ta'][-1, -1, -1])
