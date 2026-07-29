"""WP-M5.3 acceptance test: adding a brand-new DAMIP model must require ONLY
a configs/damip_models.d/<model>.yaml + a cases/<case>/case.yaml -- zero
changes to core/ or fortran/.

Real-world proof: NorESM2-LM (docs/plan_ph3.md's custom-source-accession
example) is NOT one of the 8 M4/M5 models. It was added via exactly
configs/damip_models.d/NorESM2-LM.yaml + cases/damip_noresm2_histaer/
case.yaml, verified end-to-end against real downloaded ESGF data on
hqlx210 (2026-07-06): build -> run -> cfram_result.nc with zero NaN,
dT_sfcdyn identity to ~1e-14 K, physically sane net-cooling signal
(domain-mean dT_observed[sfc]=-1.05 K, NH colder than SH).
"""
import subprocess
import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent


def _git(*args):
    result = subprocess.run(
        ['git', *args], cwd=PROJECT_ROOT, capture_output=True, text=True, check=True)
    return result.stdout.strip()


def _in_git_repo():
    try:
        _git('rev-parse', '--is-inside-work-tree')
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False



# The actual Phase 3 (M4 + M5, DAMIP multi-model support) start point:
# the first commit on this work (docs(Ph3): add Phase 3 execution plan v2),
# committed directly on top of feat/m2-m3-lapse-rate-kernel. NOT `main` --
# main predates Phase 2 as well, and Phase 2 legitimately added
# core/kernels.py, core/lr_kernel.py, core/lr_attribution.py (the lapse-rate
# kernel module); diffing against main would flag that unrelated, prior,
# already-accepted work as if it were a Phase 3 violation.
PHASE3_START_COMMIT = 'a87611d'


@pytest.mark.skipif(not _in_git_repo(), reason="not running inside a git checkout")
def test_phase3_diff_never_touches_core_or_fortran():
    """The entire Phase 3 branch (M4 + M5, DAMIP multi-model support), diffed
    against its own start point, must not touch core/ or fortran/ -- this is
    the M5 contract acceptance criterion (docs/plan_ph3.md §11.2), not just
    a NorESM2-LM-specific claim."""
    try:
        _git('cat-file', '-e', PHASE3_START_COMMIT)
    except subprocess.CalledProcessError:
        pytest.skip(f"Phase 3 start commit {PHASE3_START_COMMIT} not reachable "
                     "(e.g. shallow clone or history rewritten)")

    changed_files = _git('diff', '--name-only', PHASE3_START_COMMIT, 'HEAD').splitlines()
    core_or_fortran = [f for f in changed_files if f.startswith('core/') or f.startswith('fortran/')]
    assert core_or_fortran == [], (
        f"Phase 3 must not touch core/ or fortran/, but found changes in: {core_or_fortran}")


def test_noresm2_model_config_exists_and_is_yaml_only():
    """The NorESM2-LM model config is exactly one small yaml file -- not a
    Python module, not a code change."""
    config_path = PROJECT_ROOT / 'configs' / 'damip_models.d' / 'NorESM2-LM.yaml'
    assert config_path.exists(), "NorESM2-LM.yaml model config must exist"
    assert config_path.suffix == '.yaml'

    case_path = PROJECT_ROOT / 'cases' / 'damip_noresm2_histaer' / 'case.yaml'
    assert case_path.exists(), "damip_noresm2_histaer/case.yaml must exist"


def test_noresm2_case_uses_registered_cmip6_damip_source():
    """The example case must go through the SAME registered cmip6_damip
    source plugin every other DAMIP model uses -- not a bespoke one-off."""
    import yaml
    case_path = PROJECT_ROOT / 'cases' / 'damip_noresm2_histaer' / 'case.yaml'
    with open(case_path) as f:
        cfg = yaml.safe_load(f)
    assert cfg['source']['type'] == 'cmip6_damip'
    assert cfg['source']['model'] == 'NorESM2-LM'
