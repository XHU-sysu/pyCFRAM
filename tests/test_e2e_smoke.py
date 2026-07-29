"""End-to-end smoke test (docs/plan.md WP-M2.7): one CLI command,
compute_lr_kernel.py, against a tiny 8x8 case + coarse global kernel
subset (~0.5MB each, see tests/data/smoke/).
"""
import os
import shutil
import subprocess
import sys

import numpy as np
import pytest
from netCDF4 import Dataset

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SMOKE_DATA = os.path.join(PROJECT_ROOT, 'tests', 'data', 'smoke')
SMOKE_CASE_NAME = '_smoke_lr'
SMOKE_CASE_DIR = os.path.join(PROJECT_ROOT, 'cases', SMOKE_CASE_NAME)

CASE_YAML = """
case_name: smoke_lr_test
description: "M2 kernel module smoke test (tiny 8x8 subset + coarse kernels)"
grid:
  pressure_levels: [1, 5, 10, 20, 30, 50, 70, 100, 150, 200, 250, 300,
                    400, 500, 600, 700, 850, 925, 1000]
input:
  perturbed_surf: input/perturbed_surf.nc
lapse_rate:
  kernels: [SmokeCloudSat, SmokeGFDL]
  kernel_months: annual
  sky: all-sky
"""

KERNEL_NAMES = ['SmokeCloudSat', 'SmokeGFDL']


@pytest.fixture
def smoke_case():
    os.makedirs(os.path.join(SMOKE_CASE_DIR, 'input'), exist_ok=True)
    os.makedirs(os.path.join(SMOKE_CASE_DIR, 'output'), exist_ok=True)
    with open(os.path.join(SMOKE_CASE_DIR, 'case.yaml'), 'w') as f:
        f.write(CASE_YAML)
    shutil.copy(os.path.join(SMOKE_DATA, 'cfram_result_mini.nc'),
                os.path.join(SMOKE_CASE_DIR, 'output', 'cfram_result.nc'))
    shutil.copy(os.path.join(SMOKE_DATA, 'input', 'perturbed_surf_mini.nc'),
                os.path.join(SMOKE_CASE_DIR, 'input', 'perturbed_surf.nc'))

    staged_kernel_dirs = []
    for kname in KERNEL_NAMES:
        dst_dir = os.path.join(PROJECT_ROOT, 'data', 'kernels', kname)
        os.makedirs(dst_dir, exist_ok=True)
        staged_kernel_dirs.append(dst_dir)
        shutil.copy(
            os.path.join(SMOKE_DATA, 'kernels', kname, 'TOA_%s_Kerns.nc' % kname),
            os.path.join(dst_dir, 'TOA_%s_Kerns.nc' % kname))

    yield SMOKE_CASE_DIR

    shutil.rmtree(SMOKE_CASE_DIR, ignore_errors=True)
    for d in staged_kernel_dirs:
        shutil.rmtree(d, ignore_errors=True)


def test_compute_lr_kernel_e2e(smoke_case):
    result = subprocess.run(
        [sys.executable, os.path.join(PROJECT_ROOT, 'scripts', 'compute_lr_kernel.py'),
         SMOKE_CASE_NAME],
        cwd=PROJECT_ROOT, capture_output=True, text=True)
    assert result.returncode == 0, "stdout:\n%s\nstderr:\n%s" % (result.stdout, result.stderr)

    out_nc = os.path.join(smoke_case, 'output', 'lr_kernel.nc')
    assert os.path.exists(out_nc)

    with Dataset(out_nc) as d:
        for kname in KERNEL_NAMES:
            for prefix in ('dR_lr_', 'dR_pl_'):
                arr = np.array(d.variables[prefix + kname][:])
                assert arr.shape == (8, 8)
                # sanity: no NaN-fill-value explosion (regression guard for
                # the masked-array bug fixed in core/kernels.py)
                assert np.all(np.abs(arr) < 1000.0), \
                    "%s%s has out-of-range values: %r" % (prefix, kname, arr)


def test_run_case_step_lr_dispatches_to_script(smoke_case):
    """run_case.py --step lr should dispatch to compute_lr_kernel.py (docs/plan.md §4)."""
    result = subprocess.run(
        [sys.executable, os.path.join(PROJECT_ROOT, 'run_case.py'), SMOKE_CASE_NAME,
         '--step', 'lr'],
        cwd=PROJECT_ROOT, capture_output=True, text=True)
    assert result.returncode == 0, "stdout:\n%s\nstderr:\n%s" % (result.stdout, result.stderr)
    assert os.path.exists(os.path.join(smoke_case, 'output', 'lr_kernel.nc'))
