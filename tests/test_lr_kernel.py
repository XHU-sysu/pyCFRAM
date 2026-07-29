"""Unit tests for core/lr_kernel.py -- pure numpy/scipy, no xesmf/climkern
dependency (see docs/plan.md WP-M2.2 Done criteria)."""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.lr_kernel import (
    tropopause_pa, layer_dp_troposphere, delta_R_lr, delta_R_planck,
    interp_to_kernel_plev, extract_and_apply,
)


def test_tropopause_pa_equator_and_pole():
    lat = np.array([0.0, 90.0, -90.0])
    tropo = tropopause_pa(lat)
    np.testing.assert_allclose(tropo[0], 1e4)   # 100 hPa at equator
    np.testing.assert_allclose(tropo[1], 3e4)   # 300 hPa at poles
    np.testing.assert_allclose(tropo[2], 3e4)


def test_layer_dp_sums_to_ps_minus_tropo():
    plev_pa = np.array([100000., 92500., 85000., 70000., 50000., 30000.,
                         10000., 5000., 1000.])  # descending, surface-first
    lat = np.linspace(-89, 89, 12)
    tropo = np.broadcast_to(tropopause_pa(lat)[:, None], (12, 8))
    ps = np.full((12, 8), 101300.0)
    dp = layer_dp_troposphere(plev_pa, ps, tropo)
    total = dp.sum(axis=0)
    np.testing.assert_allclose(total, ps - tropo, atol=1.0)


def test_layer_dp_ascending_plev_matches_descending():
    """Order-independence: reversing plev order should reindex dp, not
    change the physical layer thicknesses."""
    plev_desc = np.array([100000., 70000., 50000., 20000., 5000.])
    plev_asc = plev_desc[::-1]
    lat = np.array([0.0, 45.0])
    tropo = np.broadcast_to(tropopause_pa(lat)[:, None], (2, 3))
    ps = np.full((2, 3), 101300.0)
    dp_desc = layer_dp_troposphere(plev_desc, ps, tropo)
    dp_asc = layer_dp_troposphere(plev_asc, ps, tropo)
    np.testing.assert_allclose(dp_desc, dp_asc[::-1])


def test_uniform_warming_gives_zero_lr_feedback():
    """ΔT(p) ≡ ΔTS everywhere -> ΔR_LR ≡ 0 to machine precision."""
    nk, nlat, nlon = 10, 4, 5
    plev_pa = np.linspace(100000, 5000, nk)
    dTs = np.full((nlat, nlon), 2.5)
    dTa = np.broadcast_to(dTs[None, ...], (nk, nlat, nlon)).copy()
    K = np.random.RandomState(0).uniform(-2, 2, size=(nk, nlat, nlon))
    dp = np.random.RandomState(1).uniform(0, 5000, size=(nk, nlat, nlon))
    dR_lr = delta_R_lr(dTa, dTs, K, dp)
    np.testing.assert_allclose(dR_lr, 0.0, atol=1e-10)


def test_constant_kernel_analytic_match():
    """K == -1 everywhere -> ΔR_LR = -Σ(ΔT-ΔTS)*dp/1e4, hand-computed."""
    nk, nlat, nlon = 3, 2, 2
    dTa = np.array([[[3.0, 1.0], [2.0, 0.5]],
                     [[1.0, 0.0], [1.5, 0.5]],
                     [[0.0, -1.0], [0.5, -0.5]]])
    dTs = np.array([[1.0, 0.0], [1.0, 0.0]])
    K = -np.ones((nk, nlat, nlon))
    dp = np.full((nk, nlat, nlon), 10000.0)  # 100 hPa each layer
    dR_lr = delta_R_lr(dTa, dTs, K, dp)
    expected = -np.sum(dTa - dTs[None, ...], axis=0) * (10000.0 / 1e4)
    np.testing.assert_allclose(dR_lr, expected)


def test_planck_feedback_uniform_warming():
    """Uniform warming: ΔR_PL should equal (K_ts + Σ K_t*dp/1e4) * ΔTS."""
    nk, nlat, nlon = 4, 2, 2
    dTs = np.full((nlat, nlon), 3.0)
    dTa = np.broadcast_to(dTs[None, ...], (nk, nlat, nlon)).copy()
    K_lw_t = np.full((nk, nlat, nlon), -1.5)
    K_lw_ts = np.full((nlat, nlon), -2.0)
    dp = np.full((nk, nlat, nlon), 10000.0)
    dR_pl = delta_R_planck(dTa, dTs, K_lw_t, K_lw_ts, dp)
    expected = (K_lw_ts + np.sum(K_lw_t * dp / 1e4, axis=0)) * dTs
    np.testing.assert_allclose(dR_pl, expected)


def test_interp_masks_underground_as_nan_before_interp():
    """Levels below ps should be NaN and propagate into neighbouring
    interpolated kernel levels (docs/plan.md §2.3 易错点 3)."""
    plev_native_pa = np.array([100000., 92500., 85000., 70000., 50000.])
    dTa = np.array([[[5.0]], [[4.0]], [[3.0]], [[2.0]], [[1.0]]])  # (5,1,1)
    ps = np.array([[90000.0]])  # level 0 (1000hPa) is underground
    plev_kernel_pa = np.array([100000., 92500., 85000., 70000., 50000.])
    out = interp_to_kernel_plev(dTa, plev_native_pa, plev_kernel_pa, ps)
    # the underground level (100000 Pa > ps) must be NaN
    assert np.isnan(out[0, 0, 0])


def test_extract_and_apply_uniform_warming_zero_lr():
    """End-to-end smoke test: uniform warming through the full pipeline
    (regridded kernel not required -- use a tiny synthetic KernelSet)."""
    class _FakeKernelSet:
        plev_pa = np.array([100000., 70000., 50000., 30000., 10000.])

        def select(self, month='annual', sky='all-sky'):
            nk, nlat, nlon = 5, 3, 4
            K_lw_t = np.random.RandomState(2).uniform(-2, 0, (nk, nlat, nlon))
            K_lw_ts = np.random.RandomState(3).uniform(-2, 0, (nlat, nlon))
            return K_lw_t, K_lw_ts

    nlat, nlon = 3, 4
    lat = np.linspace(-60, 60, nlat)
    dT_skin = np.full((nlat, nlon), 4.0)
    plev_native_pa = np.array([100000., 85000., 70000., 50000., 30000., 10000.])
    dT_atm = np.broadcast_to(dT_skin[None, ...],
                              (plev_native_pa.size, nlat, nlon)).copy()
    ps = np.full((nlat, nlon), 101300.0)

    result = extract_and_apply(dT_atm, dT_skin, ps, plev_native_pa, lat,
                                _FakeKernelSet(), month='annual', sky='all-sky')
    np.testing.assert_allclose(result['dR_lr'], 0.0, atol=1e-8)
