"""Unit tests for core/lr_attribution.py (M3 route i)."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.lr_attribution import attribute_lr, additivity_residual, DEFAULT_TERMS


class _FakeKernelSet:
    plev_pa = np.array([100000., 70000., 50000., 30000., 10000.])

    def select(self, month='annual', sky='all-sky'):
        nk, nlat, nlon = 5, 3, 4
        K_lw_t = np.random.RandomState(7).uniform(-2, 0, (nk, nlat, nlon))
        K_lw_ts = np.random.RandomState(8).uniform(-2, 0, (nlat, nlon))
        return K_lw_t, K_lw_ts


def _make_dT_terms(nterm, nlat, nlon, plev_native_pa, seed=0):
    rs = np.random.RandomState(seed)
    terms = {}
    for i in range(nterm):
        atm = rs.uniform(-2, 2, (plev_native_pa.size, nlat, nlon))
        skin = rs.uniform(-2, 2, (nlat, nlon))
        terms['term%d' % i] = (atm, skin)
    return terms


def test_attribute_lr_returns_all_terms():
    nlat, nlon = 3, 4
    lat = np.linspace(-60, 60, nlat)
    plev_native_pa = np.array([100000., 85000., 70000., 50000., 30000., 10000.])
    ps = np.full((nlat, nlon), 101300.0)
    dT_terms = _make_dT_terms(3, nlat, nlon, plev_native_pa)

    out = attribute_lr(dT_terms, ps, plev_native_pa, lat, _FakeKernelSet())
    assert set(out.keys()) == set(dT_terms.keys())
    for term, result in out.items():
        assert 'dR_lr' in result and 'dR_pl' in result
        assert result['dR_lr'].shape == (nlat, nlon)


def test_additivity_residual_zero_when_terms_sum_to_total():
    nlat, nlon = 2, 2
    total = np.array([[5.0, 3.0], [1.0, -2.0]])
    by_term = {
        'a': np.array([[2.0, 1.0], [0.5, -1.0]]),
        'b': np.array([[3.0, 2.0], [0.5, -1.0]]),
    }
    resid = additivity_residual(total, by_term)
    np.testing.assert_allclose(resid, 0.0, atol=1e-12)


def test_additivity_residual_nonzero_reports_gap():
    nlat, nlon = 2, 2
    total = np.array([[10.0, 0.0], [0.0, 0.0]])
    by_term = {'a': np.array([[3.0, 0.0], [0.0, 0.0]])}
    resid = additivity_residual(total, by_term)
    np.testing.assert_allclose(resid, [[7.0, 0.0], [0.0, 0.0]])


def test_default_terms_nonempty():
    assert len(DEFAULT_TERMS) > 0
    assert 'ts' not in DEFAULT_TERMS  # ts has no vertical structure to decompose
