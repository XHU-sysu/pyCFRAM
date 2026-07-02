"""Unit tests for core/kernels.py -- pure numpy/scipy/netCDF4, no xesmf
required for the base tests (see docs/plan.md WP-M2.1)."""
import os
import sys

import numpy as np
import pytest
from netCDF4 import Dataset

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.kernels import KernelSet, _bilinear_regrid_2d


@pytest.fixture
def synthetic_kernel_nc(tmp_path):
    """A tiny synthetic ClimKern-format kernel file: 2 months, 3 plev,
    global 10x8 grid, with a known analytic field so regrid/select can be
    checked exactly."""
    path = tmp_path / "TOA_Synthetic_Kerns.nc"
    nlat, nlon, nplev, ntime = 10, 8, 3, 2
    lat = np.linspace(-81, 81, nlat)
    lon = np.linspace(0, 315, nlon)  # 45-deg spacing, global periodic
    plev = np.array([850.0, 500.0, 200.0])  # hPa

    with Dataset(path, 'w') as nc:
        nc.createDimension('time', ntime)
        nc.createDimension('plev', nplev)
        nc.createDimension('lat', nlat)
        nc.createDimension('lon', nlon)
        nc.createVariable('lat', 'f8', ('lat',))[:] = lat
        nc.createVariable('lon', 'f8', ('lon',))[:] = lon
        v = nc.createVariable('plev', 'f8', ('plev',))
        v[:] = plev
        v.units = 'hPa'

        # lw_t: simple separable field lat + lon/100 + 10*plev_index + 100*time
        LA, LO = np.meshgrid(lat, lon, indexing='ij')
        lw_t = np.empty((ntime, nplev, nlat, nlon))
        for t in range(ntime):
            for k in range(nplev):
                lw_t[t, k] = -1.0 - 0.01 * (LA + LO) - k - 10 * t
        nc.createVariable('lw_t', 'f8', ('time', 'plev', 'lat', 'lon'))[:] = lw_t

        lw_ts = np.empty((ntime, nlat, nlon))
        for t in range(ntime):
            lw_ts[t] = -2.0 - 0.02 * (LA + LO) - 10 * t
        nc.createVariable('lw_ts', 'f8', ('time', 'lat', 'lon'))[:] = lw_ts

    return str(path), lat, lon, plev


def test_kernelset_reads_dims_and_vars(synthetic_kernel_nc):
    path, lat, lon, plev = synthetic_kernel_nc
    ks = KernelSet(path, name='Synthetic')
    np.testing.assert_allclose(ks.lat, lat)
    np.testing.assert_allclose(ks.lon, lon)
    np.testing.assert_allclose(ks.plev_pa, plev * 100.0)  # hPa -> Pa
    assert ks.nmonth == 2
    assert 'lw_t' in ks.vars and ks.vars['lw_t'].shape == (2, 3, 10, 8)
    assert 'lw_ts' in ks.vars and ks.vars['lw_ts'].shape == (2, 10, 8)


def test_kernelset_select_annual_mean(synthetic_kernel_nc):
    path, lat, lon, plev = synthetic_kernel_nc
    ks = KernelSet(path, name='Synthetic')
    K_lw_t, K_lw_ts = ks.select(month='annual', sky='all-sky')
    expected_lw_t = ks.vars['lw_t'].mean(axis=0)
    expected_lw_ts = ks.vars['lw_ts'].mean(axis=0)
    np.testing.assert_allclose(K_lw_t, expected_lw_t)
    np.testing.assert_allclose(K_lw_ts, expected_lw_ts)


def test_kernelset_select_single_month(synthetic_kernel_nc):
    path, lat, lon, plev = synthetic_kernel_nc
    ks = KernelSet(path, name='Synthetic')
    K_lw_t, _ = ks.select(month=1, sky='all-sky')
    np.testing.assert_allclose(K_lw_t, ks.vars['lw_t'][0])


def test_regrid_identity_on_same_grid(synthetic_kernel_nc):
    """Regridding onto the kernel's own native grid should reproduce the
    original field (bilinear interpolation at the source nodes is exact)."""
    path, lat, lon, plev = synthetic_kernel_nc
    ks = KernelSet(path, name='Synthetic')
    ks2 = ks.regrid_to(lat, lon)
    np.testing.assert_allclose(ks2.vars['lw_t'], ks.vars['lw_t'], atol=1e-8)


def test_regrid_to_finer_grid_matches_analytic_field(synthetic_kernel_nc):
    """The synthetic field is bilinear-exact (linear in lat, linear in
    lon), so regridding to a different grid should reproduce the analytic
    formula exactly (away from the poles, where clipping kicks in)."""
    path, lat, lon, plev = synthetic_kernel_nc
    ks = KernelSet(path, name='Synthetic')
    lat_tgt = np.linspace(-60, 60, 13)
    # Stay away from the 315->360(=0) wrap seam: the synthetic field is a
    # linear ramp in lon (not actually periodic), so periodic wrapping is
    # only exact strictly inside the source lon range.
    lon_tgt = np.linspace(0, 300, 16)
    ks2 = ks.regrid_to(lat_tgt, lon_tgt)

    LA, LO = np.meshgrid(lat_tgt, lon_tgt, indexing='ij')
    expected_ts_t0 = -2.0 - 0.02 * (LA + LO)
    np.testing.assert_allclose(ks2.vars['lw_ts'][0], expected_ts_t0, atol=0.05)


def test_bilinear_regrid_periodic_wrap():
    """Longitude wrap: interpolating near lon=0/360 boundary should not
    show a discontinuity for a smooth periodic field."""
    lat_src = np.linspace(-80, 80, 9)
    lon_src = np.linspace(0, 315, 8)  # 45-deg spacing
    LA, LO = np.meshgrid(lat_src, lon_src, indexing='ij')
    data = np.cos(np.deg2rad(LO))  # periodic in lon

    lon_tgt = np.array([358.0, 1.0, 2.0])  # straddles the wrap point
    lat_tgt = np.array([0.0])
    out = _bilinear_regrid_2d(lat_src, lon_src, data, lat_tgt, lon_tgt)
    expected = np.cos(np.deg2rad(lon_tgt))
    np.testing.assert_allclose(out[0], expected, atol=0.05)


def test_regrid_vs_xesmf():
    """Cross-check native bilinear regrid against xesmf (docs/plan.md
    WP-M2.1 Done: corr>0.999, max|Δ|<1% of kernel value range). Only runs
    inside the pycfram-kern env; skipped elsewhere."""
    xe = pytest.importorskip("xesmf")
    xr = pytest.importorskip("xarray")

    lat_src = np.linspace(-85, 85, 18)
    lon_src = np.linspace(0, 340, 18)
    LA, LO = np.meshgrid(lat_src, lon_src, indexing='ij')
    data = np.sin(np.deg2rad(LA)) * np.cos(np.deg2rad(LO)) + 0.1 * LA

    lat_tgt = np.linspace(-89, 89, 40)
    lon_tgt = np.linspace(0, 358, 60)

    native = _bilinear_regrid_2d(lat_src, lon_src, data, lat_tgt, lon_tgt)

    src_da = xr.DataArray(data, dims=('lat', 'lon'),
                           coords={'lat': lat_src, 'lon': lon_src})
    tgt_grid = xr.Dataset({'lat': (['lat'], lat_tgt), 'lon': (['lon'], lon_tgt)})
    regridder = xe.Regridder(src_da, tgt_grid, method='bilinear',
                              periodic=True, extrap_method='nearest_s2d')
    xesmf_out = regridder(src_da).values

    interior = slice(2, -2)  # exclude polar rows where clip/extrap differ
    corr = np.corrcoef(native[interior].ravel(), xesmf_out[interior].ravel())[0, 1]
    assert corr > 0.999
    value_range = data.max() - data.min()
    assert np.nanmax(np.abs(native[interior] - xesmf_out[interior])) < 0.01 * value_range
