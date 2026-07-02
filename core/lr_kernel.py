"""Radiative-kernel Lapse-Rate / Planck decomposition (native implementation).

Reimplements the algorithm behind ClimKern's ``calc_T_feedbacks`` (see
``docs/plan.md`` §1.3/§4.2 for the full derivation) using only numpy/scipy,
so it can run on any machine without the xesmf/esmpy dependency stack.

Sign/units convention matches ClimKern throughout:
    - pressures in Pa
    - kernels in W m^-2 K^-1 per 100 hPa
    - dp/1e4 converts a Pa layer thickness into "per 100 hPa" units
"""
import numpy as np
from scipy.interpolate import interp1d


def tropopause_pa(lat_1d):
    """Diagnostic tropopause pressure (Pa), following climkern.util.make_tropo.

    tropo = 3e4 - 2e4*cos(lat) -> 1e4 Pa (100 hPa) at the equator,
    3e4 Pa (300 hPa) at the poles.

    Parameters
    ----------
    lat_1d : (nlat,) array, degrees

    Returns
    -------
    (nlat,) array, Pa
    """
    lat_1d = np.asarray(lat_1d, dtype=np.float64)
    return 3e4 - 2e4 * np.cos(np.deg2rad(lat_1d))


def layer_dp_troposphere(plev_pa, ps_pa, tropo_pa):
    """Tropospheric layer thickness (Pa) on a fixed pressure grid.

    Reimplements climkern.util.get_dp(layer='troposphere'):
      - interfaces = midpoints between adjacent plev; bottom bound = ps,
        top bound = 0 (TOA)
      - interfaces are clipped to [tropo, ps] (this is what makes
        underground layers and layers above the tropopause automatically
        get dp=0, i.e. excluded from any layer sum)
      - dp = |interface_k - interface_{k+1}|

    Handles either pressure-axis direction (ascending or descending);
    output dp is returned in the same order as the input plev_pa.

    Parameters
    ----------
    plev_pa : (nk,) array, Pa, monotonic (either direction)
    ps_pa, tropo_pa : broadcastable arrays, Pa (e.g. (lat,lon) or scalar)

    Returns
    -------
    dp : (nk,) + broadcast(ps_pa, tropo_pa).shape, Pa
    """
    plev_pa = np.asarray(plev_pa, dtype=np.float64)
    nk = plev_pa.size
    ps_pa = np.asarray(ps_pa, dtype=np.float64)
    tropo_pa = np.asarray(tropo_pa, dtype=np.float64)
    out_shape = np.broadcast(ps_pa, tropo_pa).shape
    ps_b = np.broadcast_to(ps_pa, out_shape)
    tropo_b = np.broadcast_to(tropo_pa, out_shape)

    # Sort descending (surface-first) so bottom bound = ps, top bound = 0,
    # matching climkern's "if plev[0] > plev[-1]" branch. Unsort at the end.
    order = np.argsort(-plev_pa)
    plev_sorted = plev_pa[order]

    interfaces = np.empty((nk + 1,) + out_shape, dtype=np.float64)
    interfaces[0] = ps_b
    if nk > 1:
        mids = 0.5 * (plev_sorted[:-1] + plev_sorted[1:])
        interfaces[1:nk] = mids[:, None, None] if len(out_shape) == 2 else \
            mids.reshape((nk - 1,) + (1,) * len(out_shape))
    interfaces[nk] = 0.0

    lo = tropo_b[None, ...]
    hi = ps_b[None, ...]
    interfaces = np.clip(interfaces, lo, hi)

    dp_sorted = np.abs(interfaces[:-1] - interfaces[1:])  # (nk,) + out_shape

    dp = np.empty_like(dp_sorted)
    dp[order] = dp_sorted
    return dp


def delta_R_lr(dTa_klev, dTs, K_lw_t, dp):
    """Lapse-rate TOA LW radiative perturbation (W/m^2).

    ΔR_LR = Σ_k K_lw_t(k) * (ΔT(k) - ΔTS) * dp(k)/1e4

    NaNs in the (ΔT - ΔTS) term (e.g. from underground extrapolation) are
    zeroed at this product step -- NOT before -- matching ClimKern's
    ``.fillna(0)`` placement (see docs/plan.md §2.3 "易错点 3").

    Parameters
    ----------
    dTa_klev : (nk,lat,lon) K, ΔT interpolated onto the kernel's plev grid
    dTs      : (lat,lon) K, ΔTS (skin)
    K_lw_t   : (nk,lat,lon) W/m^2/K/100hPa
    dp       : (nk,lat,lon) Pa

    Returns
    -------
    (lat,lon) W/m^2
    """
    nonuniform = np.nan_to_num(dTa_klev - dTs[None, ...])
    return np.nansum(K_lw_t * nonuniform * dp / 1e4, axis=0)


def delta_R_planck(dTa_klev, dTs, K_lw_t, K_lw_ts, dp):
    """Planck (vertically-uniform warming) TOA LW radiative perturbation.

    ΔR_PL = K_lw_ts * ΔTS + Σ_k K_lw_t(k) * ΔTS * dp(k)/1e4

    Parameters
    ----------
    dTa_klev : (nk,lat,lon) -- only used for its NaN mask / shape (mirrors
               ClimKern's xr.broadcast(diff_ts, diff_ta) then .fillna(0))
    dTs      : (lat,lon) K
    K_lw_t   : (nk,lat,lon)
    K_lw_ts  : (lat,lon)
    dp       : (nk,lat,lon)

    Returns
    -------
    (lat,lon) W/m^2
    """
    dTs_klev = np.broadcast_to(dTs[None, ...], dTa_klev.shape).copy()
    dTs_klev = np.where(np.isnan(dTa_klev), 0.0, dTs_klev)
    return K_lw_ts * dTs + np.nansum(K_lw_t * dTs_klev * dp / 1e4, axis=0)


def interp_to_kernel_plev(dTa_native, plev_native_pa, plev_kernel_pa, ps_pa,
                           underground_nan=True):
    """Interpolate ΔT(p) from the native CFRAM plev grid to the kernel's
    plev grid, preserving NaN propagation the way ClimKern's
    ``interp_like(..., fill_value="extrapolate")`` does.

    Order of operations (docs/plan.md §2.3 易错点 3, must not be reordered):
      1. Mask native ΔT below the surface (plev_native > ps) as NaN.
      2. Linearly interpolate (with extrapolation) onto the kernel's plev
         grid. NaNs poison the interpolation interval(s) touching them
         (this is the intended "diffusion" — do not fillna before this
         step).
      3. Leave any remaining NaN in place; the caller (delta_R_lr /
         delta_R_planck) zeroes them at the final product step.

    Parameters
    ----------
    dTa_native : (nk_native, lat, lon) K
    plev_native_pa : (nk_native,) Pa
    plev_kernel_pa : (nk_kernel,) Pa
    ps_pa : (lat, lon) Pa
    underground_nan : bool
        If True (default), mask native levels below ps as NaN before
        interpolating (step 1). Set False if the input is already masked
        (e.g. CFRAM output already has -999 fill handled upstream).

    Returns
    -------
    (nk_kernel, lat, lon) K, on plev_kernel_pa's grid.
    """
    dTa_native = np.array(dTa_native, dtype=np.float64, copy=True)
    plev_native_pa = np.asarray(plev_native_pa, dtype=np.float64)
    plev_kernel_pa = np.asarray(plev_kernel_pa, dtype=np.float64)
    nlat, nlon = dTa_native.shape[1:]

    if underground_nan:
        underground = plev_native_pa[:, None, None] > ps_pa[None, :, :]
        dTa_native = np.where(underground, np.nan, dTa_native)

    # interp1d requires ascending x; sort native plev ascending.
    order_n = np.argsort(plev_native_pa)
    x = plev_native_pa[order_n]
    y = dTa_native[order_n]

    out = np.empty((plev_kernel_pa.size, nlat, nlon), dtype=np.float64)
    # Vectorise over (lat,lon) via reshape; interp1d over axis=0.
    y2 = y.reshape(x.size, nlat * nlon)
    # NaN handling: interp1d with NaNs in y produces NaN for any output
    # point whose bracketing interval touches a NaN sample -- exactly the
    # "propagate" semantics we want. bounds_error=False + fill_value=
    # "extrapolate" mirrors ClimKern's interp_like(...,'extrapolate').
    f = interp1d(x, y2, axis=0, kind='linear', bounds_error=False,
                 fill_value='extrapolate')
    out_flat = f(plev_kernel_pa)
    out = out_flat.reshape(plev_kernel_pa.size, nlat, nlon)
    return out


def extract_and_apply(dT_atm, dT_skin, ps_pa, plev_native_pa, lat_1d,
                       kernelset, month='annual', sky='all-sky'):
    """End-to-end: native ΔT(p)/ΔTS -> ΔR_LR, ΔR_PL for one kernel.

    Mirrors ClimKern's calc_T_feedbacks pipeline (docs/plan.md §1.3/§4.2):
      1. ΔT(p) masked underground -> NaN (interp_to_kernel_plev, step 1)
      2. interpolate to kernel's plev grid, NaN-propagating (step 2)
      3. dp on kernel's plev grid via layer_dp_troposphere
      4. ΔR_LR / ΔR_PL via delta_R_lr / delta_R_planck (NaN -> 0 only here)

    Parameters
    ----------
    dT_atm : (nk_native, lat, lon) K -- ΔT(p), native CFRAM grid
    dT_skin : (lat, lon) K -- ΔTS
    ps_pa : (lat, lon) Pa -- surface pressure (perturbed state, per ClimKern)
    plev_native_pa : (nk_native,) Pa -- native CFRAM pressure levels
    lat_1d : (lat,) degrees
    kernelset : core.kernels.KernelSet, already regridded to (lat,lon)
    month : 'annual' or int (1-12)
    sky : 'all-sky' or 'clear-sky'

    Returns
    -------
    dict with keys 'dR_lr', 'dR_pl' (lat,lon) W/m^2
    """
    K_lw_t, K_lw_ts = kernelset.select(month=month, sky=sky)
    plev_kernel_pa = kernelset.plev_pa

    dTa_klev = interp_to_kernel_plev(dT_atm, plev_native_pa, plev_kernel_pa,
                                      ps_pa)

    tropo_pa = np.broadcast_to(tropopause_pa(lat_1d)[:, None],
                                ps_pa.shape)
    dp = layer_dp_troposphere(plev_kernel_pa, ps_pa, tropo_pa)

    dR_lr = delta_R_lr(dTa_klev, dT_skin, K_lw_t, dp)
    dR_pl = delta_R_planck(dTa_klev, dT_skin, K_lw_t, K_lw_ts, dp)
    return {'dR_lr': dR_lr, 'dR_pl': dR_pl}
