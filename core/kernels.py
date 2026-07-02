"""Radiative kernel I/O + regridding, dependency-free (numpy/scipy/netCDF4 only).

Reads ClimKern-format kernel files (``TOA_<name>_Kerns.nc``, e.g. the
Kramer/CloudSat and GFDL temperature kernels used for M2's lapse-rate
comparison) and regrids them to an arbitrary target lat/lon grid using a
self-contained bilinear interpolator. Deliberately does **not** import
xesmf/climkern -- this module must run on any machine (see
docs/plan.md §3.3 "原则"). The xesmf-based cross-check lives in
scripts/validate_lr_vs_climkern.py only.
"""
import numpy as np
from netCDF4 import Dataset
from scipy.interpolate import RegularGridInterpolator

# Variable name candidates per "kind", in the order ClimKern exposes them.
_TA_VARS = {'all-sky': 'lw_t', 'clear-sky': 'lwclr_t'}
_TS_VARS = {'all-sky': 'lw_ts', 'clear-sky': 'lwclr_ts'}


def _find_dim(nc, candidates):
    for c in candidates:
        if c in nc.dimensions:
            return c
    raise KeyError("None of %s found in dims %s" % (candidates, list(nc.dimensions)))


def _bilinear_regrid_2d(lat_src, lon_src, data_src, lat_tgt, lon_tgt):
    """Regrid one (lat_src, lon_src) field to (lat_tgt, lon_tgt).

    Bilinear, periodic in longitude (source is always the global kernel
    grid, wrapped by padding one column each side). Latitude targets
    beyond the source range are clipped to the nearest source boundary
    before interpolating -- an approximation of ClimKern's
    ``extrap_method="nearest_s2d"`` (see docs/plan.md WP-M2.1 回退).
    """
    lat_src = np.asarray(lat_src, dtype=np.float64)
    lon_src = np.asarray(lon_src, dtype=np.float64)
    data_src = np.asarray(data_src, dtype=np.float64)

    if lat_src[0] > lat_src[-1]:
        lat_src = lat_src[::-1]
        data_src = data_src[::-1, :]

    dlon = lon_src[1] - lon_src[0]
    lon_pad = np.concatenate(([lon_src[0] - dlon], lon_src, [lon_src[-1] + dlon]))
    data_pad = np.concatenate(
        [data_src[:, -1:], data_src, data_src[:, :1]], axis=1)

    interp = RegularGridInterpolator(
        (lat_src, lon_pad), data_pad, method='linear',
        bounds_error=False, fill_value=None)

    lat_tgt = np.asarray(lat_tgt, dtype=np.float64)
    lon_tgt = np.asarray(lon_tgt, dtype=np.float64)
    lat_clip = np.clip(lat_tgt, lat_src.min(), lat_src.max())
    lon_wrapped = np.mod(lon_tgt - lon_src[0], 360.0) + lon_src[0]

    LA, LO = np.meshgrid(lat_clip, lon_wrapped, indexing='ij')
    pts = np.stack([LA.ravel(), LO.ravel()], axis=-1)
    return interp(pts).reshape(lat_tgt.size, lon_tgt.size)


class KernelSet:
    """A single radiative kernel file (e.g. TOA_CloudSat_Kerns.nc), lazily
    regridded to a target lat/lon grid.

    Attributes (post-construction, on the *kernel's native* grid)
    ----------
    name : str
    plev_pa : (nk,) Pa
    lat, lon : (nlat,), (nlon,) degrees, native kernel grid
    nmonth : int (usually 12)
    vars : dict[str, ndarray]  # (time, plev, lat, lon) or (time, lat, lon)
    """

    def __init__(self, nc_path, name=None):
        self.path = nc_path
        self.name = name or nc_path
        with Dataset(nc_path) as nc:
            time_dim = _find_dim(nc, ['time', 'month'])
            lat_dim = _find_dim(nc, ['lat', 'latitude'])
            lon_dim = _find_dim(nc, ['lon', 'longitude'])
            plev_dim = None
            for c in ['plev', 'lev', 'player', 'level']:
                if c in nc.dimensions:
                    plev_dim = c
                    break

            self.lat = np.ma.filled(nc.variables[lat_dim][:].astype(np.float64), np.nan)
            self.lon = np.ma.filled(nc.variables[lon_dim][:].astype(np.float64), np.nan)
            self.nmonth = nc.dimensions[time_dim].size

            if plev_dim is not None:
                plev = np.ma.filled(nc.variables[plev_dim][:].astype(np.float64), np.nan)
                units = getattr(nc.variables[plev_dim], 'units', 'Pa')
                if str(units).lower() in ('hpa', 'mb', 'millibars'):
                    plev = plev * 100.0
                self.plev_pa = plev
            else:
                self.plev_pa = None

            self.vars = {}
            for vname in nc.variables:
                v = nc.variables[vname]
                if vname in (time_dim, lat_dim, lon_dim, plev_dim):
                    continue
                if lat_dim in v.dimensions and lon_dim in v.dimensions:
                    # netCDF4 returns masked arrays; np.array() on a masked
                    # array silently substitutes the raw fill value (~1e36)
                    # for masked (e.g. underground) points instead of NaN --
                    # must go through np.ma.filled() to get real NaNs.
                    raw = v[:]
                    self.vars[vname] = np.ma.filled(raw.astype(np.float64), np.nan)

    def regrid_to(self, lat_tgt, lon_tgt):
        """Return a new KernelSet-like object with all variables regridded
        to (lat_tgt, lon_tgt). plev_pa / nmonth are unchanged (regridding
        is horizontal-only)."""
        out = KernelSet.__new__(KernelSet)
        out.path = self.path
        out.name = self.name
        out.lat = np.asarray(lat_tgt, dtype=np.float64)
        out.lon = np.asarray(lon_tgt, dtype=np.float64)
        out.nmonth = self.nmonth
        out.plev_pa = self.plev_pa
        out.vars = {}
        for vname, arr in self.vars.items():
            # arr shape: (time, plev, lat, lon) or (time, lat, lon)
            leading_shape = arr.shape[:-2]
            flat = arr.reshape(-1, arr.shape[-2], arr.shape[-1])
            regridded = np.empty((flat.shape[0], out.lat.size, out.lon.size))
            for i in range(flat.shape[0]):
                regridded[i] = _bilinear_regrid_2d(
                    self.lat, self.lon, flat[i], out.lat, out.lon)
            out.vars[vname] = regridded.reshape(leading_shape + (out.lat.size, out.lon.size))
        return out

    def select(self, month='annual', sky='all-sky'):
        """Return (K_lw_t, K_lw_ts) on the current grid.

        month : 'annual' (12-month mean) or int 1-12
        sky   : 'all-sky' or 'clear-sky'

        Returns
        -------
        K_lw_t  : (nk, lat, lon)  W/m^2/K/100hPa
        K_lw_ts : (lat, lon)      W/m^2/K/100hPa
        """
        ta_name = _TA_VARS[sky]
        ts_name = _TS_VARS[sky]
        if ta_name not in self.vars or ts_name not in self.vars:
            raise KeyError(
                "Kernel %s missing %s/%s (have: %s)"
                % (self.name, ta_name, ts_name, sorted(self.vars)))

        ta = self.vars[ta_name]   # (time, plev, lat, lon)
        ts = self.vars[ts_name]   # (time, lat, lon)

        if month == 'annual':
            K_lw_t = np.nanmean(ta, axis=0)
            K_lw_ts = np.nanmean(ts, axis=0)
        else:
            idx = int(month) - 1
            K_lw_t = ta[idx]
            K_lw_ts = ts[idx]
        return K_lw_t, K_lw_ts
