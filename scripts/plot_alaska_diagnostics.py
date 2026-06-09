#!/usr/bin/env python3
"""Plot Alaska wildfire response diagnostics and regional means."""

import calendar
import os
import sys
from datetime import datetime

import matplotlib
matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.config import load_case
from scripts.plot_fig3_self import get_map_extent


REGIONS = [
    ("main_fire_region", "Main fire"),
    ("extended_region", "Extended"),
    ("arctic_response_region", "Arctic response"),
]


def clean_missing(values):
    """Convert legacy sentinel values to NaN before any arithmetic."""
    arr = np.asarray(values, dtype=np.float64)
    missing = np.isclose(arr, -999.0) | (np.abs(arr) > 1.0e30)
    return np.where(missing, np.nan, arr)


def coord_name(ds, candidates):
    for name in candidates:
        if name in ds.coords or name in ds.dims:
            return name
    raise KeyError(f"none of these coordinates found: {candidates}")


def data_dir(cfg):
    return cfg["source"]["data_dir"]


def month_range(start_mmdd, end_mmdd):
    start_month = int(start_mmdd.split("-")[0])
    end_month = int(end_mmdd.split("-")[0])
    if end_month >= start_month:
        return range(start_month, end_month + 1)
    return list(range(start_month, 13)) + list(range(1, end_month + 1))


def period_bounds(year, temporal):
    sm, sd = [int(x) for x in temporal["start_mmdd"].split("-")]
    em, ed = [int(x) for x in temporal["end_mmdd"].split("-")]
    start = np.datetime64(datetime(year, sm, sd, 0, 0, 0))
    end = np.datetime64(datetime(year, em, ed, 23, 59, 59))
    return start, end


def open_month_series(root, year, month, filename, variables):
    path = os.path.join(
        root, f"era5_sl_{year}{month:02d}", filename
    )
    ds = xr.open_dataset(path, engine="h5netcdf")
    keep = [v for v in variables if v in ds]
    if len(keep) != len(variables):
        missing = sorted(set(variables) - set(keep))
        ds.close()
        raise KeyError(f"{path}: missing variables {missing}")
    out = ds[keep].load()
    for name in keep:
        out[name] = xr.where(np.isclose(out[name], -999.0), np.nan, out[name])
    ds.close()
    return out


def select_period(da, year, temporal):
    start, end = period_bounds(year, temporal)
    return da.sel(time=slice(start, end))


def load_period_variable(root, year, temporal, filename, variables):
    pieces = []
    for month in month_range(temporal["start_mmdd"], temporal["end_mmdd"]):
        pieces.append(open_month_series(root, year, month, filename, variables))
    ds = xr.concat(pieces, dim="time")
    ds = ds.sortby("time")
    return select_period(ds, year, temporal)


def mean_t2m(root, year, temporal):
    ds = load_period_variable(
        root,
        year,
        temporal,
        "data_stream-oper_stepType-instant_diag.nc",
        ["t2m"],
    )
    out = ds["t2m"].mean("time")
    ds.close()
    return out


def accum_to_wm2(da):
    arr = clean_missing(da)
    nday = arr.shape[0] // 4
    if nday < 1:
        raise ValueError(f"not enough 6-hourly samples in {da.name}")
    arr = arr[:nday * 4].reshape((nday, 4) + arr.shape[1:])
    # Downloader stores per-hour accumulation equivalents; convert daily
    # total energy to daily-mean W m-2.
    return arr.sum(axis=1) * 6.0 / 86400.0


def mean_surface_sw_and_albedo(root, year, temporal):
    ds = load_period_variable(
        root,
        year,
        temporal,
        "data_stream-oper_stepType-accum.nc",
        ["ssrd", "ssr"],
    )
    ssrd = accum_to_wm2(ds["ssrd"])
    ssr = accum_to_wm2(ds["ssr"])
    with np.errstate(divide="ignore", invalid="ignore"):
        albedo = (ssrd - ssr) / ssrd
    albedo = np.where(np.isfinite(albedo), albedo, np.nan)
    albedo = np.clip(albedo, 0.0, 1.0)
    coords = {
        "latitude": ds["latitude"],
        "longitude": ds["longitude"],
    }
    ssr_mean = xr.DataArray(
        np.nanmean(ssr, axis=0),
        dims=("latitude", "longitude"),
        coords=coords,
        name="ssr",
    )
    alb_mean = xr.DataArray(
        np.nanmean(albedo, axis=0),
        dims=("latitude", "longitude"),
        coords=coords,
        name="albedo",
    )
    ds.close()
    return ssr_mean, alb_mean


def climatology(root, years, temporal, loader):
    fields = [loader(root, year, temporal) for year in years]
    return xr.concat(fields, dim="year").mean("year")


def load_era5_anomalies(cfg):
    root = data_dir(cfg)
    temporal = cfg["source"]["temporal"]
    y0, y1 = temporal["clim_years"]
    clim_years = list(range(int(y0), int(y1) + 1))
    event_year = int(temporal["event_year"])

    base_t2m = climatology(root, clim_years, temporal, mean_t2m)
    pert_t2m = mean_t2m(root, event_year, temporal)
    dt2m = pert_t2m - base_t2m

    def load_ssr(root_, year, temporal_):
        return mean_surface_sw_and_albedo(root_, year, temporal_)[0]

    def load_albedo(root_, year, temporal_):
        return mean_surface_sw_and_albedo(root_, year, temporal_)[1]

    base_ssr = climatology(root, clim_years, temporal, load_ssr)
    pert_ssr = load_ssr(root, event_year, temporal)
    dssr = pert_ssr - base_ssr

    base_alb = climatology(root, clim_years, temporal, load_albedo)
    pert_alb = load_albedo(root, event_year, temporal)
    dalb = pert_alb - base_alb
    return dt2m, dssr, dalb


def load_bc_column_anomaly(cfg):
    base = xr.open_dataset(cfg["_case_dir"] + "/input/base_pres.nc",
                           engine="h5netcdf")
    pert = xr.open_dataset(cfg["_case_dir"] + "/input/perturbed_pres.nc",
                           engine="h5netcdf")
    lat_name = coord_name(base, ["lat", "latitude"])
    lon_name = coord_name(base, ["lon", "longitude"])
    lev_name = coord_name(base, ["lev", "level", "pressure_level"])
    lev = np.asarray(base[lev_name].values, dtype=np.float64)
    bc_base_da = base["bc"]
    bc_pert_da = pert["bc"]
    extra_dims = [
        dim for dim in bc_base_da.dims
        if dim not in (lev_name, lat_name, lon_name)
    ]
    if extra_dims:
        bc_base_da = bc_base_da.mean(extra_dims)
        bc_pert_da = bc_pert_da.mean(extra_dims)
    bc_base_da = bc_base_da.transpose(lev_name, lat_name, lon_name)
    bc_pert_da = bc_pert_da.transpose(lev_name, lat_name, lon_name)
    bc_base = clean_missing(bc_base_da.values)
    bc_pert = clean_missing(bc_pert_da.values)

    if lev[0] < lev[-1]:
        lev = lev[::-1]
        bc_base = bc_base[::-1, ...]
        bc_pert = bc_pert[::-1, ...]
    p = lev * 100.0
    interior = 0.5 * (p[:-1] + p[1:])
    bottom = min(101325.0, p[0] + 0.5 * (p[0] - p[1]))
    top = max(0.0, p[-1] - 0.5 * (p[-2] - p[-1]))
    edges = np.concatenate([[bottom], interior, [top]])
    dp = np.abs(edges[:-1] - edges[1:])
    valid = np.isfinite(bc_pert) & np.isfinite(bc_base)
    burden_terms = np.where(
        valid, (bc_pert - bc_base) * dp[:, None, None] / 9.80665, np.nan,
    )
    burden = np.nansum(burden_terms, axis=0)
    burden = np.where(np.any(valid, axis=0), burden, np.nan)

    out = xr.DataArray(
        burden * 1.0e6,
        dims=(lat_name, lon_name),
        coords={lat_name: base[lat_name], lon_name: base[lon_name]},
        name="bc_column",
    )
    base.close()
    pert.close()
    return out


def shift_longitudes(lons, data):
    lons = np.asarray(lons)
    shifted = np.where(lons > 180.0, lons - 360.0, lons)
    order = np.argsort(shifted)
    return shifted[order], data[:, order]


def shift_extent(extent):
    lon = [x - 360.0 if x > 180.0 else x for x in extent[:2]]
    return lon + extent[2:]


def region_box(plot_cfg, name):
    region = plot_cfg.get(name)
    if not region:
        return None
    lon = [x - 360.0 if x > 180.0 else x for x in region["lon"]]
    return lon + region["lat"]


def add_box(ax, box, color, lw):
    if not box:
        return
    lon0, lon1, lat0, lat1 = box
    ax.plot(
        [lon0, lon1, lon1, lon0, lon0],
        [lat0, lat0, lat1, lat1, lat0],
        color=color,
        linewidth=lw,
        transform=ccrs.PlateCarree(),
    )


def symmetric_levels(field, n=13, floor=1.0):
    vmax = np.nanpercentile(np.abs(field), 98.0)
    vmax = max(float(vmax), floor)
    return np.linspace(-vmax, vmax, n)


def plot_map(ax, lons, lats, field, title, units, cmap, levels, extent, boxes,
             native=False):
    if native:
        plot_lons = np.asarray(lons)
        plot_field = field
        data_crs = ccrs.PlateCarree()
    else:
        plot_lons, plot_field = shift_longitudes(lons, field)
        data_crs = ccrs.PlateCarree()
    norm = mcolors.BoundaryNorm(levels, plt.get_cmap(cmap).N, clip=True)
    ax.set_extent(extent, crs=data_crs)
    cf = ax.pcolormesh(
        plot_lons,
        lats,
        plot_field,
        cmap=cmap,
        norm=norm,
        transform=data_crs,
        shading="nearest",
    )
    ax.coastlines(resolution="50m", linewidth=0.6)
    ax.add_feature(cfeature.BORDERS, linewidth=0.35, edgecolor="0.35")
    ax.add_feature(cfeature.LAKES, linewidth=0.25, edgecolor="0.5",
                   facecolor="none")
    add_box(ax, boxes["extended_region"], "0.25", 1.0)
    add_box(ax, boxes["arctic_response_region"], "tab:orange", 1.2)
    add_box(ax, boxes["main_fire_region"], "tab:blue", 1.6)
    ax.gridlines(linewidth=0.3, color="gray", alpha=0.5, linestyle="--")
    ax.set_title(title, fontsize=10, fontweight="bold")
    cb = plt.colorbar(cf, ax=ax, orientation="horizontal", pad=0.035,
                      fraction=0.045, extend="both")
    cb.set_label(units, fontsize=8)
    cb.ax.tick_params(labelsize=7)


def field_array(da):
    lat_name = coord_name(da, ["lat", "latitude"])
    lon_name = coord_name(da, ["lon", "longitude"])
    return (
        np.asarray(da[lat_name].values),
        np.asarray(da[lon_name].values),
        clean_missing(da.values),
    )


def area_mean(da, region):
    lat_name = coord_name(da, ["lat", "latitude"])
    lon_name = coord_name(da, ["lon", "longitude"])
    lats = np.asarray(da[lat_name].values)
    lons = np.asarray(da[lon_name].values)
    lon0, lon1 = region["lon"]
    lat0, lat1 = region["lat"]
    lon_mask = (lons >= lon0) & (lons <= lon1)
    lat_mask = (lats >= lat0) & (lats <= lat1)
    sub = clean_missing(da.values)[np.ix_(lat_mask, lon_mask)]
    weights = np.cos(np.deg2rad(lats[lat_mask]))[:, None]
    denom = np.nansum(np.where(np.isfinite(sub), weights, 0.0))
    return float(np.nansum(sub * weights) / denom) if denom > 0 else np.nan


def robust_levels(field, floor=1.0, symmetric=True):
    arr = clean_missing(field)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return np.linspace(-floor, floor, 13)
    p2, p98 = np.nanpercentile(arr, [2, 98])
    p5, p95 = np.nanpercentile(arr, [5, 95])
    if symmetric and p5 < 0 < p95:
        vmax = max(abs(float(p2)), abs(float(p98)), floor)
        return np.linspace(-vmax, vmax, 13)
    if abs(p98 - p2) < floor:
        mid = 0.5 * (float(p2) + float(p98))
        return np.linspace(mid - floor / 2.0, mid + floor / 2.0, 13)
    return np.linspace(float(p2), float(p98), 13)


def plot_maps(cfg, fields, suffix=""):
    plot_cfg = cfg.get("plot", {})
    native = suffix.endswith("strict")
    if native:
        boxes = {
            key: plot_cfg.get(key, {}).get("lon", [])
            + plot_cfg.get(key, {}).get("lat", [])
            if plot_cfg.get(key) else None
            for key, _ in REGIONS
        }
        raw_extent = get_map_extent(cfg)
        extent = raw_extent
        projection = ccrs.PlateCarree(central_longitude=210)
    else:
        boxes = {key: region_box(plot_cfg, key) for key, _ in REGIONS}
        extent = shift_extent(get_map_extent(cfg))
        projection = ccrs.LambertConformal(
            central_longitude=-150,
            central_latitude=66,
            standard_parallels=(55, 75),
        )
    fig, axes = plt.subplots(
        2, 2, figsize=(12.0, 8.8),
        subplot_kw={"projection": projection},
        squeeze=False,
    )
    meta = [
        ("t2m", "2 m Temperature Anomaly", "K", "RdBu_r", 0.25, False),
        ("bc", "BC Column Anomaly (AOD proxy)", "mg m$^{-2}$", "PuOr_r", 0.02, False),
        ("ssr", "Net Surface SW Anomaly", "W m$^{-2}$", "BrBG_r", 0.25, True),
        ("albedo", "Surface Albedo Anomaly", "fraction", "PiYG", 0.01, True),
    ]
    for ax, (key, title, units, cmap, floor, symmetric) in zip(axes.ravel(), meta):
        lats, lons, arr = field_array(fields[key])
        levels = robust_levels(arr, floor=floor, symmetric=symmetric)
        stats = np.nanpercentile(arr, [2, 5, 50, 95, 98])
        print(
            f"{key} levels {levels[0]:.4g}..{levels[-1]:.4g}; "
            f"p2/p5/p50/p95/p98={stats}",
            flush=True,
        )
        plot_map(
            ax, lons, lats, arr, title, units, cmap, levels, extent, boxes,
            native=native,
        )
    fig.suptitle(
        "Alaska Wildfire 2022: Temperature, Smoke and Surface Radiation Signals",
        fontsize=14,
        fontweight="bold",
        y=0.985,
    )
    fig.subplots_adjust(left=0.04, right=0.98, top=0.93, bottom=0.05,
                        wspace=0.12, hspace=0.18)
    out = os.path.join(
        cfg["_figures_dir"], f"alaska_response_diagnostics{suffix}.png"
    )
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out}")


def plot_bars(cfg, fields, suffix=""):
    plot_cfg = cfg.get("plot", {})
    labels = [label for _, label in REGIONS]
    metrics = [
        ("t2m", "T2m", "K", 1.0, "#d95f02"),
        ("bc", "BC column", "mg m$^{-2}$", 1.0, "#7570b3"),
        ("ssr", "Net SW", "W m$^{-2}$", 1.0, "#1b9e77"),
        ("albedo", "Albedo", "pct. points", 100.0, "#66a61e"),
    ]
    values = {
        key: [
            area_mean(fields[key], plot_cfg[region_key]) * scale
            for region_key, _ in REGIONS
        ]
        for key, _, _, scale, _ in metrics
    }

    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.0), squeeze=False)
    for ax, (key, title, units, _, color) in zip(axes.ravel(), metrics):
        vals = values[key]
        bars = ax.bar(labels, vals, color=color, alpha=0.85)
        ax.axhline(0, color="0.2", linewidth=0.8)
        ax.set_title(title, fontsize=11, fontweight="bold")
        ax.set_ylabel(units, fontsize=9)
        ax.tick_params(axis="x", labelrotation=12)
        ax.grid(axis="y", color="0.85", linewidth=0.6)
        for bar, val in zip(bars, vals):
            offset = 0.03 * (max(np.abs(vals)) if max(np.abs(vals)) > 0 else 1)
            va = "bottom" if val >= 0 else "top"
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                val + (offset if val >= 0 else -offset),
                f"{val:.2g}",
                ha="center",
                va=va,
                fontsize=8,
            )
    fig.suptitle("Regional Mean Anomalies", fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    out = os.path.join(
        cfg["_figures_dir"], f"alaska_regional_mean_bars{suffix}.png"
    )
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out}")


def main():
    case_name = sys.argv[1] if len(sys.argv) > 1 else "alaska_wf22"
    suffix = sys.argv[2] if len(sys.argv) > 2 else ""
    cfg = load_case(case_name)
    os.makedirs(cfg["_figures_dir"], exist_ok=True)

    print("Loading ERA5 T2m, SSR and albedo anomalies...")
    dt2m, dssr, dalb = load_era5_anomalies(cfg)
    print("Loading BC column anomaly from CFRAM aerosol input...")
    dbc = load_bc_column_anomaly(cfg)

    fields = {
        "t2m": dt2m,
        "ssr": dssr,
        "albedo": dalb,
        "bc": dbc,
    }
    plot_maps(cfg, fields, suffix=suffix)
    plot_bars(cfg, fields, suffix=suffix)


if __name__ == "__main__":
    main()
