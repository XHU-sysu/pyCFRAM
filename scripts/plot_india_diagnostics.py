#!/usr/bin/env python3
"""India-Bangladesh Apr 2023 wet-heat diagnostics.

Outputs to cases/india_wb23/figures:
  - wet_heat_anomaly_phases.png
  - moisture_transport_phases.png
  - regional_contrast_phases.png
"""

import os
import sys
import argparse
import numpy as np
import xarray as xr

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Rectangle

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.config import load_case, PROJECT_ROOT


PHASES = [
    ("pre", "Apr 1-16", list(range(0, 16)), "india_wb23_pre"),
    ("core", "Apr 17-20", [16, 17, 18, 19], "india_wb23_core"),
    ("post", "Apr 21-30", list(range(20, 30)), "india_wb23_post"),
]
YEARS = list(range(1991, 2021))
EVENT_YEAR = 2023
G = 9.80665


def sixhour_indices(days):
    idx = []
    for d in days:
        idx.extend([4 * d + i for i in range(4)])
    return idx


def sat_vapor_pressure_hpa(temp_c):
    return 6.112 * np.exp((17.67 * temp_c) / (temp_c + 243.5))


def wetbulb_stull_k(t2m_k, d2m_k):
    """Approximate 2m wet-bulb temperature in K from T and dewpoint."""
    t_c = t2m_k - 273.15
    td_c = d2m_k - 273.15
    rh = 100.0 * sat_vapor_pressure_hpa(td_c) / sat_vapor_pressure_hpa(t_c)
    rh = np.clip(rh, 1.0, 100.0)
    tw_c = (t_c * np.arctan(0.151977 * np.sqrt(rh + 8.313659))
            + np.arctan(t_c + rh)
            - np.arctan(rh - 1.676331)
            + 0.00391838 * rh ** 1.5 * np.arctan(0.023101 * rh)
            - 4.686035)
    return tw_c + 273.15


def phase_sl_mean(data_dir, year, days):
    ym = f"{year}04"
    path = os.path.join(data_dir, f"era5_sl_{ym}",
                        "data_stream-oper_stepType-instant_diag.nc")
    idx = sixhour_indices(days)
    with xr.open_dataset(path, engine="h5netcdf") as ds:
        t2m = ds["t2m"].isel(time=idx).mean("time").load()
        d2m = ds["d2m"].isel(time=idx).mean("time").load()
    return t2m, d2m


def phase_pl_mean(data_dir, year, days):
    ym = f"{year}04"
    idx = sixhour_indices(days)

    def read_var(short):
        path = os.path.join(data_dir, f"era5_pl_{short}_{ym}.nc")
        with xr.open_dataset(path, engine="h5netcdf") as ds:
            da = ds[short].isel(time=idx).mean("time").load()
        return da

    q = read_var("q")
    u = read_var("u")
    v = read_var("v")
    return q, u, v


def ivt_from_phase_means(q, u, v):
    pcoord = "pressure_level" if "pressure_level" in q.coords else "level"
    p = np.asarray(q[pcoord].values, dtype=np.float64) * 100.0
    order = np.argsort(p)
    p_sorted = p[order]
    qv = np.asarray(q.values, dtype=np.float64)[order]
    uv = np.asarray(u.values, dtype=np.float64)[order]
    vv = np.asarray(v.values, dtype=np.float64)[order]
    trapz = getattr(np, "trapezoid", np.trapz)
    ivtu = trapz(qv * uv, p_sorted, axis=0) / G
    ivtv = trapz(qv * vv, p_sorted, axis=0) / G
    return ivtu, ivtv


def area_mean(field, lat, lon, region):
    lat_min, lat_max = region["lat"]
    lon_min, lon_max = region["lon"]
    lat = np.asarray(lat)
    lon = np.asarray(lon)
    mask_lat = (lat >= lat_min) & (lat <= lat_max)
    mask_lon = (lon >= lon_min) & (lon <= lon_max)
    sub = field[np.ix_(mask_lat, mask_lon)]
    weights = np.cos(np.deg2rad(lat[mask_lat]))[:, None]
    valid = np.isfinite(sub)
    if not np.any(valid):
        return np.nan
    return np.nansum(sub * weights) / np.nansum(weights * valid)


def add_region_boxes(ax, cfg):
    boxes = [
        ("extend_region", "0.25", 1.0),
        ("contrast_region", "tab:orange", 1.4),
        ("key_region", "tab:blue", 1.6),
    ]
    for key, color, lw in boxes:
        reg = cfg["plot"].get(key)
        if not reg:
            continue
        lon0, lon1 = reg["lon"]
        lat0, lat1 = reg["lat"]
        ax.add_patch(Rectangle((lon0, lat0), lon1 - lon0, lat1 - lat0,
                               fill=False, edgecolor=color, linewidth=lw))


def plot_map(ax, lon, lat, field, title, levels, cmap, cfg):
    norm = mcolors.BoundaryNorm(levels, plt.get_cmap(cmap).N)
    im = ax.contourf(lon, lat, field, levels=levels, cmap=cmap, norm=norm,
                     extend="both")
    ax.set_title(title, fontsize=10, fontweight="bold")
    ax.set_xlim(cfg["plot"]["map_extent"]["lon"])
    ax.set_ylim(cfg["plot"]["map_extent"]["lat"])
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    add_region_boxes(ax, cfg)
    return im


def load_cfram_region_terms(case_name, cfg):
    from netCDF4 import Dataset

    path = os.path.join(PROJECT_ROOT, "cases", case_name, "output", "cfram_result.nc")
    if not os.path.exists(path):
        return None
    terms = {
        "WV": "dT_q",
        "Cloud": "dT_cloud",
        "Aerosol": "dT_aerosol",
        "Surface": "dT_sfcdyn",
        "AtmDyn": "dT_atmdyn",
        "Total": "dT_observed",
    }
    nc = Dataset(path)
    lat = np.asarray(nc.variables["lat"][:])
    lon = np.asarray(nc.variables["lon"][:])
    out = {"wet": {}, "dry": {}}
    for label, var in terms.items():
        if var not in nc.variables:
            continue
        arr = np.asarray(nc.variables[var][-1, :, :], dtype=np.float64)
        arr = np.where(np.abs(arr) > 900, np.nan, arr)
        out["wet"][label] = area_mean(arr, lat, lon, cfg["plot"]["key_region"])
        out["dry"][label] = area_mean(arr, lat, lon, cfg["plot"]["contrast_region"])
    nc.close()
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--case", default="india_wb23",
                    help="Base case with plot regions and output figures dir")
    args = ap.parse_args()

    cfg = load_case(args.case)
    data_dir = cfg["source"]["data_dir"]
    outdir = cfg["_figures_dir"]
    os.makedirs(outdir, exist_ok=True)

    phase_data = []
    for key, label, days, phase_case in PHASES:
        print(f"Phase {key}: {label}")
        t2m_ev, d2m_ev = phase_sl_mean(data_dir, EVENT_YEAR, days)
        tw_ev = wetbulb_stull_k(t2m_ev, d2m_ev)
        q_ev, u_ev, v_ev = phase_pl_mean(data_dir, EVENT_YEAR, days)
        ivtu_ev, ivtv_ev = ivt_from_phase_means(q_ev, u_ev, v_ev)
        q925_ev = q_ev.sel(pressure_level=925, method="nearest")

        t2m_clim_sum = None
        tw_clim_sum = None
        ivtu_clim_sum = None
        ivtv_clim_sum = None
        q925_clim_sum = None
        for year in YEARS:
            t2m_y, d2m_y = phase_sl_mean(data_dir, year, days)
            tw_y = wetbulb_stull_k(t2m_y, d2m_y)
            q_y, u_y, v_y = phase_pl_mean(data_dir, year, days)
            iu_y, iv_y = ivt_from_phase_means(q_y, u_y, v_y)
            q925_y = q_y.sel(pressure_level=925, method="nearest")
            t2m_clim_sum = t2m_y if t2m_clim_sum is None else t2m_clim_sum + t2m_y
            tw_clim_sum = tw_y if tw_clim_sum is None else tw_clim_sum + tw_y
            ivtu_clim_sum = iu_y if ivtu_clim_sum is None else ivtu_clim_sum + iu_y
            ivtv_clim_sum = iv_y if ivtv_clim_sum is None else ivtv_clim_sum + iv_y
            q925_clim_sum = q925_y if q925_clim_sum is None else q925_clim_sum + q925_y

        n = float(len(YEARS))
        t2m_clim = t2m_clim_sum / n
        tw_clim = tw_clim_sum / n
        ivtu_clim = ivtu_clim_sum / n
        ivtv_clim = ivtv_clim_sum / n
        q925_clim = q925_clim_sum / n

        lon = np.asarray(t2m_ev.longitude.values)
        lat = np.asarray(t2m_ev.latitude.values)
        ivtmag_ev = np.sqrt(ivtu_ev ** 2 + ivtv_ev ** 2)
        ivtmag_clim = np.sqrt(ivtu_clim ** 2 + ivtv_clim ** 2)
        data = {
            "key": key,
            "label": label,
            "case": phase_case,
            "lat": lat,
            "lon": lon,
            "t2m_anom": np.asarray(t2m_ev - t2m_clim),
            "tw_anom": np.asarray(tw_ev - tw_clim),
            "q925_anom": np.asarray((q925_ev - q925_clim) * 1000.0),
            "ivtu_anom": np.asarray(ivtu_ev - ivtu_clim),
            "ivtv_anom": np.asarray(ivtv_ev - ivtv_clim),
            "ivtmag_anom": np.asarray(ivtmag_ev - ivtmag_clim),
        }
        data["region"] = {
            "wet_tw": area_mean(data["tw_anom"], lat, lon, cfg["plot"]["key_region"]),
            "dry_tw": area_mean(data["tw_anom"], lat, lon, cfg["plot"]["contrast_region"]),
            "wet_t2m": area_mean(data["t2m_anom"], lat, lon, cfg["plot"]["key_region"]),
            "dry_t2m": area_mean(data["t2m_anom"], lat, lon, cfg["plot"]["contrast_region"]),
            "wet_q925": area_mean(data["q925_anom"], lat, lon, cfg["plot"]["key_region"]),
            "dry_q925": area_mean(data["q925_anom"], lat, lon, cfg["plot"]["contrast_region"]),
            "wet_ivt": area_mean(data["ivtmag_anom"], lat, lon, cfg["plot"]["key_region"]),
            "dry_ivt": area_mean(data["ivtmag_anom"], lat, lon, cfg["plot"]["contrast_region"]),
        }
        data["cfram"] = load_cfram_region_terms(phase_case, cfg)
        phase_data.append(data)

    # Wet-heat maps.
    fig, axes = plt.subplots(2, 3, figsize=(15, 7), constrained_layout=True)
    levels_t = np.arange(-6, 6.5, 0.5)
    for col, data in enumerate(phase_data):
        im = plot_map(axes[0, col], data["lon"], data["lat"], data["t2m_anom"],
                      f"T2m anomaly {data['label']}", levels_t, "RdBu_r", cfg)
        im = plot_map(axes[1, col], data["lon"], data["lat"], data["tw_anom"],
                      f"Tw anomaly {data['label']}", levels_t, "RdBu_r", cfg)
    fig.colorbar(im, ax=axes.ravel().tolist(), orientation="horizontal",
                 shrink=0.75, pad=0.04, label="K")
    out = os.path.join(outdir, "wet_heat_anomaly_phases.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print("Saved", out)

    # IVT maps.
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2), constrained_layout=True)
    levels_ivt = np.arange(-160, 180, 20)
    for col, data in enumerate(phase_data):
        ax = axes[col]
        im = plot_map(ax, data["lon"], data["lat"], data["ivtmag_anom"],
                      f"IVT anomaly {data['label']}", levels_ivt, "BrBG", cfg)
        step = 8
        ax.quiver(data["lon"][::step], data["lat"][::step],
                  data["ivtu_anom"][::step, ::step],
                  data["ivtv_anom"][::step, ::step],
                  scale=900, width=0.002, color="0.15")
    fig.colorbar(im, ax=axes.ravel().tolist(), orientation="horizontal",
                 shrink=0.75, pad=0.04, label="kg m-1 s-1")
    out = os.path.join(outdir, "moisture_transport_phases.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print("Saved", out)

    # Regional contrast.
    fig, axes = plt.subplots(3, 2, figsize=(13, 12), constrained_layout=True)
    axes = axes.ravel()
    x = np.arange(len(PHASES))
    width = 0.36
    metrics = [
        ("Tw", "tw", "K"),
        ("T2m", "t2m", "K"),
        ("q925", "q925", "g kg-1"),
        ("IVT", "ivt", "kg m-1 s-1"),
    ]
    for ax, (name, key, unit) in zip(axes[:4], metrics):
        wet = [d["region"][f"wet_{key}"] for d in phase_data]
        dry = [d["region"][f"dry_{key}"] for d in phase_data]
        ax.bar(x - width / 2, wet, width, color="tab:blue", label="Wet-heat region")
        ax.bar(x + width / 2, dry, width, color="tab:orange", label="Dry-heat region")
        ax.axhline(0, color="0.4", linewidth=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels([d["label"] for d in phase_data])
        ax.set_ylabel(unit)
        ax.set_title(f"{name} anomaly")
        ax.grid(axis="y", alpha=0.3)
        ax.legend(fontsize=8)

    terms = ["WV", "Cloud", "Aerosol", "Surface", "AtmDyn", "Total"]
    colors = plt.get_cmap("tab10").colors
    term_x = np.arange(len(terms))
    phase_width = 0.24
    for ax, region_key, title in [
            (axes[4], "wet", "CFRAM: wet-heat region"),
            (axes[5], "dry", "CFRAM: dry-heat region")]:
        for i, data in enumerate(phase_data):
            values = []
            for term in terms:
                if data["cfram"] is None:
                    values.append(np.nan)
                else:
                    values.append(data["cfram"][region_key].get(term, np.nan))
            ax.bar(term_x + (i - 1) * phase_width, values, phase_width,
                   color=colors[i], label=data["label"])
        ax.axhline(0, color="0.4", linewidth=0.8)
        ax.set_xticks(term_x)
        ax.set_xticklabels(terms, rotation=25)
        ax.set_ylabel("Surface contribution (K)")
        ax.set_title(title)
        ax.grid(axis="y", alpha=0.3)
        ax.legend(fontsize=8)
    fig.suptitle("Wet-heat vs dry-heat regional contrast", fontsize=15)
    out = os.path.join(outdir, "regional_contrast_phases.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print("Saved", out)


if __name__ == "__main__":
    main()
