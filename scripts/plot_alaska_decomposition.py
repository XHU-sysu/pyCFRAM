#!/usr/bin/env python3
"""Plot Alaska wildfire CFRAM decomposition maps."""

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.config import load_case
from scripts.plot_fig3_self import (
    CMAP,
    LEVELS_DT,
    PLOT_ROWS,
    get_map_extent,
    load_case_data,
)


REGION_KEYS = [
    ("extended_region", "0.25"),
    ("arctic_response_region", "tab:orange"),
    ("main_fire_region", "tab:blue"),
]


def region_box(plot_cfg, name):
    region = plot_cfg.get(name)
    if not region:
        return None
    lon = [x - 360.0 if x > 180.0 else x for x in region["lon"]]
    return lon + region["lat"]


def shift_longitudes(lons, data):
    shifted = np.where(np.asarray(lons) > 180.0, np.asarray(lons) - 360.0, lons)
    order = np.argsort(shifted)
    return shifted[order], data[:, order]


def shift_extent(extent):
    lon = [x - 360.0 if x > 180.0 else x for x in extent[:2]]
    return lon + extent[2:]


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


def plot_alaska_panel(ax, lons, lats, field, title, norm, cmap, extent, boxes):
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    cf = ax.contourf(
        lons,
        lats,
        field,
        levels=LEVELS_DT,
        cmap=cmap,
        norm=norm,
        transform=ccrs.PlateCarree(),
        extend="both",
    )
    ax.coastlines(resolution="50m", linewidth=0.6)
    ax.add_feature(cfeature.BORDERS, linewidth=0.35, edgecolor="0.35")
    ax.add_feature(cfeature.LAKES, linewidth=0.25, edgecolor="0.5",
                   facecolor="none")
    add_box(ax, boxes["extended_region"], "0.25", 1.0)
    add_box(ax, boxes["arctic_response_region"], "tab:orange", 1.2)
    add_box(ax, boxes["main_fire_region"], "tab:blue", 1.6)
    gl = ax.gridlines(
        draw_labels=False,
        linewidth=0.3,
        color="gray",
        alpha=0.5,
        linestyle="--",
    )
    ax.set_title(title, fontsize=9, fontweight="bold")
    return cf


def main():
    case_name = sys.argv[1] if len(sys.argv) > 1 else "alaska_wf22"
    cfg = load_case(case_name)
    lats, lons, data = load_case_data(case_name)
    plot_cfg = cfg.get("plot", {})
    extent = shift_extent(get_map_extent(cfg))

    boxes = {key: region_box(plot_cfg, key) for key, _ in REGION_KEYS}
    norm = mcolors.BoundaryNorm(LEVELS_DT, CMAP.N, clip=True)
    labels = "abcdefghijklmnopqrstuvwxyz"

    projection = ccrs.LambertConformal(
        central_longitude=-150,
        central_latitude=66,
        standard_parallels=(55, 75),
    )
    fig, axes = plt.subplots(
        len(PLOT_ROWS),
        1,
        figsize=(7.0, 14.0),
        subplot_kw={"projection": projection},
        squeeze=False,
    )

    for row, (term, term_label) in enumerate(PLOT_ROWS):
        ax = axes[row, 0]
        field = data.get(term, np.full((len(lats), len(lons)), np.nan))
        field = np.clip(np.nan_to_num(field, nan=0.0), -20, 20)
        plot_lons, field = shift_longitudes(lons, field)
        title = f"({labels[row]}) {term_label}"
        plot_alaska_panel(
            ax,
            plot_lons,
            lats,
            field,
            title,
            norm,
            CMAP,
            extent,
            boxes,
        )

    fig.suptitle(
        "Alaska Wildfire 2022: CFRAM Surface Temperature Decomposition",
        fontsize=14,
        fontweight="bold",
        y=0.995,
    )
    fig.subplots_adjust(
        bottom=0.055,
        top=0.965,
        left=0.08,
        right=0.95,
        hspace=0.25,
    )
    cbar_ax = fig.add_axes([0.15, 0.02, 0.7, 0.014])
    colorbar = fig.colorbar(
        plt.cm.ScalarMappable(norm=norm, cmap=CMAP),
        cax=cbar_ax,
        orientation="horizontal",
        ticks=LEVELS_DT,
    )
    colorbar.set_label("Partial Temperature Contribution (K)", fontsize=10)
    colorbar.ax.tick_params(labelsize=8)

    outdir = cfg["_figures_dir"]
    os.makedirs(outdir, exist_ok=True)
    outpath = os.path.join(outdir, "fig3_decomposition.png")
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {outpath}")


if __name__ == "__main__":
    main()
