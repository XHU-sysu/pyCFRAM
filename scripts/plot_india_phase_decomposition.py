#!/usr/bin/env python3
"""Plot India-Bangladesh CFRAM decomposition for the three April 2023 phases."""

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.config import load_case
from scripts.plot_fig3_self import (
    CMAP,
    LEVELS_DT,
    PLOT_ROWS,
    get_map_extent,
    load_case_data,
    plot_panel,
)


PHASES = [
    ("india_wb23_pre", "Pre-event: Apr 1-16"),
    ("india_wb23_core", "Core: Apr 17-20"),
    ("india_wb23_post", "Post-event: Apr 21-30"),
]


def region_box(plot_cfg, name):
    region = plot_cfg.get(name)
    return region["lon"] + region["lat"] if region else None


def main():
    norm = mcolors.BoundaryNorm(LEVELS_DT, CMAP.N, clip=True)
    nrows = len(PLOT_ROWS)
    ncols = len(PHASES)
    labels = "abcdefghijklmnopqrstuvwxyz"

    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(15, 18),
        subplot_kw={"projection": ccrs.PlateCarree()},
        squeeze=False,
    )

    for col, (case_name, phase_label) in enumerate(PHASES):
        cfg = load_case(case_name)
        lats, lons, data = load_case_data(case_name)
        plot_cfg = cfg.get("plot", {})
        extent = get_map_extent(cfg)
        key_region = region_box(plot_cfg, "key_region")
        extend_region = region_box(plot_cfg, "extend_region")
        contrast_region = region_box(plot_cfg, "contrast_region")

        for row, (term, term_label) in enumerate(PLOT_ROWS):
            field = data.get(term, np.full((len(lats), len(lons)), np.nan))
            field = np.clip(np.nan_to_num(field, nan=0.0), -20, 20)
            panel = labels[row * ncols + col]
            title = f"({panel}) {term_label} | {phase_label}"
            plot_panel(
                axes[row, col],
                lons,
                lats,
                field,
                title,
                norm,
                CMAP,
                extent,
                key_region,
                extend_region,
                contrast_region,
            )

    fig.suptitle(
        "India-Bangladesh Wet-Heat Event: Three-Phase CFRAM Decomposition",
        fontsize=15,
        fontweight="bold",
        y=0.995,
    )
    fig.subplots_adjust(
        bottom=0.055,
        top=0.965,
        left=0.045,
        right=0.97,
        hspace=0.25,
        wspace=0.12,
    )
    cbar_ax = fig.add_axes([0.17, 0.02, 0.66, 0.014])
    colorbar = fig.colorbar(
        plt.cm.ScalarMappable(norm=norm, cmap=CMAP),
        cax=cbar_ax,
        orientation="horizontal",
        ticks=LEVELS_DT,
    )
    colorbar.set_label("Partial Temperature Contribution (K)", fontsize=10)
    colorbar.ax.tick_params(labelsize=8)

    outdir = load_case("india_wb23")["_figures_dir"]
    os.makedirs(outdir, exist_ok=True)
    outpath = os.path.join(outdir, "fig3_decomposition_phases.png")
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {outpath}")


if __name__ == "__main__":
    main()
