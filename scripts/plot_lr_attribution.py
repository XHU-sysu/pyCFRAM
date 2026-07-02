#!/usr/bin/env python3
"""Plot M3 per-process Lapse-Rate attribution (docs/plan.md WP-M3.3).

Reads cases/<case>/output/lr_attribution.nc (from
scripts/compute_lr_attribution.py) and produces:
  fig_lr_attribution_<kernel>.png   -- per-process ΔR_LR maps, one kernel
  fig_lr_zonal_profile.png          -- zonal-mean ΔR_LR by process (both
                                        kernels), analogous to the paper's
                                        Fig 3a zonal-profile style

Usage:
    python scripts/plot_lr_attribution.py cesm2_4xco2_official
"""
import os
import sys
import numpy as np
from netCDF4 import Dataset
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.config import load_case

FILL = -999.0
LEVELS = np.arange(-5, 5.1, 0.5)
CMAP = plt.get_cmap('RdBu_r')

TERM_COLORS = {
    'q':       '#1f77b4',
    'co2':     '#d62728',
    'o3':      '#9467bd',
    'solar':   '#8c564b',
    'albedo':  '#e377c2',
    'cloud':   '#7f7f7f',
    'aerosol': '#bcbd22',
    'lhflx':   '#17becf',
    'shflx':   '#1abc9c',
    'atmdyn':  '#2ca02c',
    'ocndyn':  '#34495e',
}


def _mask_fill(arr):
    return np.where(np.abs(arr) > 900.0, np.nan, arr)


def area_weighted_zonal_mean(field, lat, lon):
    """(lat,lon) -> (lat,) zonal mean (simple lon average; area weighting
    is already implicit since all lon bins share the same cos(lat))."""
    return np.nanmean(field, axis=1)


def main():
    if len(sys.argv) < 2:
        print("Usage: plot_lr_attribution.py <case>")
        sys.exit(1)
    case = sys.argv[1]
    cfg = load_case(case)

    attr_nc = os.path.join(cfg['_output_dir'], 'lr_attribution.nc')
    if not os.path.exists(attr_nc):
        sys.exit("Missing %s -- run scripts/compute_lr_attribution.py %s" % (attr_nc, case))
    os.makedirs(cfg['_figures_dir'], exist_ok=True)

    with Dataset(attr_nc) as d:
        lat = np.array(d.variables['lat'][:])
        lon = np.array(d.variables['lon'][:])
        terms_attr = d.terms.split(',') if hasattr(d, 'terms') else []
        by_kernel = {}
        for vname in d.variables:
            if not vname.startswith('dR_lr_from_'):
                continue
            rest = vname[len('dR_lr_from_'):]
            # rest = "<term>_<kernel>"; kernel names are known, terms are
            # in terms_attr -- split on the longest matching term prefix.
            for term in sorted(terms_attr, key=len, reverse=True):
                if rest.startswith(term + '_'):
                    kname = rest[len(term) + 1:]
                    by_kernel.setdefault(kname, {})[term] = _mask_fill(
                        np.array(d.variables[vname][:]))
                    break

    # ---- per-process maps, one figure per kernel ----
    for kname, term_fields in by_kernel.items():
        terms = sorted(term_fields.keys())
        ncols = 3
        nrows = int(np.ceil(len(terms) / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 2.6 * nrows),
                                  subplot_kw={'projection': ccrs.PlateCarree()},
                                  squeeze=False)
        norm = mcolors.BoundaryNorm(LEVELS, CMAP.N, clip=True)
        for i, term in enumerate(terms):
            ax = axes[i // ncols, i % ncols]
            field = term_fields[term]
            ax.set_global()
            cf = ax.contourf(lon, lat, np.nan_to_num(field, nan=0.0), levels=LEVELS,
                              cmap=CMAP, norm=norm, transform=ccrs.PlateCarree(), extend='both')
            ax.coastlines(resolution='110m', linewidth=0.4)
            ax.add_feature(cfeature.BORDERS, linewidth=0.2, edgecolor='gray')
            ax.set_title(term, fontsize=9, fontweight='bold')
        for j in range(len(terms), nrows * ncols):
            axes[j // ncols, j % ncols].axis('off')
        fig.colorbar(cf, ax=axes, orientation='horizontal', fraction=0.03, pad=0.05,
                     label='ΔR_LR from process, W m$^{-2}$')
        fig.suptitle("%s: per-process lapse-rate attribution (%s kernel)" % (case, kname),
                     fontsize=12)
        out = os.path.join(cfg['_figures_dir'], 'fig_lr_attribution_%s.png' % kname)
        fig.savefig(out, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print("Wrote %s" % out)

    # ---- zonal-mean profile (paper Fig 3a style) ----
    fig, axes = plt.subplots(1, len(by_kernel), figsize=(6 * len(by_kernel), 5), squeeze=False)
    for k_i, (kname, term_fields) in enumerate(sorted(by_kernel.items())):
        ax = axes[0, k_i]
        for term in sorted(term_fields.keys()):
            zm = area_weighted_zonal_mean(term_fields[term], lat, lon)
            ax.plot(zm, lat, label=term, color=TERM_COLORS.get(term), linewidth=1.3)
        ax.axvline(0, color='k', linewidth=0.5)
        ax.set_xlabel('ΔR_LR (W m$^{-2}$)')
        ax.set_ylabel('Latitude')
        ax.set_title('%s kernel' % kname)
        ax.legend(fontsize=7, ncol=2, loc='best')
        ax.grid(alpha=0.3)
    fig.suptitle("%s: zonal-mean per-process lapse-rate attribution" % case)
    out = os.path.join(cfg['_figures_dir'], 'fig_lr_zonal_profile.png')
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print("Wrote %s" % out)


if __name__ == '__main__':
    main()
