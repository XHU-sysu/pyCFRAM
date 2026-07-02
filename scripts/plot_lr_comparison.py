#!/usr/bin/env python3
"""Plot native vs ClimKern Lapse-Rate ΔR_LR comparison (docs/plan.md WP-M2.5).

Reads cases/<case>/output/lr_kernel.nc (native) and
cases/<case>/output/lr_climkern_ref.nc (ClimKern reference, produced by
scripts/validate_lr_vs_climkern.py) -- no xesmf/climkern import needed
here, just the two NetCDF outputs.

Produces:
  fig_lr_comparison_<kernel>.png   -- native / ClimKern / diff triptych, per kernel
  fig_lr_kramer_vs_gfdl.png        -- native CloudSat(Kramer) vs GFDL difference

Usage:
    python scripts/plot_lr_comparison.py cesm2_4xco2_official
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
LEVELS = np.arange(-20, 20.1, 2)
CMAP = plt.get_cmap('RdBu_r')
DIFF_LEVELS = np.arange(-5, 5.1, 0.5)


def _mask_fill(arr):
    return np.where(np.abs(arr) > 900.0, np.nan, arr)


def area_weighted_corr(a, b, lat):
    w = np.cos(np.deg2rad(lat))
    w2d = np.broadcast_to(w[:, None], a.shape)
    mask = ~np.isnan(a) & ~np.isnan(b)
    aw, bw, ww = a[mask], b[mask], w2d[mask]
    ww = ww / ww.sum()
    ma, mb = np.sum(aw * ww), np.sum(bw * ww)
    cov = np.sum(ww * (aw - ma) * (bw - mb))
    va, vb = np.sum(ww * (aw - ma) ** 2), np.sum(ww * (bw - mb) ** 2)
    return cov / np.sqrt(va * vb)


def plot_panel(ax, lon, lat, field, title, levels, cmap):
    norm = mcolors.BoundaryNorm(levels, cmap.N, clip=True)
    ax.set_global()
    cf = ax.contourf(lon, lat, np.nan_to_num(field, nan=0.0), levels=levels,
                      cmap=cmap, norm=norm, transform=ccrs.PlateCarree(), extend='both')
    ax.coastlines(resolution='110m', linewidth=0.5)
    ax.add_feature(cfeature.BORDERS, linewidth=0.3, edgecolor='gray')
    ax.set_title(title, fontsize=9, fontweight='bold')
    return cf


def main():
    if len(sys.argv) < 2:
        print("Usage: plot_lr_comparison.py <case>")
        sys.exit(1)
    case = sys.argv[1]
    cfg = load_case(case)

    native_nc = os.path.join(cfg['_output_dir'], 'lr_kernel.nc')
    ref_nc = os.path.join(cfg['_output_dir'], 'lr_climkern_ref.nc')
    if not os.path.exists(native_nc):
        sys.exit("Missing %s -- run scripts/compute_lr_kernel.py %s" % (native_nc, case))
    if not os.path.exists(ref_nc):
        sys.exit("Missing %s -- run scripts/validate_lr_vs_climkern.py %s" % (ref_nc, case))

    os.makedirs(cfg['_figures_dir'], exist_ok=True)

    with Dataset(native_nc) as d:
        lat = np.array(d.variables['lat'][:])
        lon = np.array(d.variables['lon'][:])
        native_vars = {v[6:]: _mask_fill(np.array(d.variables[v][:]))
                        for v in d.variables if v.startswith('dR_lr_')}
    with Dataset(ref_nc) as d:
        ref_vars = {v[6:]: _mask_fill(np.array(d.variables[v][:]))
                    for v in d.variables if v.startswith('dR_lr_')}

    kernel_names = sorted(set(native_vars) & set(ref_vars))
    if not kernel_names:
        sys.exit("No matching kernels between %s and %s" % (native_nc, ref_nc))

    for kname in kernel_names:
        native = native_vars[kname]
        ref = ref_vars[kname]
        diff = native - ref
        corr = area_weighted_corr(native, ref, lat)

        fig, axes = plt.subplots(1, 3, figsize=(16, 4),
                                  subplot_kw={'projection': ccrs.PlateCarree()})
        cf0 = plot_panel(axes[0], lon, lat, native, "Native ΔR_LR (%s)" % kname, LEVELS, CMAP)
        cf1 = plot_panel(axes[1], lon, lat, ref, "ClimKern ΔR_LR (%s)" % kname, LEVELS, CMAP)
        cf2 = plot_panel(axes[2], lon, lat, diff,
                          "Native − ClimKern (corr=%.3f)" % corr, DIFF_LEVELS, CMAP)
        fig.colorbar(cf0, ax=axes[0:2], orientation='horizontal', fraction=0.05, pad=0.08,
                     label='W m$^{-2}$')
        fig.colorbar(cf2, ax=axes[2], orientation='horizontal', fraction=0.05, pad=0.08,
                     label='W m$^{-2}$')
        fig.suptitle("%s: %s lapse-rate kernel comparison" % (case, kname), fontsize=11)
        out = os.path.join(cfg['_figures_dir'], 'fig_lr_comparison_%s.png' % kname)
        fig.savefig(out, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print("Wrote %s (corr=%.4f)" % (out, corr))

    if 'CloudSat' in native_vars and 'GFDL' in native_vars:
        diff_kernels = native_vars['CloudSat'] - native_vars['GFDL']
        fig, ax = plt.subplots(1, 1, figsize=(8, 4.5),
                                subplot_kw={'projection': ccrs.PlateCarree()})
        cf = plot_panel(ax, lon, lat, diff_kernels,
                         "Native ΔR_LR: CloudSat(Kramer) − GFDL", DIFF_LEVELS, CMAP)
        fig.colorbar(cf, ax=ax, orientation='horizontal', fraction=0.05, pad=0.08,
                     label='W m$^{-2}$')
        out = os.path.join(cfg['_figures_dir'], 'fig_lr_kramer_vs_gfdl.png')
        fig.savefig(out, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print("Wrote %s" % out)


if __name__ == '__main__':
    main()
