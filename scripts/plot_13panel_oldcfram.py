#!/usr/bin/env python3
"""Plot 13-panel surface dT decomposition from OLD Fortran CFRAM raw output
(xiaominghu CESM2 4xCO2 collab data). Generates BOTH global and north-polar
figures in one run.

Mirrors plot_13panel_global.py / plot_13panel_polar.py layout/colorbar, but
reads partial_T_1.grd directly (13 terms × 18 levs × 192 × 288, float32).

Usage:
    python3 scripts/plot_13panel_oldcfram.py
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.path as mpath
import cartopy.crs as ccrs

OLD_DIR = '/disk/r074/lzhenn/workspace/ust-jumper/pyCFRAM/raw_data/cesm2_cmip6/CFRAM_xiaominghu'
OUT_DIR = '/disk/r074/lzhenn/workspace/ust-jumper/pyCFRAM/cases/cesm2_4xco2_official_17p_old/figures'

IX, IY = 288, 192
ZD18, NLEV, NTERM = 18, 17, 13
SFC_IDX = 17                                # last level = surface row (TOA→sfc)

# OLD CFRAM record order in partial_T_1.grd
TERM_ORDER = ['albedo','wv','cloud','cloud_sw','cloud_lw','co2','o3','solar',
              'dyn','atm_dyn','sfc_dyn','shflx','lhflx']

# (panel_label, OLD_term_name, panel_letter) — mirror pyCFRAM panel order
PANELS = [
    ('AL',   'albedo',   'a'),
    ('WV',   'wv',       'b'),
    ('CLD',  'cloud',    'c'),
    ('CLDS', 'cloud_sw', 'd'),
    ('CLDL', 'cloud_lw', 'e'),
    ('CO2',  'co2',      'f'),
    ('O3',   'o3',       'g'),
    ('SR',   'solar',    'h'),
    ('DYN',  'dyn',      'i'),    # pyCFRAM "DYN" panel maps to OLD dt_dyn
    ('ATM',  'atm_dyn',  'j'),
    ('OCH',  'sfc_dyn',  'k'),    # pyCFRAM "OCH" maps to OLD dt_sfc_dyn
    ('SH',   'shflx',    'l'),
    ('LH',   'lhflx',    'm'),
]

# Color scheme: identical to pyCFRAM 13-panel scripts
LEVELS = [-15.0, -7.0, -2.0, -0.5, 0.0, 0.5, 2.0, 7.0, 15.0]
_BIN_COLORS = [
    '#053061', '#2166ac', '#92c5de', '#d1e5f0',
    '#fddbc7', '#f4a582', '#b2182b', '#67001f',
]
CMAP = mcolors.ListedColormap(_BIN_COLORS)
CMAP.set_under('#021A40'); CMAP.set_over('#400015')

# CESM2 f09 grid (matches xiaominghu data)
LAT = np.linspace(-90 + 90.0/IY, 90 - 90.0/IY, IY)
LON = np.linspace(0, 360 - 360.0/IX, IX)


def load_old_surface():
    """Read partial_T_1.grd → dict of {term: 2D sfc field}."""
    path = os.path.join(OLD_DIR, 'partial_T_1.grd')
    raw = np.fromfile(path, dtype=np.float32).reshape((NTERM, ZD18, IY, IX))
    raw = np.where(np.abs(raw) > 900, np.nan, raw.astype(np.float64))
    out = {}
    for i, term in enumerate(TERM_ORDER):
        out[term] = raw[i, SFC_IDX, :, :]    # surface row
    return out


def plot_global(data, out_png):
    fig = plt.figure(figsize=(20, 14))
    proj = ccrs.PlateCarree(central_longitude=0)
    for idx, (label, term, letter) in enumerate(PANELS):
        row, col = divmod(idx, 4)
        ax = fig.add_subplot(4, 4, row*4 + col + 1, projection=proj)
        field = data[term]
        cf = ax.contourf(LON, LAT, field, levels=LEVELS, colors=_BIN_COLORS,
                         extend='both', transform=ccrs.PlateCarree())
        cf.cmap.set_under('#021A40'); cf.cmap.set_over('#400015')
        ax.coastlines(resolution='110m', linewidth=0.4, color='k')
        ax.set_global()
        ax.set_title(f'({letter}) {label}', fontsize=11, loc='left')
        gl = ax.gridlines(draw_labels=True, linewidth=0.3, alpha=0.4)
        gl.top_labels = False; gl.right_labels = False
        gl.xlabel_style = {'size': 7}; gl.ylabel_style = {'size': 7}
        dm = np.nanmean(field)
        ax.text(0.02, 0.02, f'mean={dm:+.2f}K',
                transform=ax.transAxes, fontsize=8,
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
    cbar_ax = fig.add_axes([0.15, 0.04, 0.7, 0.015])
    cb = fig.colorbar(cf, cax=cbar_ax, orientation='horizontal',
                      ticks=LEVELS, extend='both')
    cb.set_label('Partial Temperature dT (K)', fontsize=11)
    fig.suptitle('OLD Fortran CFRAM 13-panel surface dT decomposition: '
                 'CESM2 4xCO2 (xiaominghu collab, with O3=huss corruption)',
                 fontsize=13, y=0.97)
    plt.tight_layout(rect=[0, 0.07, 1, 0.95])
    plt.savefig(out_png, dpi=150, bbox_inches='tight')
    print('Saved:', out_png)
    plt.close(fig)


def _make_circle_boundary():
    """Circular boundary for polar stereographic axes."""
    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    return mpath.Path(verts*radius + center)


def plot_polar(data, out_png):
    """North polar stereographic — same style as pyCFRAM plot_13panel_polar.py:
    contourf + circle boundary + figsize (20, 22)."""
    LAT_LIMIT = 60
    proj = ccrs.NorthPolarStereo(central_longitude=0)
    fig = plt.figure(figsize=(20, 22))
    circle = _make_circle_boundary()
    # Append cyclic lon column
    lon_w = np.append(LON, LON[0] + 360.0)
    cf = None
    for idx, (label, term, letter) in enumerate(PANELS):
        ax = fig.add_subplot(4, 4, idx+1, projection=proj)
        ax.set_extent([-180, 180, LAT_LIMIT, 90], crs=ccrs.PlateCarree())
        ax.set_boundary(circle, transform=ax.transAxes)
        field = data[term]
        # Append cyclic column to data
        field_w = np.concatenate([field, field[:, :1]], axis=1)
        if not np.all(np.isnan(field_w)):
            # Use pcolormesh (NOT contourf): shapely 2.x + cartopy 0.22
            # MultiPolygon bug can silently drop contour bins on polar projection.
            norm = mcolors.BoundaryNorm(LEVELS, len(_BIN_COLORS), clip=False)
            cf = ax.pcolormesh(lon_w, LAT, field_w,
                               norm=norm, cmap=CMAP, shading='auto',
                               transform=ccrs.PlateCarree())
        ax.coastlines(resolution='110m', linewidth=0.4, color='k')
        ax.gridlines(draw_labels=False, linewidth=0.3, alpha=0.4)
        ax.set_title(f'({letter}) {label}', fontsize=12, loc='left')
    cbar_ax = fig.add_axes([0.15, 0.06, 0.7, 0.012])
    fig.colorbar(cf, cax=cbar_ax, orientation='horizontal',
                 ticks=LEVELS, extend='both',
                 label='Partial Temperature dT (K)')
    fig.suptitle('OLD Fortran CFRAM 13-panel surface dT decomposition: '
                 'CESM2 4xCO2 (xiaominghu collab, O3=huss corruption)\n'
                 f'North polar stereographic, |lat| ≥ {LAT_LIMIT}°',
                 fontsize=14, y=0.97)
    plt.tight_layout(rect=[0, 0.08, 1, 0.94])
    plt.savefig(out_png, dpi=130, bbox_inches='tight')
    print('Saved:', out_png)
    plt.close(fig)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    data = load_old_surface()
    plot_global(data, os.path.join(OUT_DIR, 'fig_13panel_global_oldcfram.png'))
    plot_polar(data,  os.path.join(OUT_DIR, 'fig_13panel_north_polar_oldcfram.png'))


if __name__ == '__main__':
    main()
