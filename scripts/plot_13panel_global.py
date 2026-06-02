#!/usr/bin/env python3
"""Plot 13-panel global surface dT decomposition matching OLD Fortran CFRAM layout.

Layout:
  (a) AL   (b) WV   (c) CLD   (d) CLDS
  (e) CLDL (f) CO2  (g) O3    (h) SR
  (i) DYN  (j) ATM  (k) OCH   (l) SH
  (m) LH

Usage:
    python3 scripts/plot_13panel_global.py <case_name>
"""
import os, sys
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

# Mirror the OLD CFRAM colorbar (asymmetric, with overflow arrows). Each bin
# uses one tier deeper than stock RdBu_r — (7, 15] gets the deepest red that
# would otherwise be extend-only, etc. Center two bins use very-pale tints
# (sign-preserving near-white) so panels with frc≈0 render essentially blank.
LEVELS = [-15.0, -7.0, -2.0, -0.5, 0.0, 0.5, 2.0, 7.0, 15.0]
_BIN_COLORS = [
    '#053061',   # (-15, -7]    deepest blue  (was extend < -15)
    '#2166ac',   # (-7, -2]     2nd deepest   (was (-15, -7])
    '#92c5de',   # (-2, -0.5]   maintain
    '#d1e5f0',   # (-0.5, 0]    very pale blue
    '#fddbc7',   # (0, 0.5]     very pale red
    '#f4a582',   # (0.5, 2]     maintain
    '#b2182b',   # (2, 7]       2nd deepest   (was (7, 15])
    '#67001f',   # (7, 15]      deepest red   (was extend > 15)
]
CMAP = mcolors.ListedColormap(_BIN_COLORS)
CMAP.set_under('#021A40')   # < -15: even darker navy
CMAP.set_over('#400015')    # > 15:  even darker red

# (label, var_or_callable, panel_letter)
PANELS = [
    ('AL',   lambda d: d['albedo'],                                'a'),
    ('WV',   lambda d: d['q'],                                     'b'),
    ('CLD',  lambda d: d['cloud'],                                 'c'),
    ('CLDS', lambda d: d['cloud_sw'],                              'd'),
    ('CLDL', lambda d: d['cloud_lw'],                              'e'),
    ('CO2',  lambda d: d['co2'],                                   'f'),
    ('O3',   lambda d: d['o3'],                                    'g'),
    ('SR',   lambda d: d['solar'],                                 'h'),
    ('DYN',  lambda d: d['dyn'],                                   'i'),
    ('ATM',  lambda d: d['atmdyn'],                                'j'),
    ('OCH',  lambda d: d['ocndyn'],                                'k'),
    ('SH',   lambda d: d['shflx'],                                 'l'),
    ('LH',   lambda d: d['lhflx'],                                 'm'),
]


def nansum(arrs):
    out = np.zeros_like(arrs[0])
    cnt = np.zeros_like(arrs[0])
    for a in arrs:
        m = ~np.isnan(a)
        out[m] += a[m]
        cnt[m] += 1
    out[cnt == 0] = np.nan
    return out


def load_surface(nc, varname):
    if varname not in nc.variables:
        return None
    arr = np.array(nc.variables[varname][-1, :, :], dtype=np.float64)
    arr = np.where(np.abs(arr) > 900, np.nan, arr)
    return arr


def main(case_name):
    cfg = load_case(case_name)
    out_file = os.path.join(cfg['_output_dir'], 'cfram_result.nc')
    if not os.path.exists(out_file):
        sys.exit('output not found: ' + out_file)

    nc = Dataset(out_file)
    lats = np.array(nc.variables['lat'][:])
    lons = np.array(nc.variables['lon'][:])

    var_map = {
        'albedo':   'dT_albedo', 'q':       'dT_q',
        'cloud':    'dT_cloud',  'cloud_sw':'dT_cloud_sw',
        'cloud_lw': 'dT_cloud_lw','co2':    'dT_co2',
        'o3':       'dT_o3',     'solar':   'dT_solar',
        'atmdyn':   'dT_atmdyn', 'sfcdyn':  'dT_sfcdyn',
        'ocndyn':   'dT_ocndyn', 'shflx':   'dT_shflx',
        'lhflx':    'dT_lhflx',  'dyn':     'dT_dry',
    }
    data = {}
    for k, v in var_map.items():
        a = load_surface(nc, v)
        data[k] = a if a is not None else np.full((len(lats), len(lons)), np.nan)
    nc.close()

    # Plot
    fig = plt.figure(figsize=(20, 14))
    proj = ccrs.PlateCarree(central_longitude=0)

    for idx, (label, getter, letter) in enumerate(PANELS):
        row, col = divmod(idx, 4)
        ax = fig.add_subplot(4, 4, row * 4 + col + 1, projection=proj)
        field = getter(data)
        # Use explicit `colors=` (8 entries for 8 bins between 9 LEVELS).
        # Center bins are very-pale (#fddbc7 / #d1e5f0); outer bins shifted
        # one tier deeper than stock RdBu_r — see _BIN_COLORS at top of file.
        cf = ax.contourf(lons, lats, field, levels=LEVELS, colors=_BIN_COLORS,
                         extend='both', transform=ccrs.PlateCarree())
        cf.cmap.set_under('#021A40')
        cf.cmap.set_over('#400015')
        ax.coastlines(resolution='110m', linewidth=0.4, color='k')
        ax.set_global()
        ax.set_title(f'({letter}) {label}', fontsize=11, loc='left')
        gl = ax.gridlines(draw_labels=True, linewidth=0.3, alpha=0.4)
        gl.top_labels = False
        gl.right_labels = False
        gl.xlabel_style = {'size': 7}
        gl.ylabel_style = {'size': 7}

        dm = np.nanmean(field)
        ax.text(0.02, 0.02, f'mean={dm:+.2f}K',
                transform=ax.transAxes, fontsize=8,
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))

    # Colorbar
    cbar_ax = fig.add_axes([0.15, 0.04, 0.7, 0.015])
    cb = fig.colorbar(cf, cax=cbar_ax, orientation='horizontal',
                      ticks=LEVELS, extend='both')
    cb.set_label('Partial Temperature dT (K)', fontsize=11)

    fig.suptitle(f'pyCFRAM 13-panel surface dT decomposition: {case_name} (4×CO2)',
                 fontsize=14, y=0.97)

    plt.tight_layout(rect=[0, 0.07, 1, 0.95])
    fig_dir = cfg['_figures_dir']
    os.makedirs(fig_dir, exist_ok=True)
    suffix = 'pycfram-' + cfg.get('radiation', {}).get('scheme', 'unknown').lower()
    out_png = os.path.join(fig_dir, f'fig_13panel_global_{suffix}.png')
    plt.savefig(out_png, dpi=150, bbox_inches='tight')
    print('Saved:', out_png)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit('usage: plot_13panel_global.py <case_name>')
    main(sys.argv[1])
