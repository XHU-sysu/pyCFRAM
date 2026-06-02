#!/usr/bin/env python3
"""Plot 13-panel surface dT decomposition on north polar stereographic projection.

Mirrors raw/north.jpg layout: 4 cols × 4 rows (last row has only LH panel).
Uses same color levels and panel order as plot_13panel_global.py.

Usage:
    python3 scripts/plot_13panel_polar.py <case_name>
"""
import os, sys
import numpy as np
from netCDF4 import Dataset
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.path as mpath
import cartopy.crs as ccrs
import cartopy.feature as cfeature

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.config import load_case

# Match OLD CFRAM raw/north.jpg colorbar exactly. 9 boundaries → 8 internal bins.
# CRITICAL: the two center bins (-0.5, 0] and (0, 0.5] must be TRUE WHITE so
# panels with frc≈0 (e.g., dT_o3/dT_solar in cases where o3_b==o3_w) render
# blank instead of light-blue (which would happen with stock RdBu_r BoundaryNorm
# since RdBu_r at center index 0.375 is light blue, not white).
LEVELS = [-15.0, -7.0, -2.0, -0.5, 0.0, 0.5, 2.0, 7.0, 15.0]
# 8 internal bin colors. Each bin one tier DEEPER than previous default, so
# (7, 15] uses what used to be the extend color (very dark red), and the
# previous extend gets pushed to an even darker shade. (0, 0.5] and (-0.5, 0]
# use very pale tints (near-white but with sign info), per user spec.
_BIN_COLORS = [
    '#053061',   # (-15, -7]    deepest blue  (was extend < -15)
    '#2166ac',   # (-7, -2]     2nd deepest   (was (-15, -7])
    '#92c5de',   # (-2, -0.5]   maintain (3rd deepest)
    '#d1e5f0',   # (-0.5, 0]    very pale blue
    '#fddbc7',   # (0, 0.5]     very pale red
    '#f4a582',   # (0.5, 2]     maintain (3rd deepest)
    '#b2182b',   # (2, 7]       2nd deepest   (was (7, 15])
    '#67001f',   # (7, 15]      deepest red   (was extend > 15)
]
CMAP = mcolors.ListedColormap(_BIN_COLORS)
CMAP.set_under('#021A40')   # values < -15: even darker navy
CMAP.set_over('#400015')    # values > 15:  even darker red

# Same panel order as plot_13panel_global.py
PANELS = [
    ('AL',   'dT_albedo',   'a'),
    ('WV',   'dT_q',        'b'),
    ('CLD',  'dT_cloud',    'c'),
    ('CLDS', 'dT_cloud_sw', 'd'),
    ('CLDL', 'dT_cloud_lw', 'e'),
    ('CO2',  'dT_co2',      'f'),
    ('O3',   'dT_o3',       'g'),
    ('SR',   'dT_solar',    'h'),
    ('DYN',  'dT_dry',      'i'),
    ('ATM',  'dT_atmdyn',   'j'),
    ('OCH',  'dT_ocndyn',   'k'),
    ('SH',   'dT_shflx',    'l'),
    ('LH',   'dT_lhflx',    'm'),
]

# Latitude extent for polar view (degrees, abs value)
LAT_LIMIT = 60   # show region poleward of 60°N (north pole) or 60°S (south)


def load_surface(nc_path):
    """Load surface row of all panel variables. Returns dict with arrays + coords.

    Wraps longitude by appending lon[0]+360 column to eliminate the visible
    seam at the prime meridian on polar projections.
    """
    nc = Dataset(nc_path)
    lat = np.array(nc.variables['lat'][:])
    lon = np.array(nc.variables['lon'][:])
    # Append wrap column
    lon_w = np.append(lon, lon[0] + 360.0)
    out = {'lat': lat, 'lon': lon_w}
    for _, var, _ in PANELS:
        if var in nc.variables:
            arr = np.array(nc.variables[var][-1, :, :], dtype=np.float64)
            arr = np.where(np.abs(arr) > 900, np.nan, arr)
            arr_w = np.concatenate([arr, arr[:, :1]], axis=1)
        else:
            arr_w = np.full((len(lat), len(lon_w)), np.nan)
        out[var] = arr_w
    nc.close()
    return out


def make_circle_boundary():
    """Create a circular boundary path for polar stereographic axes."""
    theta = np.linspace(0, 2 * np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    return mpath.Path(verts * radius + center)


def plot_polar(data, case_name, hemisphere='north', out_path=None):
    """Plot 13-panel polar stereographic figure.

    hemisphere: 'north' (NorthPolarStereo, lat > 40°N) or
                'south' (SouthPolarStereo, lat < -40°S).
    """
    if hemisphere == 'north':
        proj = ccrs.NorthPolarStereo(central_longitude=0)
        lat_min, lat_max = LAT_LIMIT, 90
        title_suffix = 'North polar'
    elif hemisphere == 'south':
        proj = ccrs.SouthPolarStereo(central_longitude=0)
        lat_min, lat_max = -90, -LAT_LIMIT
        title_suffix = 'South polar'
    else:
        raise ValueError('hemisphere must be north or south')

    fig = plt.figure(figsize=(20, 22))
    circle = make_circle_boundary()
    cf = None

    for idx, (label, var, letter) in enumerate(PANELS):
        ax = fig.add_subplot(4, 4, idx + 1, projection=proj)
        field = data[var]
        # Restrict to hemisphere extent
        ax.set_extent([-180, 180, lat_min, lat_max], crs=ccrs.PlateCarree())
        ax.set_boundary(circle, transform=ax.transAxes)
        if not np.all(np.isnan(field)):
            # Use pcolormesh (NOT contourf): shapely 2.x + cartopy 0.22
            # MultiPolygon API has a latent bug that for SOME polar columns
            # silently drops contour bins, leaving large white patches even
            # when the data is correct (observed on RRTMG cesm2_4xco2_official_17p
            # panel (b) WV, where Arctic dT_q is uniformly +4..+8 K but contourf
            # rendered as huge white area). pcolormesh sidesteps the issue.
            norm = mcolors.BoundaryNorm(LEVELS, len(_BIN_COLORS), clip=False)
            cf = ax.pcolormesh(data['lon'], data['lat'], field,
                               norm=norm, cmap=CMAP, shading='auto',
                               transform=ccrs.PlateCarree())
        ax.coastlines(resolution='110m', linewidth=0.4, color='k')
        gl = ax.gridlines(draw_labels=False, linewidth=0.3, alpha=0.4)
        ax.set_title(f'({letter}) {label}', fontsize=12, loc='left')

    # Colorbar
    cbar_ax = fig.add_axes([0.15, 0.06, 0.7, 0.012])
    fig.colorbar(cf, cax=cbar_ax, orientation='horizontal',
                 ticks=LEVELS, extend='both',
                 label='Partial Temperature dT (K)')
    fig.suptitle(f'pyCFRAM 13-panel surface dT decomposition: {case_name} (4×CO2)\n'
                 f'{title_suffix} stereographic, |lat| ≥ {LAT_LIMIT}°',
                 fontsize=14, y=0.97)
    plt.tight_layout(rect=[0, 0.08, 1, 0.94])

    if out_path:
        plt.savefig(out_path, dpi=130, bbox_inches='tight')
        print('Saved:', out_path)
    plt.close(fig)


def main(case_name):
    cfg = load_case(case_name)
    nc_path = os.path.join(cfg['_output_dir'], 'cfram_result.nc')
    if not os.path.exists(nc_path):
        sys.exit('Missing %s' % nc_path)

    data = load_surface(nc_path)
    fig_dir = cfg['_figures_dir']
    os.makedirs(fig_dir, exist_ok=True)

    suffix = 'pycfram-' + cfg.get('radiation', {}).get('scheme', 'unknown').lower()
    # Both hemispheres
    plot_polar(data, case_name, 'north',
               os.path.join(fig_dir, f'fig_13panel_north_polar_{suffix}.png'))
    plot_polar(data, case_name, 'south',
               os.path.join(fig_dir, f'fig_13panel_south_polar_{suffix}.png'))


if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit('usage: plot_13panel_polar.py <case_name>')
    main(sys.argv[1])
