#!/usr/bin/env python3
"""Workaround for shapely 2.x / cartopy 0.22 MultiPolygon bug on hqlx204.

Plot 13-panel north polar surface dT decomposition from a pre-extracted
2D surface slice NC (~11 MB) instead of the full 380 MB cfram_result.nc.

Usage:
    python3 scripts/plot_13panel_polar_from_slice.py /tmp/cesm2_rrtmg_sfc.nc <out_path>
"""
import sys
import os
import numpy as np
from netCDF4 import Dataset
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeat

LEVELS_DT = [-15, -7, -3, -1, -0.5, 0.5, 1, 3, 7, 15]
CMAP = plt.get_cmap('RdBu_r', len(LEVELS_DT)-1)

# 13 panels matching plot_13panel_polar.py layout
PANELS = [
    ('AL',   'dT_albedo'),
    ('WV',   'dT_q'),
    ('CLD',  'dT_cloud'),
    ('CLDS', 'dT_cloud_sw'),
    ('CLDL', 'dT_cloud_lw'),
    ('CO2',  'dT_co2'),
    ('O3',   'dT_o3'),
    ('SR',   'dT_solar'),
    ('DYN',  'dT_dry'),
    ('ATM',  'dT_atmdyn'),
    ('OCH',  'dT_ocndyn'),
    ('SH',   'dT_shflx'),
    ('LH',   'dT_lhflx'),
]


def main():
    if len(sys.argv) < 3:
        sys.exit('usage: plot_13panel_polar_from_slice.py <slice.nc> <out.png>')
    slice_path = sys.argv[1]
    out_path = sys.argv[2]
    case_label = os.path.basename(slice_path).replace('_sfc.nc', '').replace('.nc', '')

    ds = Dataset(slice_path)
    lats = np.array(ds['lat'][:])
    lons = np.array(ds['lon'][:])

    norm = mcolors.BoundaryNorm(LEVELS_DT, CMAP.N, clip=True)

    # 4 columns × 4 rows layout, last row first panel only
    fig = plt.figure(figsize=(15, 16))
    for idx, (label, key) in enumerate(PANELS):
        ax = fig.add_subplot(4, 4, idx+1,
                             projection=ccrs.NorthPolarStereo())
        ax.set_extent([-180, 180, 60, 90], crs=ccrs.PlateCarree())
        try:
            data = np.array(ds[key][:])
        except KeyError:
            ax.set_title(f'({chr(97+idx)}) {label}  (MISSING)', fontsize=10)
            ax.coastlines(linewidth=0.4)
            continue
        data = np.clip(np.where(np.abs(data) > 900, np.nan, data), -20, 20)
        # Use pcolormesh (workaround for shapely contourf bug)
        # Need to add cyclic point for proper wrap-around on polar projection
        lon2d, lat2d = np.meshgrid(lons, lats)
        # Append last lon = first lon + 360 to close
        lon_w = np.concatenate([lons, [lons[0] + 360]])
        data_w = np.concatenate([data, data[:, :1]], axis=1)
        lon2d_w, lat2d_w = np.meshgrid(lon_w, lats)
        im = ax.pcolormesh(lon2d_w, lat2d_w, data_w, transform=ccrs.PlateCarree(),
                           norm=norm, cmap=CMAP, shading='auto')
        ax.coastlines(linewidth=0.4)
        ax.gridlines(linewidth=0.3, alpha=0.4)
        ax.set_title(f'({chr(97+idx)}) {label}', fontsize=10)

    fig.suptitle(f'pyCFRAM 13-panel surface dT decomposition: {case_label}\n'
                 f'North polar stereographic, |lat| ≥ 60°  (workaround: pcolormesh)',
                 fontsize=12)
    plt.subplots_adjust(bottom=0.08, top=0.94, left=0.05, right=0.95,
                        hspace=0.25, wspace=0.1)
    cbar_ax = fig.add_axes([0.15, 0.04, 0.7, 0.012])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=CMAP), cax=cbar_ax,
                      orientation='horizontal', ticks=LEVELS_DT)
    cb.set_label('Partial Temperature dT (K)', fontsize=10)

    fig.savefig(out_path, dpi=140, bbox_inches='tight')
    print(f'Saved: {out_path}')
    ds.close()


if __name__ == '__main__':
    main()
