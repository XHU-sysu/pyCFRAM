#!/usr/bin/env python3
"""Region-averaged vertical CFRAM decomposition profile.

Same six terms as the surface figure (plot_alaska_decomposition.py):
  Water Vapor (dT_q), Cloud (dT_cloud), Aerosol (dT_aerosol),
  Surface Process (dT_sfcdyn — already = ocndyn+lhflx+shflx, no double count),
  Atm. Dynamics (dT_atmdyn), Total (dT_observed).

Plots region-averaged dT(level) vs pressure. The bottom point is the SURFACE
SKIN row (~1013 hPa), so the bottom of each curve matches the surface map; the
next point up is the 1000 hPa AIR layer (these can differ strongly — surface
skin vs boundary-layer air decoupling).

Usage: python3 scripts/plot_alaska_profile.py <case> [region_key] [land|ocean|all]
  region_key default = main_fire_region; surface mask default = land.
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from netCDF4 import Dataset

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.config import load_case

TERMS = [
    ('Water Vapor',    ['dT_q'],        'tab:green'),
    ('Cloud',          ['dT_cloud'],    'tab:blue'),
    ('Aerosol',        ['dT_aerosol'],  'tab:brown'),
    # sfcdyn already = ocndyn+lhflx+shflx (verified ~1e-14); use it alone.
    ('Surface Process',['dT_sfcdyn'],   'tab:red'),
    ('Atm. Dynamics',  ['dT_atmdyn'],   'tab:purple'),
    ('Total',          ['dT_observed'], 'black'),
]


def land_mask(lat, lon, yi, xi):
    """Boolean land mask over the (yi,xi) region using Natural Earth polygons."""
    import cartopy.io.shapereader as shp
    from shapely.geometry import Point
    from shapely.prepared import prep
    from shapely.ops import unary_union
    L = prep(unary_union(list(shp.Reader(shp.natural_earth(
        resolution='110m', category='physical', name='land')).geometries())))
    lon2 = np.where(lon > 180, lon - 360, lon)
    mask = np.zeros((len(lat), len(lon)), dtype=bool)
    for i in yi:
        for j in xi:
            mask[i, j] = L.contains(Point(lon2[j], lat[i]))
    return mask


def main():
    case = sys.argv[1] if len(sys.argv) > 1 else 'alaska_wf22_peak'
    region_key = sys.argv[2] if len(sys.argv) > 2 else 'main_fire_region'
    surf = sys.argv[3] if len(sys.argv) > 3 else 'land'
    cfg = load_case(case)
    nc = Dataset(os.path.join(cfg['_output_dir'], 'cfram_result.nc'))
    lat = np.array(nc.variables['lat'][:])
    lon = np.array(nc.variables['lon'][:])
    lev = np.array(nc.variables['lev'][:])      # last entry = surface skin

    reg = cfg['plot'].get(region_key) or cfg['plot']['map_extent']
    (lo0, lo1), (la0, la1) = reg['lon'], reg['lat']
    yi = np.where((lat >= la0) & (lat <= la1))[0]
    xi = np.where((lon >= lo0) & (lon <= lo1))[0]

    mask = np.zeros((len(lat), len(lon)), dtype=bool)
    mask[np.ix_(yi, xi)] = True
    if surf in ('land', 'ocean'):
        lm = land_mask(lat, lon, yi, xi)
        mask = mask & (lm if surf == 'land' else ~lm)
    print('%s region (%s): %d cells' % (region_key, surf, int(mask.sum())))

    # Per-level above-surface mask. CFRAM terms are filled with -999 below the
    # surface (where plev > ps), whereas dT_observed carries ERA5's below-ground
    # EXTRAPOLATED temperature (a finite, non-physical value). Build the
    # above-surface mask from a CFRAM term and apply it to ALL variables
    # (including observed) so that, at every level, all curves are averaged over
    # the SAME above-ground cells: this drops observed's below-ground
    # extrapolation and keeps the terms + Total comparable (level-by-level
    # closure). Cell count still shrinks toward the surface (terrain) — reported
    # below; for a fixed population across levels use only cells valid at 1000 hPa.
    refname = 'dT_sfcdyn' if 'dT_sfcdyn' in nc.variables else 'dT_q'
    ref = np.array(nc.variables[refname][:], dtype=np.float64)
    abovesurf = np.abs(ref) <= 900   # (lev, lat, lon): True = above surface
    nvalid = np.array([int((mask & abovesurf[k]).sum()) for k in range(len(lev))])
    print('above-surface cells per level: min=%d (near surface), max=%d (free trop), skin=%d'
          % (int(nvalid[:-1].min()), int(nvalid[:-1].max()), int(nvalid[-1])))

    def prof(vname):
        a = np.array(nc.variables[vname][:], dtype=np.float64)
        a = np.where(np.abs(a) > 900, np.nan, a)
        out = np.full(len(lev), np.nan)
        for k in range(len(lev)):
            m = mask & abovesurf[k]
            if m.any():
                out[k] = np.nanmean(a[k][m])
        return out

    natm = len(lev) - 1
    # pressure axis: surface skin (bottom) then atmosphere 1000->TOA
    p_full = np.concatenate([[lev[-1]], lev[:natm]])

    fig, ax = plt.subplots(figsize=(6.5, 7.8))
    for label, vnames, color in TERMS:
        pr = np.zeros(len(lev))
        for v in vnames:
            if v in nc.variables:
                pr = pr + np.nan_to_num(prof(v))
        prof_full = np.concatenate([[pr[-1]], pr[:natm]])   # skin first
        ax.plot(prof_full, p_full, color=color, lw=2.0, marker='o', ms=3, label=label)

    ax.axhline(lev[-1], color='0.7', lw=0.6, ls=':')   # surface skin level marker
    ax.text(ax.get_xlim()[0], lev[-1], ' surface skin', va='bottom', fontsize=7, color='0.4')
    ax.axvline(0, color='0.5', lw=0.8)
    ax.set_yscale('log')
    ax.set_ylim(1025, 50)
    ax.set_yticks([1013, 1000, 850, 700, 500, 300, 200, 100, 50])
    ax.set_yticklabels(['skin', 1000, 850, 700, 500, 300, 200, 100, 50])
    ax.set_xlabel('Partial temperature contribution (K)')
    ax.set_ylabel('Pressure (hPa)')
    ax.set_title('%s: %s-region (%s) CFRAM vertical profile\n'
                 'bottom point = surface skin (matches the map)'
                 % (cfg.get('case_name', case), region_key, surf), fontsize=10)
    ax.legend(fontsize=9, loc='best')
    ax.grid(alpha=0.3)
    nc.close()

    tag = '%s_%s' % (region_key, surf)
    out = os.path.join(cfg['_figures_dir'], 'fig_profile_%s.png' % tag)
    fig.savefig(out, dpi=140, bbox_inches='tight')
    print('Saved:', out)


if __name__ == '__main__':
    main()
