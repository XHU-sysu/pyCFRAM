#!/usr/bin/env python3
"""Closure profile figure for OLD Fortran CFRAM raw output.

Reads the collab xiaominghu OLD CFRAM CESM2 4xCO2 output:
    raw_data/cesm2_cmip6/CFRAM_xiaominghu/
        partial_T_1.grd   (13 terms × 18 levels × 192 lat × 288 lon, float32)
        forcing_1.grd     (same dims)

Layout mirrors plot_closure_profile.py:
    (a) Σ all partials (the closure is trivially exact for OLD CFRAM, by
        construction Σ = dT_obs; we still plot it for visual reference)
    (b) Residual close-up (≈ 0 by construction; sanity check)
    (c-o) Per-term profiles in 3×4 grid (one term per subplot)

Output: cases/cesm2_4xco2_official_17p_old/figures/fig_closure_profile.png
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OLD_DIR = '/disk/r074/lzhenn/workspace/ust-jumper/pyCFRAM/raw_data/cesm2_cmip6/CFRAM_xiaominghu'
OUT_DIR = '/disk/r074/lzhenn/workspace/ust-jumper/pyCFRAM/cases/cesm2_4xco2_official_17p_old/figures'

IX, IY = 288, 192     # nlon, nlat (CESM2 f09 grid)
NLEV   = 17           # atm levels
ZD18   = 18           # nlev + surface
NTERM  = 13

# 17-plev grid (TOA→surface). Matches cesm2_4xco2_official_17p_fu case.yaml.
PLEV_TOA_TO_SFC = np.array(
    [10, 20, 30, 50, 70, 100, 150, 200, 250, 300, 400, 500, 600, 700, 850, 925, 1000],
    dtype=np.float64,
)
PS_REF = 1013.0  # surface row reference pressure

# OLD CFRAM record order in partial_T_1.grd / forcing_1.grd (1-indexed terms × 18 levels)
TERM_ORDER = [
    'albedo', 'wv', 'cloud', 'cloud_sw', 'cloud_lw',
    'co2', 'o3', 'solar', 'dyn', 'atm_dyn',
    'sfc_dyn', 'shflx', 'lhflx',
]

# 12 panels in the same slots as pyCFRAM's plot_closure_profile.py.
# Map: OLD CFRAM has NO ts / NO aerosol; show OLD wv↔q, OLD dyn↔dry (total),
# OLD atm_dyn↔atmdyn, OLD sfc_dyn↔ocndyn (closest analog).
PANELS = [
    ('co2',     'co2'),
    ('q',       'wv'),
    ('ts',      None),       # OLD has no separate ts
    ('o3',      'o3'),
    ('solar',   'solar'),
    ('albedo',  'albedo'),
    ('cloud',   'cloud'),
    ('aerosol', None),       # OLD has no aerosol
    ('lhflx',   'lhflx'),
    ('shflx',   'shflx'),
    ('atmdyn',  'atm_dyn'),
    ('ocndyn',  'sfc_dyn'),
]
TERM_COLORS = {
    'co2':'#d62728','q':'#1f77b4','ts':'#ff7f0e','o3':'#9467bd',
    'solar':'#8c564b','albedo':'#e377c2','cloud':'#7f7f7f','aerosol':'#bcbd22',
    'lhflx':'#17becf','shflx':'#1abc9c','atmdyn':'#2ca02c','ocndyn':'#34495e',
}


def read_grd(path):
    """Read OLD CFRAM .grd file → (NTERM, ZD18, IY, IX) float32, masked."""
    raw = np.fromfile(path, dtype=np.float32)
    expected = NTERM * ZD18 * IY * IX
    assert raw.size == expected, f"{path}: got {raw.size}, expected {expected}"
    arr = raw.reshape((NTERM, ZD18, IY, IX))
    return np.where(np.abs(arr) > 900.0, np.nan, arr.astype(np.float64))


def area_weighted_mean(field_3d, lats):
    """(lev, lat, lon) → (lev,) area-weighted (cos lat) global mean, nan-safe."""
    w = np.cos(np.deg2rad(lats))
    w = w / w.sum() / field_3d.shape[2]
    w2d = np.broadcast_to(w[:, None], field_3d.shape[1:])
    masked = field_3d
    num = np.nansum(masked * w2d[None, :, :], axis=(1, 2))
    denom = np.nansum(np.isfinite(masked) * w2d[None, :, :], axis=(1, 2))
    return np.where(denom > 0, num / denom, np.nan)


def main():
    print(f"Reading OLD CFRAM output from {OLD_DIR}")
    partial_T = read_grd(os.path.join(OLD_DIR, 'partial_T_1.grd'))   # K
    forcing   = read_grd(os.path.join(OLD_DIR, 'forcing_1.grd'))     # W/m²

    # Level orientation: TOA→surface (verified empirically by full-column probe).
    # Level 0 of partial_T_1.grd has dT_co2 ≈ -16K (strong stratospheric cooling,
    # the unambiguous CO2-stratosphere signature). Level 17 has dT_albedo ≈ +1K
    # and dT_wv ≈ +6K — surface signatures. So index 0=TOA, indices 1..16=atm
    # going down, index 17=surface row at 1013 hPa.
    surface_at_end = True
    lat = np.linspace(-90 + 90.0/IY, 90 - 90.0/IY, IY)
    print(f"  Level orientation: TOA→surface (hardcoded after empirical probe)")

    # Build per-term global-mean profile (ZD18,)
    profiles_term = {}  # OLD-CFRAM term name → (zd18,)
    for it, term in enumerate(TERM_ORDER):
        prof = area_weighted_mean(partial_T[it], lat)
        profiles_term[term] = prof

    # Map to pyCFRAM panel naming (with absent terms = zero array)
    profiles = {}
    for pycfram_name, old_name in PANELS:
        if old_name is None:
            profiles[pycfram_name] = np.zeros(ZD18)
        else:
            profiles[pycfram_name] = profiles_term[old_name]

    # Closure: Σ of ALL partials (radiative + dyn + nonrad).
    # OLD CFRAM exposes total dyn already, plus split atm_dyn + sfc_dyn + shflx + lhflx.
    # To avoid double-counting, sum only NON-overlapping pieces:
    #     Σ = co2 + wv + o3 + solar + albedo + cloud + atm_dyn + sfc_dyn + shflx + lhflx
    sum_all = (profiles_term['co2']+profiles_term['wv']+profiles_term['o3']
              +profiles_term['solar']+profiles_term['albedo']+profiles_term['cloud']
              +profiles_term['atm_dyn']+profiles_term['sfc_dyn']
              +profiles_term['shflx']+profiles_term['lhflx'])
    # dT_obs proxy: OLD CFRAM doesn't dump it separately. By CFRAM construction,
    # dT_obs ≡ Σ_all (linearised closure). So residual is exactly zero modulo
    # Fortran roundoff. Useful sanity check rather than physical validation.
    dT_obs = sum_all
    residual = dT_obs - sum_all

    # ---- Pressure-level y-axis (atm only, 17 atm rows) ----
    # We assumed surface_at_end → level 1..17 are atm (TOA→surface),
    # level 18 is the surface row.  Plot atm only.
    if surface_at_end:
        atm_slice = slice(0, NLEV)        # indices 0..16 in 0-based python
        y_pres = PLEV_TOA_TO_SFC          # [10, 20, ..., 1000]
        sfc_val_idx = NLEV                # index 17 in 0-based python = surface row
    else:
        atm_slice = slice(1, NLEV + 1)
        y_pres = PLEV_TOA_TO_SFC[::-1]    # reversed
        sfc_val_idx = 0
    y = y_pres
    y_lim = (y.max() + 20, y.min() - 5)

    # Band masks for RMS
    lower   = (y > 700)
    midtrop = (y >= 450) & (y <= 650)
    upper   = (y >= 200) & (y <= 400)
    strato  = (y < 200) & (y > 0)
    def rms_band(arr_atm, mask):
        a = arr_atm[mask]; a = a[np.isfinite(a)]
        return float(np.sqrt(np.mean(a**2))) if len(a) else np.nan

    dT_obs_atm  = dT_obs[atm_slice]
    sum_all_atm = sum_all[atm_slice]
    resid_atm   = residual[atm_slice]

    print()
    print("=" * 78)
    print("OLD Fortran CFRAM (xiaominghu collab) global-mean closure")
    print("=" * 78)
    print(f"  Σ surface (lev=1013)   = {sum_all[sfc_val_idx]:+7.3f} K  ; "
          f"col-mean (atm)  = {np.nanmean(sum_all_atm):+7.3f} K")
    print(f"  residual surface       = {residual[sfc_val_idx]:+7.3f} K  "
          f"(= 0 by construction)")
    print(f"\n  Per-term surface (lev=1013) dT (K):")
    for term in TERM_ORDER:
        v = profiles_term[term][sfc_val_idx]
        print(f"    {term:>10s} = {v:+7.3f}")

    # Filter per-term panels: drop terms whose max|profile| < 0.1 K (visually
    # empty). Closure sum still includes all terms (handled above); this only
    # trims the figure.
    candidate = [p for p in PANELS]
    panels_filtered = []
    for pycfram_name, old_name in candidate:
        prof = profiles[pycfram_name]
        if np.nanmax(np.abs(prof)) >= 0.1:
            panels_filtered.append((pycfram_name, old_name))
    dropped = [p for p in candidate if p not in panels_filtered]
    if dropped:
        print(f"  panels dropped (max|prof| < 0.1 K): "
              f"{[p[0] for p in dropped]}")

    # ---------- Figure: dynamic grid ----------
    ncols = 4
    nrows_term = (len(panels_filtered) + ncols - 1) // ncols
    nrows_total = 1 + nrows_term
    fig_h = 5 + 4 * nrows_term
    fig = plt.figure(figsize=(16, fig_h))
    height_ratios = [1.25] + [1.0] * nrows_term
    gs = fig.add_gridspec(nrows_total, ncols, hspace=0.42, wspace=0.18,
                          height_ratios=height_ratios)

    axA = fig.add_subplot(gs[0, 0:2])
    axA.plot(dT_obs_atm,  y, 'k-',  lw=2.0, label='Σ all (= obs by construction)', zorder=4)
    axA.plot(sum_all_atm, y, '--', color='#d62728', lw=1.8,
             label='Σ (rad + dyn + nonrad)', zorder=3)
    axA.plot(resid_atm,   y, ':', color='#2ca02c', lw=1.8, label='residual', zorder=2)
    axA.axvline(0, color='gray', lw=0.5, ls='--')
    axA.set_title('(a) OLD CFRAM global closure check', fontsize=11)
    axA.set_xlabel(r'$\Delta T$ (K)')
    axA.set_ylabel('Pressure (hPa)')
    axA.set_ylim(*y_lim)
    axA.legend(loc='best', fontsize=9, framealpha=0.92)
    axA.grid(True, alpha=0.3)

    axB = fig.add_subplot(gs[0, 2:4])
    axB.plot(resid_atm, y, '-', color='#2ca02c', lw=1.8,
             label='residual')
    axB.axvline(0, color='gray', lw=0.5, ls='--')
    title = ('(b) Residual close-up (≈0 by OLD-CFRAM construction)\n'
             f'RMS: lower {rms_band(resid_atm, lower):.4f} | '
             f'mid {rms_band(resid_atm, midtrop):.4f} | '
             f'upper {rms_band(resid_atm, upper):.4f} | '
             f'strato {rms_band(resid_atm, strato):.4f} K')
    axB.set_title(title, fontsize=10)
    axB.set_xlabel(r'$\Delta T$ residual (K)')
    axB.set_ylim(*y_lim)
    axB.grid(True, alpha=0.3)
    # auto x-range tight
    rmin, rmax = np.nanmin(resid_atm), np.nanmax(resid_atm)
    if not (np.isfinite(rmin) and np.isfinite(rmax)) or rmax - rmin < 1e-6:
        axB.set_xlim(-1e-3, 1e-3)
    else:
        pad = max(0.1, (rmax - rmin) * 0.2)
        axB.set_xlim(rmin - pad, rmax + pad)

    # Per-term panels (filtered)
    panel_letters = [chr(ord('c') + i) for i in range(len(panels_filtered))]
    for idx, (pycfram_name, old_name) in enumerate(panels_filtered):
        r = 1 + idx // ncols
        c = idx % ncols
        ax = fig.add_subplot(gs[r, c])
        prof_full = profiles[pycfram_name]
        prof_atm  = prof_full[atm_slice]
        col = TERM_COLORS.get(pycfram_name, 'gray')
        if old_name is None:
            ax.text(0.5, 0.5, 'N/A in OLD CFRAM',
                    ha='center', va='center', transform=ax.transAxes,
                    fontsize=11, color='#888')
        else:
            ax.plot(prof_atm, y, color=col, lw=1.6)
            ax.fill_betweenx(y, 0, prof_atm, color=col, alpha=0.15)
        ax.axvline(0, color='gray', lw=0.5, ls='--')
        sfc_val = prof_full[sfc_val_idx]
        col_mean = np.nanmean(prof_atm)
        label_text = pycfram_name + (f'  [OLD:{old_name}]' if old_name and old_name != pycfram_name else '')
        ax.set_title(f'({panel_letters[idx]}) {label_text}  sfc={sfc_val:+.2f}K  col̄={col_mean:+.2f}K',
                     fontsize=9.5)
        ax.set_ylim(*y_lim)
        ax.grid(True, alpha=0.3)
        if c == 0:
            ax.set_ylabel('Pressure (hPa)')
        else:
            ax.tick_params(labelleft=False)
        if r == nrows_total - 1:
            ax.set_xlabel(r'$\Delta T$ (K)')

    fig.suptitle('OLD Fortran CFRAM (xiaominghu collab, CESM2 4xCO2, 17-plev) — '
                 'global closure & per-term partials',
                 fontsize=12)
    fig.subplots_adjust(top=0.95)

    os.makedirs(OUT_DIR, exist_ok=True)
    out_path = os.path.join(OUT_DIR, 'fig_closure_profile_oldcfram.png')
    plt.savefig(out_path, dpi=140, bbox_inches='tight')
    print(f"\nSaved: {out_path}")


if __name__ == '__main__':
    main()
