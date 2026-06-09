#!/usr/bin/env python3
"""Global-mean vertical profile closure check.

Reads cases/<case>/output/cfram_result.nc, computes the area-weighted
(cos(lat)) global mean of every dT_* term, sums them per the CFRAM
closure identity

    dT_obs  ≈  Σ_radiative  +  (lhflx + shflx)  +  atmdyn  +  sfcdyn

and plots:
    (a) dT_observed   vs   Σ_all      vs   residual
    (b) Per-term profile decomposition (the contributing pieces)
    (c) Residual close-up + RMS by altitude band

Usage:
    python3 scripts/plot_closure_profile.py cesm2_4xco2_official
"""
import os
import sys
import argparse
import numpy as np
from netCDF4 import Dataset
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.config import load_case


# Closure identity — OLD CFRAM (Lu/Cai) convention:
#
#     dT_obs  =  Σ_rad_no_ts  +  dT_lhflx + dT_shflx  +  dT_atmdyn_resid + dT_ocndyn
#
# Key differences from pyCFRAM's native decomposition:
#   - NO dT_ts in the rad sum (OLD CFRAM doesn't decompose ts as a partial;
#     the ts effect is absorbed into atm_dyn).
#   - dT_atmdyn is computed as the TRUE RESIDUAL:
#       dT_atmdyn_resid = dT_obs − (Σ_rad_no_ts + dT_lhflx + dT_shflx + dT_ocndyn)
#     This guarantees Σ = dT_obs (closure exact, by construction).
#   - dT_ts remains in the per-panel display as a diagnostic but is NOT in Σ.
#
# Net effect: pyCFRAM's native dT_ts and frc_full-based dT_atmdyn are regrouped
# into a single OLD-style "atm_dyn residual". This reveals the true closure
# residual as the OLD-style atm_dyn value (compare with OLD CFRAM atm_dyn).
RAD_TERMS = ['co2', 'q', 'o3', 'solar', 'albedo', 'cloud', 'aerosol']  # NO ts (OLD style)
NONRAD_TERMS = ['lhflx', 'shflx']
DYN_TERMS = ['atmdyn', 'ocndyn']  # atmdyn is recomputed as residual below
DIAG_TERMS = ['ts']  # diagnostic only, NOT summed (kept for visualization)

# Per-term plotting colors
TERM_COLORS = {
    'co2':     '#d62728',
    'q':       '#1f77b4',
    'ts':      '#ff7f0e',
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


def area_weighted_mean(field_3d, lats):
    """Area-weighted (cos lat) global mean over (lat, lon) → 1D (lev,).

    field_3d shape: (lev, lat, lon).
    Ignores fill values (-999) by masking.
    """
    w = np.cos(np.deg2rad(lats))  # (lat,)
    w = w / w.sum() / field_3d.shape[2]  # normalise: integral = 1
    w2d = np.broadcast_to(w[:, None], field_3d.shape[1:])  # (lat, lon)
    mask = np.abs(field_3d) > 900.0
    masked = np.where(mask, np.nan, field_3d)
    # nan-safe weighted mean
    num = np.nansum(masked * w2d[None, :, :], axis=(1, 2))
    denom = np.sum(~mask * w2d[None, :, :], axis=(1, 2))
    return np.where(denom > 0, num / denom, np.nan)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('case', help='case name (cases/<case>/)')
    parser.add_argument('--out', default=None, help='output PNG path')
    args = parser.parse_args()

    cfg = load_case(args.case)
    nc_path = os.path.join(cfg['_output_dir'], 'cfram_result.nc')
    if not os.path.exists(nc_path):
        sys.exit("Missing %s — run --step run first." % nc_path)

    ds = Dataset(nc_path)
    lev = np.array(ds['lev'][:])
    lats = np.array(ds['lat'][:])
    nlev_atm = len(lev) - 1   # last index is surface

    def get_profile(name):
        if name not in ds.variables:
            return None
        return area_weighted_mean(np.array(ds[name][:]), lats)

    # 1D (nlev_atm+1,) profiles for each term
    dT_obs = get_profile('dT_observed')
    if dT_obs is None:
        sys.exit("dT_observed missing from %s" % nc_path)

    profiles = {}
    for t in RAD_TERMS + NONRAD_TERMS + DYN_TERMS + DIAG_TERMS:
        p = get_profile('dT_' + t)
        if p is None:
            print("  WARN: dT_%s not in NC, treating as 0" % t)
            p = np.zeros_like(dT_obs)
        profiles[t] = p

    # Stash pyCFRAM's native dT_atmdyn for diagnostic comparison; then OVERRIDE
    # profiles['atmdyn'] with the OLD CFRAM-style residual so closure is exact.
    profiles['atmdyn_native'] = profiles['atmdyn'].copy()
    sum_rad = sum(profiles[t] for t in RAD_TERMS)            # no ts
    sum_nonrad = sum(profiles[t] for t in NONRAD_TERMS)
    profiles['atmdyn'] = (dT_obs - sum_rad - sum_nonrad
                          - profiles['ocndyn'])               # OLD-style residual
    sum_dyn = profiles['atmdyn'] + profiles['ocndyn']
    sum_all = sum_rad + sum_nonrad + sum_dyn
    residual = dT_obs - sum_all                               # ≡ 0 by construction

    # Scalar summary line
    plev_atm = lev[:nlev_atm]
    lower = (plev_atm > 700)
    midtrop = (plev_atm >= 450) & (plev_atm <= 650)
    upper = (plev_atm >= 200) & (plev_atm <= 400)
    strato = (plev_atm < 200) & (plev_atm > 0)

    def rms_band(arr, mask):
        a = arr[:nlev_atm][mask]
        a = a[np.isfinite(a)]
        return float(np.sqrt(np.mean(a**2))) if len(a) else np.nan

    print("=" * 78)
    print("Global area-weighted closure check: %s" % args.case)
    print("=" * 78)
    print(f"  dT_obs     surface = {dT_obs[-1]:+7.3f} K  ;  col-mean (atm) = {np.nanmean(dT_obs[:nlev_atm]):+7.3f} K")
    print(f"  Σ all      surface = {sum_all[-1]:+7.3f} K  ;  col-mean       = {np.nanmean(sum_all[:nlev_atm]):+7.3f} K")
    print(f"  residual   surface = {residual[-1]:+7.3f} K  ;  col-mean       = {np.nanmean(residual[:nlev_atm]):+7.3f} K")
    print()
    print(f"  RMS residual by band (K):")
    print(f"    lower (>700 hPa)     = {rms_band(residual, lower):.4f}")
    print(f"    mid-trop (450-650)   = {rms_band(residual, midtrop):.4f}")
    print(f"    upper-trop (200-400) = {rms_band(residual, upper):.4f}")
    print(f"    strato (<200 hPa)    = {rms_band(residual, strato):.4f}")

    atm = slice(0, nlev_atm)
    # invert y axis (pressure decreasing upward); surface at bottom
    y = plev_atm
    y_lim = (y.max() + 20, y.min() - 5)

    # Per-term panel list, then filter: drop terms where max|profile| < 0.1 K
    # (saves space for visually empty terms — typically o3/solar/aerosol when
    # base==warm). Closure sum still includes all terms; this is display-only.
    # DIAG_TERMS (ts) is intentionally NOT in the panel list so the figure
    # mirrors OLD CFRAM's layout exactly (panels c..j = co2, q, albedo, cloud,
    # lhflx, shflx, atmdyn, ocndyn — no ts, since OLD CFRAM has no ts partial).
    candidate_terms = RAD_TERMS + NONRAD_TERMS + DYN_TERMS
    all_terms = [t for t in candidate_terms
                 if np.nanmax(np.abs(profiles[t])) >= 0.1]
    dropped = [t for t in candidate_terms if t not in all_terms]
    if dropped:
        print(f"  panels dropped (max|prof| < 0.1 K): {dropped}")

    # Layout: row 0 = closure overview (a/b spans full width 2 cols each);
    # rows 1+ = per-term panels (4 cols, ceil(N/4) rows).
    ncols = 4
    nrows_term = (len(all_terms) + ncols - 1) // ncols
    nrows_total = 1 + nrows_term
    fig_h = 5 + 4 * nrows_term  # ~4 inches per term row
    fig = plt.figure(figsize=(16, fig_h))
    height_ratios = [1.25] + [1.0] * nrows_term
    gs = fig.add_gridspec(nrows_total, ncols, hspace=0.42, wspace=0.18,
                          height_ratios=height_ratios)

    # Panel (a): closure check (spans cols 0-1)
    axA = fig.add_subplot(gs[0, 0:2])
    axA.plot(dT_obs[atm], y, 'k-', lw=2.0, label='dT_observed', zorder=4)
    axA.plot(sum_all[atm], y, '--', color='#d62728', lw=1.8,
             label=r'Σ (rad_no_ts + nonrad + atmdyn_resid + ocndyn)', zorder=3)
    axA.plot(residual[atm], y, ':', color='#2ca02c', lw=1.8,
             label='residual = obs − Σ  (≈0 by construction)', zorder=2)
    axA.axvline(0, color='gray', lw=0.5, ls='--')
    axA.set_title('(a) OLD-CFRAM-style closure (atmdyn = residual)', fontsize=11)
    axA.set_xlabel(r'$\Delta T$ (K)')
    axA.set_ylabel('Pressure (hPa)')
    axA.set_ylim(*y_lim)
    axA.legend(loc='best', fontsize=9, framealpha=0.92)
    axA.grid(True, alpha=0.3)

    # Panel (b): residual close-up (spans cols 2-3) — renamed from (c)
    axB = fig.add_subplot(gs[0, 2:4])
    axB.plot(residual[atm], y, '-', color='#2ca02c', lw=1.8,
             label='residual = obs − Σ')
    axB.axvline(0, color='gray', lw=0.5, ls='--')
    title = ('(b) Residual close-up\n'
             f'RMS: lower {rms_band(residual, lower):.3f} | '
             f'mid {rms_band(residual, midtrop):.3f} | '
             f'upper {rms_band(residual, upper):.3f} | '
             f'strato {rms_band(residual, strato):.3f} K')
    axB.set_title(title, fontsize=10)
    axB.set_xlabel(r'$\Delta T$ residual (K)')
    axB.set_ylim(*y_lim)
    axB.grid(True, alpha=0.3)
    rmin, rmax = np.nanmin(residual[atm]), np.nanmax(residual[atm])
    pad = max(0.1, (rmax - rmin) * 0.2)
    axB.set_xlim(rmin - pad, rmax + pad)

    # Per-term panels (filtered): rows 1+, dynamic grid
    panel_letters = [chr(ord('c') + i) for i in range(len(all_terms))]
    for idx, t in enumerate(all_terms):
        r = 1 + idx // ncols
        c = idx % ncols
        ax = fig.add_subplot(gs[r, c])
        prof = profiles[t][atm]
        col = TERM_COLORS.get(t, 'gray')
        ax.plot(prof, y, color=col, lw=1.6)
        ax.fill_betweenx(y, 0, prof, color=col, alpha=0.15)
        # For atmdyn panel: also overlay pyCFRAM's native (frc_full-based) value
        # for diagnostic comparison vs the OLD-style residual.
        if t == 'atmdyn':
            ax.plot(profiles['atmdyn_native'][atm], y,
                    color='#888', lw=1.0, ls='--', label='pyCFRAM native')
            ax.legend(loc='best', fontsize=7, framealpha=0.85)
        ax.axvline(0, color='gray', lw=0.5, ls='--')
        col_mean = np.nanmean(prof)
        sfc_val = profiles[t][-1]  # surface value (last lev index)
        suffix = ' (resid)' if t == 'atmdyn' else ''
        ax.set_title(f'({panel_letters[idx]}) {t}{suffix}  sfc={sfc_val:+.2f}K  col̄={col_mean:+.2f}K',
                     fontsize=9.0)
        ax.set_ylim(*y_lim)
        ax.grid(True, alpha=0.3)
        if c == 0:
            ax.set_ylabel('Pressure (hPa)')
        else:
            ax.tick_params(labelleft=False)
        if r == nrows_total - 1:
            ax.set_xlabel(r'$\Delta T$ (K)')

    case_name = cfg.get('case_name', args.case)
    desc = cfg.get('description', '')
    fig.suptitle(f'pyCFRAM closure: {case_name}  ({desc})', fontsize=12)
    # tight_layout conflicts with gridspec; use subplots_adjust instead
    fig.subplots_adjust(top=0.95)

    fig_dir = cfg['_figures_dir']
    os.makedirs(fig_dir, exist_ok=True)
    suffix = 'pycfram-' + cfg.get('radiation', {}).get('scheme', 'unknown').lower()
    out_path = args.out or os.path.join(fig_dir, f'fig_closure_profile_{suffix}.png')
    plt.savefig(out_path, dpi=140, bbox_inches='tight')
    print(f'\nSaved: {out_path}')
    ds.close()


if __name__ == '__main__':
    main()
