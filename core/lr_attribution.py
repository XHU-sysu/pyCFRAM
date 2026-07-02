"""M3 route i: per-physical-process Lapse-Rate attribution.

Reuses the M2 kernel machinery (core/kernels.py, core/lr_kernel.py) but
applies it to each CFRAM partial temperature response `dT_X(p)` in turn,
instead of the total `dT_observed`. Answers "which physical process
shaped the lapse-rate structure" by treating each `dT_X(p)` as if it were
its own mini temperature response and asking what TOA LW perturbation its
*non-uniform* part would produce (docs/plan.md §7.0 route i).

No new physics: 0 additional radiative-transfer or kernel code, purely a
loop over `core.lr_kernel.extract_and_apply` (docs/plan.md §7.0 "奥卡姆最
省").
"""
import numpy as np

from core.lr_kernel import extract_and_apply

# Terms available in cfram_result.nc that have a physically meaningful
# vertical dT_X(p) + skin dT_X profile (mirrors scripts/plot_closure_profile.py
# RAD_TERMS/NONRAD_TERMS/DYN_TERMS, minus 'ts' which has no vertical structure
# of its own to decompose into lapse-rate vs Planck).
DEFAULT_TERMS = ['q', 'co2', 'o3', 'solar', 'albedo', 'cloud', 'aerosol',
                  'lhflx', 'shflx', 'atmdyn', 'ocndyn']


def attribute_lr(dT_terms, ps_pa, plev_native_pa, lat_1d, kernelset,
                  month='annual', sky='all-sky'):
    """Per-process ΔR_LR / ΔR_PL attribution.

    Parameters
    ----------
    dT_terms : dict[str, (nk_native, lat, lon) K] -> atmospheric part of
        each dT_X, already split from the skin row (see
        scripts/compute_lr_attribution.py for how this is extracted from
        cfram_result.nc).
    dT_skin_terms : passed via a second dict argument is avoided here --
        callers pass `dT_terms[name]` as (atm, skin) tuples instead (see
        below) to keep this function single-purpose.

    Returns
    -------
    dict[str, dict] : term -> {'dR_lr': (lat,lon), 'dR_pl': (lat,lon)}
    """
    out = {}
    for name, (dT_atm, dT_skin) in dT_terms.items():
        out[name] = extract_and_apply(
            dT_atm, dT_skin, ps_pa, plev_native_pa, lat_1d, kernelset,
            month=month, sky=sky)
    return out


def additivity_residual(dR_lr_total, dR_lr_by_term):
    """ΔR_LR^total − Σ_X ΔR_LR^{(X)} (docs/plan.md WP-M3.1 Done / R7).

    Not expected to be zero -- CFRAM's first-order Taylor decomposition is
    linear in each dT_X but the kernel application (interpolation +
    integration) is itself linear, so the residual reflects genuine
    non-additivity already present in the CFRAM dT_X terms themselves
    (aerosol/cloud coupling etc., same character as the aerosol/cloud
    additivity residuals documented in session_log.md), not a bug in this
    module.
    """
    total_from_terms = np.zeros_like(dR_lr_total)
    for term_field in dR_lr_by_term.values():
        total_from_terms = total_from_terms + term_field
    return dR_lr_total - total_from_terms
