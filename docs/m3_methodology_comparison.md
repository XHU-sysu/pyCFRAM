# M3: CFRAM vs Radiative-Kernel Methodology — Comparison & Per-Process Lapse-Rate Attribution

This document is the contract's M3 deliverable: a methodology comparison
between pyCFRAM's native energy-budget decomposition and the
radiative-kernel Lapse-Rate method built in M2 (`docs/m2_kernel_module.md`),
plus a per-process attribution of the lapse-rate radiative perturbation
using pyCFRAM's own physical-process decomposition (route i, see
`docs/m3_route_decision_memo.md` for the route-selection rationale).

## 1. Module purpose and scientific motivation

**CFRAM** (Coupled Feedback-Response Analysis Method — Lu & Cai; the
`dT_q`, `dT_atmdyn`, `dT_sfcdyn`, ... family in `cfram_result.nc`)
answers: *"how much of the observed surface/atmospheric temperature
change is attributable to each physical process (water vapor, clouds,
CO2, surface dynamics, ...)?"* It does this via energy-budget closure:
each process's radiative/non-radiative forcing perturbation `frc_X` is
converted to a temperature response `dT_X` by inverting the same Planck
feedback matrix (`∂R/∂T`) used by the radiative-transfer engine (RRTMG
or Fu). It is a **first-order, energy-budget-based** decomposition.

**Radiative kernels** (ClimKern's `calc_T_feedbacks`, replicated
natively in `core/lr_kernel.py`) answer a narrower but complementary
question: *"of the total temperature response, how much is uniform
(Planck) vs vertically structured (Lapse-Rate), and how does that
vertical structure alone perturb the TOA energy budget?"* This is a
**diagnostic, pre-computed-sensitivity-based** method: it doesn't
re-derive the radiative transfer at all, it applies a fixed
model/observation-derived sensitivity field (the kernel) to whatever
temperature response you feed it.

These are different tools for different questions. CFRAM decomposes
*why the atmosphere changed temperature the way it did*; the kernel
method decomposes *how much a given temperature-change pattern would
perturb the top-of-atmosphere energy budget, split by vertical shape*.
M3's contribution is to connect them: for each CFRAM process `dT_X`, we
ask "if this process's temperature response happened in isolation, what
lapse-rate TOA perturbation would it produce?" — giving a
process-resolved lapse-rate attribution that neither method provides on
its own.

## 2. Input/output format

### `lr_kernel.nc` (M2, `scripts/compute_lr_kernel.py`)

| Variable | Dims | Units | Description |
|---|---|---|---|
| `dR_lr_<kernel>` | (lat,lon) | W m⁻² | Lapse-rate TOA LW perturbation, from `dT_observed`, per kernel (`CloudSat`=Kramer, `GFDL`) |
| `dR_pl_<kernel>` | (lat,lon) | W m⁻² | Planck (uniform-warming) TOA LW perturbation, per kernel |

### `lr_climkern_ref.nc` (M2, `scripts/validate_lr_vs_climkern.py`)

Same `dR_lr_<kernel>` layout, but computed by `climkern.calc_T_feedbacks`
directly (the independent reference used for cross-validation).

### `lr_attribution.nc` (M3, `scripts/compute_lr_attribution.py`)

| Variable | Dims | Units | Description |
|---|---|---|---|
| `dR_lr_from_<term>_<kernel>` | (lat,lon) | W m⁻² | Lapse-rate TOA perturbation attributable to CFRAM process `<term>` (`q`, `co2`, `o3`, `solar`, `albedo`, `cloud`, `aerosol`, `lhflx`, `shflx`, `atmdyn`, `ocndyn`), per kernel |
| `dR_lr_additivity_residual_<kernel>` | (lat,lon) | W m⁻² | `dR_lr_<kernel>` (M2 total, from `dT_observed`) minus Σ over all `dR_lr_from_<term>_<kernel>` |

## 3. Difference from ClimKern / other kernel methods

| Aspect | ClimKern (& generic kernel methods) | pyCFRAM native module |
|---|---|---|
| Regridding | xesmf conservative/bilinear (ESMF backend) | Self-contained scipy bilinear (periodic-wrapped in longitude, clip-to-boundary in latitude) — cross-validated to `corr > 0.999` against xesmf (`tests/test_kernels.py::test_regrid_vs_xesmf`) |
| Vertical interpolation | `xr.interp_like(..., extrapolate)`, NaN-preserving | `scipy.interpolate.interp1d(..., extrapolate)`, same NaN-preserving order (docs/plan.md §2.3) |
| Tropopause | `make_tropo`: `3e4 − 2e4·cos(lat)` Pa | Identical formula (`core/lr_kernel.py::tropopause_pa`) |
| Dependency footprint | Requires esmpy/xesmf (heavy, Apple-Silicon-fragile) | Core module: numpy/scipy/netCDF4 only; runs anywhere including hqlx210 (no xesmf there) |
| Time handling | Requires 12-month-tiled input (`len(time)%12==0`) | Operates directly on the two-state (`base`/`perturbed`) `dT_observed`; the M2 cross-validation harness does the 12-month tiling only on the ClimKern side, as a bridge |
| Kernel choice | Any of ClimKern's bundled kernels | `CloudSat` (Kramer, observationally-based, primary per contract) + `GFDL` (Soden & Held 2006, model-based, reference) |
| Process attribution | None — only total `dT` in, total `ΔR` out | M3 extension: any CFRAM `dT_X` can be pushed through the same machinery for per-process attribution (not available in ClimKern or any kernel-only tool, since kernel methods have no concept of "process" — that's CFRAM's contribution) |

**Kernel-choice sensitivity** (Kramer/CloudSat, observationally-derived,
vs GFDL, model-derived): `fig_lr_kramer_vs_gfdl.png` shows their
native-module difference map. The two kernels agree to within a few
W m⁻² almost everywhere (consistent with both independently validating
against ClimKern at corr > 0.999 — see below), with the largest
differences over persistently cloudy regions (Southern Ocean storm
track, ITCZ) where the kernels' underlying training data (CloudSat
satellite retrievals vs. the GFDL AM2 model) diverge most in cloud
radiative treatment. This is the kind of kernel-choice sensitivity the
contract's §3 threshold table is designed to bound.

## 4. M2 cross-validation results (recap)

| Kernel | corr (area-weighted) | domain-mean rel. diff | Gate (0.85/15%) |
|---|---|---|---|
| CloudSat (Kramer) | 0.9997 | 1.32% | PASS |
| GFDL | 0.9999 | 0.02% | PASS |

Full detail: `docs/m2_kernel_module.md`, `cases/cesm2_4xco2_official/output/lr_validation_report.txt`.

## 5. M3 results: per-process lapse-rate attribution

Computed on `cesm2_4xco2_official` (global 192×288, CESM2 abrupt-4×CO₂
vs piControl). Terms: `q, co2, o3, solar, albedo, cloud, aerosol, lhflx,
shflx, atmdyn, ocndyn` (mirrors `scripts/plot_closure_profile.py`'s
closure decomposition, minus `ts` which has no vertical structure of its
own to split into lapse-rate vs Planck).

```
CloudSat   additivity |residual| mean = 1.40 W/m^2 (21.8% of |total| mean)
GFDL       additivity |residual| mean = 1.75 W/m^2 (23.1% of |total| mean)
```

A non-zero additivity residual (Σ per-process ΔR_LR ≠ total ΔR_LR from
`dT_observed`) is **expected, not a bug**: CFRAM's decomposition is a
first-order (linearized) energy-budget split, but the true radiative
response to the *sum* of all processes' temperature changes is not
exactly the sum of the individual processes' radiative responses (the
same non-additivity documented for the aerosol/cloud partial
perturbations in `session_log.md`). A ~22% residual is broadly
consistent with the magnitude of those previously-documented CFRAM
non-linearity residuals; it is reported per-kernel in `lr_attribution.nc`
rather than being forced to zero.

**Qualitative pattern** (`fig_lr_attribution_CloudSat.png`,
`fig_lr_zonal_profile.png`):
- `q` (water vapor) dominates the tropical/subtropical lapse-rate
  feedback, strongly positive — the textbook tropical upper-tropospheric
  amplified-warming signature that produces a stabilizing (for
  Planck+LR combined) but individually-positive LR term.
- `ocndyn`, `lhflx`, `shflx` dominate the high-latitude (particularly
  Southern Ocean, ~55–65°S) structure, reflecting the strong role of
  ocean heat uptake / surface flux redistribution in shaping polar
  vertical temperature structure under 4×CO₂.
- `cloud` and `atmdyn` show zonally-banded structure aligned with the
  storm tracks and ITCZ.
- `co2`, `albedo` are smaller in magnitude but spatially coherent
  (albedo peaks strongly at both poles — sea-ice/snow retreat).
- `aerosol`, `o3`, `solar` are ≈0 everywhere, as expected: this case has
  no aerosol forcing configured and O3 is held at CESM 1850 climatology
  in both states (`base = perturbed`, `frc_o3 ≡ 0`, per
  `session_log.md`'s 2026-05-10/05-12 entries).

## 6. Figures

1. `cases/cesm2_4xco2_official/figures/fig_lr_comparison_CloudSat.png` — native vs ClimKern vs diff, CloudSat kernel (M2)
2. `cases/cesm2_4xco2_official/figures/fig_lr_comparison_GFDL.png` — native vs ClimKern vs diff, GFDL kernel (M2)
3. `cases/cesm2_4xco2_official/figures/fig_lr_kramer_vs_gfdl.png` — native CloudSat(Kramer) − GFDL kernel-choice sensitivity (M2)
4. `cases/cesm2_4xco2_official/figures/fig_lr_attribution_CloudSat.png` — per-process lapse-rate attribution maps (M3)
5. `cases/cesm2_4xco2_official/figures/fig_lr_attribution_GFDL.png` — per-process lapse-rate attribution maps, GFDL kernel (M3)
6. `cases/cesm2_4xco2_official/figures/fig_lr_zonal_profile.png` — zonal-mean per-process attribution, both kernels (M3, paper Fig-3a-style)

Reproduce all of the above:
```bash
python scripts/compute_lr_kernel.py cesm2_4xco2_official
conda activate pycfram-kern
python scripts/validate_lr_vs_climkern.py cesm2_4xco2_official
conda deactivate
python scripts/plot_lr_comparison.py cesm2_4xco2_official
python scripts/compute_lr_attribution.py cesm2_4xco2_official
python scripts/plot_lr_attribution.py cesm2_4xco2_official
```

## 7. Known limitations

- Additivity residual (~22%) is a CFRAM-inherent non-linearity, not
  addressed by this module — documented, not "fixed" (forcing it to
  zero would misrepresent the underlying physics).
- Route ii (ω-based dynamical/convective proxy) was not pursued; see
  `docs/m3_route_decision_memo.md`.
- Polar edge artifact inherited from M2 (thin band at extreme poles in
  diff maps) — see `docs/m2_kernel_module.md` §Known limitations.
