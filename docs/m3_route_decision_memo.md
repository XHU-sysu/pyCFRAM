# M3 Route Decision Memo

**Date**: 2026-07-02
**Decision**: Route i — attribute the lapse-rate radiative perturbation to
pyCFRAM's existing physical-process temperature decomposition
(`dT_q`, `dT_atmdyn`, `dT_sfcdyn`, `dT_cloud`, `dT_ocndyn`, `dT_co2`,
`dT_o3`, `dT_albedo`, `dT_solar`, `dT_aerosol`, `dT_lhflx`, `dT_shflx`).

## Options considered (docs/plan.md §7.0)

- **Route i**: treat each CFRAM partial temperature response `dT_X(p)`
  as its own mini "temperature response", extract its non-uniform
  component `dT_X(p) − dT_X(skin)`, and push it through the exact same
  kernel machinery built for M2 (`core/lr_kernel.py::extract_and_apply`).
  Zero new physics, zero new data — pure reuse of quantities already
  computed and already validated (M2 passed at corr ≥ 0.9997).
- **Route ii**: bring in an omega (ω, vertical pressure velocity)-based
  convective/dynamical proxy to split the non-uniform warming into
  convective vs large-scale-dynamical contributions, independent of
  CFRAM's own process decomposition.

## Why route i

1. **Data availability**: `grep`/`find` across both the local and
   hqlx210 `raw_data/` trees for `wap`/`omega` (the CMIP6 vertical
   velocity variable) returns nothing — the CESM2 4×CO2 CMIP6 build
   pipeline (`scripts/build_cesm2_official.py`) never ingested it. Route
   ii would require a new data-acquisition + build step before any
   analysis code could even start, which is out of scope for what's
   already a working, validated M2 module.
2. **Occam's razor** (docs/plan.md's own framing, "奥卡姆最省"): route i
   answers the contract's actual question — "which physical process
   shaped the lapse-rate structure" — using quantities pyCFRAM *already
   produces and has already closure-validated* (the `dT_sfcdyn =
   dT_ocndyn + dT_lhflx + dT_shflx` identity residual is at machine
   precision post-fedcf4d, see WP-M2.R). No new assumptions, no new
   physics, no new failure modes.
3. **Consistency with M2**: route i reuses the exact same
   `core/lr_kernel.py::extract_and_apply` function that was
   just cross-validated against ClimKern to corr ≥ 0.9997 — so any
   attribution result inherits that same validated numerical machinery,
   rather than introducing an independently-validated (or
   unvalidated) omega-diagnostic pathway.
4. **The plan's own default is route i**, with route ii explicitly
   marked "可选 / stretch" and gated on client sign-off that isn't
   available in this execution context.

## What was verified before deciding

- `core/lr_attribution.py` + `scripts/compute_lr_attribution.py`
  (route i's implementation) were run end-to-end against the real,
  freshly-rerun `cesm2_4xco2_official` case:
  ```
  CloudSat   additivity |residual| mean = 1.40 W/m^2 (21.8% of |total| mean)
  GFDL       additivity |residual| mean = 1.75 W/m^2 (23.1% of |total| mean)
  ```
  A ~22% additivity residual (Σ per-process ΔR_LR ≠ total ΔR_LR from
  `dT_observed`) is expected and documented (docs/plan.md R7): CFRAM's
  first-order Taylor decomposition is linear in each `dT_X`, but the
  underlying radiative transfer is not perfectly linear across
  processes — the same character of residual already documented for
  the aerosol/cloud additivity checks in `session_log.md`. This is not
  a bug in the kernel-attribution module; it is the CFRAM
  decomposition's own inherent non-linearity, made newly visible by
  putting a kernel on top of it.

## Outcome

Proceeding with route i as the sole M3 deliverable. Route ii is
recorded here as a known, well-scoped future extension (would need: (a)
`wap` CMIP6 variable ingestion in `build_cesm2_official.py`, (b) a
convective/dynamical partitioning threshold on ω, (c) its own
cross-validation against an independent proxy) — out of scope for this
milestone per the fallback in docs/plan.md §7.0 ("若甲方要 ii 但 ω 数据
不可得 → 退回 i").
