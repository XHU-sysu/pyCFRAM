# Changelog

All notable changes to pyCFRAM are recorded here.
The format loosely follows [Keep a Changelog](https://keepachangelog.com/)
and this project uses date-based milestones (no semver tags yet).

## [Unreleased]

### Added
- BSD 3-Clause `LICENSE` with explicit third-party notices for upstream
  RRTMG and the GOCART optical-property tables.
- GitHub Actions CI (`.github/workflows/ci.yml`) for syntax compilation
  and Ruff lint on every push / PR to `main`.
- `CHANGELOG.md` and the `technical_notes_en.md` concise English
  counterpart of `technical_notes_zh.md`.

---

## 2026-06-01 — sfcdyn forcing fix + Fu/RRTMG surface-row unbundle + India case

### Fixed
- **`sfcdyn = ocndyn + lhflx + shflx` identity restored** (the headline
  fix). The nonrad loader in `run_parallel_python.py` was copying the full
  `lhflx`/`shflx` vertical profile onto the atmosphere rows of the forcing
  vector (`frc_full[:NLEV]`), while `dT_ocndyn` subtracts only the surface
  flux. As a result `dT_lhflx`/`dT_shflx` carried a spurious atmospheric
  Planck response that `ocndyn` never removed, breaking the design identity
  with a structured residual at every level (global mean ~2.15 K, max
  ~190 K on `cesm2_4xco2_official`) — the discrepancy reported in
  `raw/SFCDYN-问题.docx`. Now the flux is placed on the **surface row only**;
  the identity holds to machine precision (max |residual| 3.6×10⁻¹⁴ K).
  Diagnostic: `dT_sfcdyn/dT_ocndyn` is column-constant (both clean
  surface-forcing responses) while `dT_lhflx/dT_sfcdyn` was not — proving
  the leak was confined to the two flux terms.

### Changed
- **Surface-row Planck perturbation unbundled** in both engines. RRTMG
  `drdt.f90`: the surface row perturbs `ts` only (reverted an earlier
  attempt that also bumped `t_atm[nlayer]` — RRTMG's `tint_cfram(nlayer)=ts`
  is independent of the bottom layer T). Fu `cfram_fu_1col.f90`:
  `pt(nv1) = pt(nv)` so the skin-T perturbation no longer pins the bottom
  atmosphere temperature. Skin T and bottom-layer T are now physically
  decoupled in the Jacobian for both engines.
- **Plot scripts tag the radiation engine** in output filenames:
  `_pycfram-{rrtmg,fu}` for pyCFRAM runs, `_oldcfram` for OLD-CFRAM
  comparison plots (`plot_13panel_global.py`, `plot_13panel_polar.py`,
  `plot_closure_profile.py`).
- `plot_13panel_polar.py` switched to `pcolormesh` to dodge a
  shapely-2.x / cartopy-0.22 `MultiPolygon` contourf rendering bug that
  produced white artifacts over the pole.
- Documentation (`README`, `algorithm_spec` §8.1, `input_spec`,
  `technical_notes_{en,zh}`): the dynamics terms `atmdyn/sfcdyn/ocndyn/dry`
  are forcing-based responses to `frc_full` (warm-state ΔR), **not**
  residuals; documented `dry = atmdyn + sfcdyn` and the surface-flux
  identity above.

### Added
- **OLD-CFRAM-style closure plotting.** `plot_closure_profile.py` (and
  `plot_closure_profile_oldcfram.py`) drop `ts` from the closure sum and
  define the plotted `atmdyn` as the closure residual `dT_obs − Σrad −
  Σnonrad − ocndyn`; panels with `max|prof| < 0.1 K` are filtered.
- **Contract M1 case ② India–Bangladesh wet-heat 2023** —
  `cases/india_wb23/case.yaml` (warm core Apr 17–20, CO2 385→421.3 ppmv).
  ERA5 6-hourly April download pipeline `scripts/download_era5_india23.py`
  (CDS, 2-var-pair requests under the per-request cost cap) and
  `download_era5_india23_ncar.py` (NCAR RDA mirror). See `docs/cases.md`.
- `fortran/fu_helpers.f` split out of the Fu source tree.
- Plot helpers `plot_13panel_oldcfram.py`, `plot_13panel_polar_from_slice.py`
  (surface-slice fallback when the full NC trips the cartopy polar bug).

### Notes
- `hqlx74` confirmed unusable for the Python runner (glibc 2.17 < the
  `GLIBC_2.25` the conda `netCDF4` wheel needs). Production re-runs use
  `hqlx220` (build) + `hqlx205`/`204`/`221–223` (ifort binary via NFS).

---

## 2026-05-10 — `radiation:` schema + climlab single-column validation

### Added
- **`radiation:` block in `case.yaml`** — declarative radiation-engine selection
  (`scheme: fu | rrtmg`, extensible via `core/config.RADIATION_SCHEMES`).
  Replaces the engine-specific `run.executable` field (legacy still honoured).
- **`radiation.output_terms` filter** — list of `dT/frc` rows to write to
  `cfram_result.nc`. Idealized cases (clear-sky RCE) use this to drop the
  always-zero cloud / aerosol / o3 / albedo rows. Derived dynamics terms
  (`dry / atmdyn / sfcdyn / ocndyn / lhflx / shflx`) are *not* filtered —
  they are physically required for closure.
- **Auto-cap `nproc` by grid size** — single-cell cases (climlab RCE, etc.)
  short-circuit `multiprocessing.Pool` and call `process_column` synchronously,
  saving ~0.5-1 s of pool spin-up.
- **`get_plev()` auto-derives from input NetCDF** — third resolution path
  alongside `case.yaml grid.pressure_levels` and `defaults.yaml`. A climlab
  30-level RCE column or any custom grid Just Works without yaml duplication.
- **CliMLab 4×CO2 single-column validation case** —
  `experiments/climlab_validation/run_rce_4xco2.py` runs clear-sky RCE at
  348 ppm and 1392 ppm CO2 (climlab default 30-sigma grid) and writes
  pyCFRAM standard 1×1 input directly. New cases `cases/climlab_4xco2/`
  (RRTMG) and `cases/climlab_4xco2_fu/` (Fu) decompose the ΔTs ≈ +4.6 K
  ECS-equivalent response. Verified closure
  `dT_co2 + dT_q + dT_dry ≈ dT_observed` within 1.2 % residual on both
  engines — a clean independent sanity check of pyCFRAM's CFRAM math.

### Changed
- All existing case.yaml migrated from `run.executable` to
  `radiation.scheme`: eh13, eh22, cesm2_4xco2{,_fu,_official,_official_fu,
  _official_17p_fu}, climlab_solar{,_fu}.

---

## 2026-05-10 — Fu radiation dual-MC fix + CESM2 4×CO2 workflow

### Added
- **Fu radiation engine** (`cfram_fu_1col`, `cfram_fu_1col.f90`,
  `fortran/Fu/cas_fu_radiation.f`, `fu_helpers.f`, `para.file`) as a
  drop-in alternative to the RRTMG engine. Selected per-case via
  `case.yaml` `run.executable`. Single Fu binary infers `nlev` at
  runtime from `data_prep/plev.dat` size, so 17/19/37-plev cases share
  one binary.
- **Dual MC sub-column overlap patterns.** Apple-to-apple with the
  OLD Fortran CFRAM (`raw/CFRAM.zip` GW-base.f / GW-warm.f / GW-cloud.f):
  the radiation calls now distinguish a `base_no_cloud` pattern (drawn
  from `cc_base`) from a `warm_no_cloud` pattern (drawn from `cc_warm`).
  Cases 0, 2, 3, 4, 5, 6, 7 + drdt use the base pattern; cases 1
  (warm), 8 (cloud), 9 (full) use the warm pattern. Fixes a +1.8 K
  global-mean mirror flip in CLDL/CLDS that pre-fix pyCFRAM showed vs
  the OLD reference (root cause: the warm-cloud RT call was using a
  base-cc-consistent sub-column realization, mis-aligned with the
  cc_warm field actually loaded into the radiation arrays).
- **Intel toolchain support** in the makefile (`make TOOLCHAIN=intel`,
  default on hqlx220/204; `gnu` for Mac local). LAPACK linked LP64 to
  match `-i8` integer kind.
- **CESM2 4×CO2 workflow.** New cases
  `cesm2_4xco2_official_{,fu,17p_fu}` plus helper scripts:
  `scripts/build_cesm2_official.py` (CMIP6 hybrid→plev),
  `scripts/subset_to_17p.py` (19→17 plev to match OLD CFRAM exact grid),
  `scripts/inject_cesm_o3.py` (CESM 1850 climatology O3 — Phase A: base
  = warm), `scripts/mask_subsurface_layers.py` (mask layers below
  ps_warm), `scripts/expand_nonrad_to_column.py` (per-layer
  distribution of column-integrated fluxes).
- **13-panel plotting.** `scripts/plot_13panel_global.py` (global
  surface dT) and `scripts/plot_13panel_polar.py` (north + south polar)
  matching the OLD CFRAM `north.jpg` reference layout.
- `data/cesm2_cmip6_source.py` (CMIP6 raw NetCDF reader +
  hybrid-sigma → plev mass-conserving re-projection).
- `core/heating_profile.py` (column-vertical distribution of
  surface-bulk lhflx/shflx for cases without per-layer tendencies).
- `experiments/climlab_validation/` (CliMLab idealized RCE
  cross-validation with pyCFRAM Fortran path).

### Fixed
- **Fu radiation MC sub-column overlap mismatch** (the major fix this
  release). Single-column DP validation at (159, 144) on the
  collaborator's CESM2 dataset now reproduces all 13 dT terms
  (al/wv/cld/cld_lw/cld_sw/co2/o3/solar/dyn/atm_dyn/sfc_dyn/sh/lh)
  within ≤ 0.003 K of the reference `partial_T_1.grd`. Global re-run
  on `cesm2_4xco2_official_17p_fu` (200 cores, 3 min) gives
  `dT_cloud_lw` mean −0.15 K (was +1.37 K pre-fix → 1.8 K mirror).
- `core/config.py:get_plev` now accepts a `case_cfg` and supports
  per-case `grid.pressure_levels` override (used by 17-plev variants).
- `fortran/cfram_rrtmg.f90`: minor cleanup (legacy full-field driver
  retained but no longer the primary entry point).

### Investigated and documented (not a code change)
- **Collaborator's `partial_T_1.grd` benchmark uses a corrupted O3
  input.** md5 audit found `o3_base.dat == o3_warm.dat == hus_base.dat`
  (md5 `1cefb325...`, byte-identical), i.e. the OLD CFRAM reference
  read water-vapor values into the O3 channel, ~5 orders of magnitude
  too high. `frc_o3 = 0` because base = warm by accident, but the
  absolute radiation budget (and therefore drdt and frc_full) is
  perturbed by the mis-typed absorber. This propagates 5–15 % residual
  bias into DYN/ATM/OCH panels — explaining residual visual differences
  vs `north.jpg` after the dual-MC fix. pyCFRAM's `cesm2_4xco2_*` cases
  use physically correct O3, so they should not be expected to
  bit-match the OLD reference at the global field level. See
  `session_log.md` 2026-05-10 entry.

### Removed
- `scripts/extract_full_field.py` (functionality merged into
  `run_parallel_python.py` per-cell extraction).

### Removed
- `core/aerosol_optics.py` — dead module with no imports anywhere in the
  codebase.
- `core/radiation.py`, `core/decomposition.py`, `core/planck_matrix.py`,
  `core/cfram_runner.py`, `core/fortran_io.py` — residue of an earlier
  pure-Python-port attempt. These modules imported constants
  (`CLOUD_RE_ICE`, `CH4_PPMV`, `ICLD_ON`, `ERA5_PLEVELS`, …) that no
  longer exist in `core/constants.py` and were not used by the current
  RRTMG-driven pipeline.
- `tests/` directory (only test was `test_radiation_prep.py`, which
  imported the stale `core/radiation.py` and could not run). With the
  whole package gone, the CI no longer needs `tests/` in its
  `compileall` list.
- Placeholder / gitignored-only directories `diagnostics/`, `tasks/`,
  `validation/`.

---

## 2026-04-18 — Cloud LW/SW forcing split and aerosol-species decomposition

### Added
- **Per-species aerosol perturbation in RRTMG.** The Fortran engine now
  loops over six aerosol species (`bc, ocphi, ocpho, sulf, ss, dust`),
  constructing for each species a mixed optical state where only that
  species is swapped from the base to the perturbed value (AOD-weighted
  combination for SSA and *g*). Outputs `frc_<species>.dat`; the Python
  runner applies the Planck inverse to produce `dT_<species>` in
  `cfram_result.nc`. Additivity residual `bulk − Σspecies` is 2–5%,
  interpreted as aerosol–aerosol / aerosol–cloud non-linear coupling.
- **Cloud LW/SW split forcing.** The LW and SW radiative-flux components
  of the cloud perturbation are snapshotted immediately after the cloud
  `rad_driver` call (zero extra RRTMG invocations). New variables
  `frc_cloud_lw`, `frc_cloud_sw`, `dT_cloud_lw`, `dT_cloud_sw` are
  exactly additive to the bulk cloud terms (float64 precision, residual
  ≈ 1.5 × 10⁻¹⁵).
- **Eight species-dimensioned aerosol-optics binary files** (`_spc.dat`)
  generated by `extract_full_field.py` / `run_parallel_python.py` for
  consumption by the Fortran per-species loop.
- `scripts/verify_phase{1,4,5}_spc.py` and `scripts/verify_cloud_split.py`
  verification harnesses.
- BSD/GPL-compatible documentation updates: README and algorithm_spec
  now reference Zhang et al. (2022) — first CFRAM study to include an
  explicit aerosol term.
- `docs/technical_notes_zh.md` §13.9: HPC LAPACK/BLAS troubleshooting
  (silent Fortran exit → NaN results). Diagnosed via `ldd`.

### Changed
- `docs/algorithm_spec.md`: §8 corrected (surface dT is now solved via
  the full (`nlayer+1`)² Planck inverse in Python — no longer set to
  0.0). §9 rewritten from the current Python aerosol-optics
  implementation (previous section still referenced the retired Matlab
  scripts). §10 coverage table extended with `cloud_lw/sw` and the
  6-species aerosol rows. §11 NCL preprocessing section retired.
- `docs/input_spec.md`: removed `sfcdyn` from required non-radiative
  forcings (now derived internally from energy balance); documented
  the optional path-B per-species aerosol forcing as a legacy
  validation-only input.

### Verified
- EH22 & EH13 end-to-end: per-species spatial correlation vs Wu et al.
  2025 paper data is 0.73–0.93 for BC/OC/Sulf/Seas; dust correlation
  is poor due to small magnitude and the GOCART single-bin
  approximation (known limitation).

---

## 2026-04-12 — Independent ERA5 + MERRA-2 pipeline (v3)

### Added
- `scripts/build_case_input.py` driven by `case.yaml` `source:` block:
  builds `base_/perturbed_{pres,surf}.nc` and `nonrad_forcing.nc`
  directly from ERA5 daily PL/SL and MERRA-2 `M2I3NVAER` — no
  dependence on the Wu et al. paper_data.
- `data/era5_source.py` (`ERA5DailySource`) and `data/merra2_aerosol.py`
  (6-species GOCART-aligned aerosol extraction with vertical log-p and
  horizontal bilinear interpolation).

### Fixed
- **ERA5 6-hourly accumulation ×6 correction** (root cause of an
  earlier cloud-sign flip and ~6× undercount of LH/SH forcings). CDS
  returns 1-hour increments per 6-hour step; the accumulator now
  multiplies by 6. After the fix, cloud `dT` is correctly positive and
  LH/SH match Wu et al. within <5%.
- CO₂ values per case aligned with Wu et al. climatology (EH22: 394.23 /
  416.89 ppmv; EH13: 395.54 / 396.14 ppmv).
- Warm-day indexing: case.yaml now supports an explicit `warm_days`
  list; auto-detection fallback uses a 90th-percentile `mx2t`
  threshold with ≥3 consecutive days.

### Verified
- EH22 WV corr = 0.998, LH/SH corr = 0.993.
- EH13 WV corr = 0.997, LH/SH corr = 0.992.
- Cloud/aerosol show complementary bias (path A vs path B methodology
  difference), total `dT` correlation 0.981 on both cases.

---

## 2026-04-04 — Surface processes, atmospheric dynamics, Fig.3

### Added
- Python decomposition of non-radiative forcings (`lhflx`, `shflx`,
  `sfcdyn`) via the same Planck inverse used for radiative terms.
- Residual-method atmospheric dynamics: `dT_atmdyn = dT_obs − Σ(all
  other decomposed terms)`.
- `scripts/plot_fig3_independent.py` 2-column (paper vs independent)
  spatial maps.

### Verified
- `dT_total` correlation vs paper: EH13 / EH22 = 0.981.
- Surface Process (LH+SH+sfcdyn) correlation: 0.999 / 0.998.

---

## 2026-04-04 — Full-field parallel runner

### Added
- `scripts/extract_full_field.py`: full-field NetCDF → Fortran binary
  with vectorised GOCART aerosol optics.
- `scripts/run_parallel_python.py`: `multiprocessing.Pool`-based
  per-grid-point parallel runner. Each worker writes a unique
  `tmpdir/data_prep/` and invokes `fortran/cfram_rrtmg_1col`. ~20 min
  per EH case at 80 workers on a 256-core host.
- Python-side Planck inverse read (`drdt_inv.dat`) + `dT = −A⁻¹·F`
  solve (moved out of Fortran, matching the new dT-in-Python design).

### Investigated and rejected
- **Fortran OpenMP parallelism** — RRTMG module state is not thread
  safe (`SAVE` variables in the absorption-coefficient lookup). Process
  isolation via `multiprocessing` is the only reliable path.

---

## 2026-04-04 — Single-column validation (initial CFRAM reproduction)

### Added
- `fortran/cfram_rrtmg.f90` modified to output the Planck inverse
  (`drdt_inv.dat`) and use a full `(nlayer+1)²` solve for surface
  `dT`.

### Fixed
- Fortran direct-access `recl` unit mismatch (`*2 → *8` for gfortran
  vs ifort).
- `-999` fill boundary: `nlayer:nlev → nlayer+1:nlev` (previously
  overwrote the last valid layer).

### Verified on single column (115 °E, 32 °N, EH22)
- Surface `dT` diff vs paper: 0.009 K.
- 925 hPa `dT` diff: 0.15 K.
