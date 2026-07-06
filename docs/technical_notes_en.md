# pyCFRAM Technical Notes (concise English version)

This is a condensed counterpart of [`technical_notes_zh.md`](technical_notes_zh.md).
The Chinese version is authoritative for deep background and exhaustive
troubleshooting; this file covers the essentials a new user needs to set up,
run, and interpret pyCFRAM.

---

## 1. What pyCFRAM does

pyCFRAM is a Python-orchestrated implementation of the **Climate Feedback–
Response Analysis Method (CFRAM)** that decomposes an observed temperature
change into partial contributions from individual physical processes:

```
ΔT_i = −(∂R/∂T)⁻¹ · F_i          (per process i)
ΔT_total ≈ Σ_i ΔT_i              (linearised)
```

`F_i` is the radiative forcing (W m⁻²) for process *i* obtained by perturbing
one variable at a time through RRTMG. `∂R/∂T` is the Planck response matrix,
also computed by RRTMG via +1 K layer-by-layer perturbations. The matrix
inverse is applied in Python.

**References**: Lu & Cai 2009 / Cai & Lu 2009 (CFRAM formulation);
Zhang et al. 2022 (first CFRAM study to include aerosol, "CFRAM-A" framework);
Wu et al. 2025 (the extreme-heat-event application reproduced here).

---

## 2. Architecture

```
Fortran radiation engine (per grid point) — RRTMG or Fu, picked by case.yaml
  ├── base + perturbed + 8 partial-perturbation rad_driver calls
  ├── Planck matrix (∂R/∂T) via nlayer+1 LW-only perturbations
  ├── 6 per-species aerosol perturbation calls (RRTMG only)
  ├── cloud LW/SW snapshot (zero extra RRTMG invocations)
  └── Fu only: dual MC sub-column overlap patterns (base + warm)
  → writes 17 forcings (9 bulk + 2 cloud_lw/sw + 6 species) + drdt_inv

Python (scripts/run_parallel_python.py)
  ├── dT = −(∂R/∂T)⁻¹ · F for every radiative forcing
  ├── dT for non-radiative lhflx/shflx via the same Planck matrix
  │     (forcing on the surface row only — see §7)
  ├── atmdyn/sfcdyn/ocndyn/dry from frc_full (warm-state ΔR) via the
  │     same Planck matrix:  dry = atmdyn(atm rows) + sfcdyn(surface row),
  │     sfcdyn = ocndyn + lhflx + shflx (exact identity)
  └── multiprocessing.Pool, one worker per grid point
```

**Why no OpenMP inside Fortran?** RRTMG has `SAVE` state in its absorption-
coefficient modules that is not thread-safe. Process-level isolation via
`multiprocessing` avoids the race cleanly. Each worker writes its own
`tmpdir/data_prep/`, invokes the single-column executable, and reads back
the forcings.

---

## 3. Install

Prerequisites:
- Python 3.8+ with `numpy`, `netCDF4`, `matplotlib`, `cartopy`, `scipy`, `pyyaml`
- `gfortran` ≥ 9 (or `ifort`, see `fortran/makefile.ifort`)
- LAPACK / BLAS — **on HPC clusters you must load the module** (e.g.
  `module load lapack`). A silent-exit symptom (absurdly high `pts/s`
  throughput and all-NaN results) almost always means `ldd cfram_rrtmg_1col`
  shows `not found`. See `technical_notes_zh.md` §13.9 for the full diagnosis.
- Accounts: Copernicus CDS (ERA5) + NASA Earthdata (MERRA-2) if you want to
  run the end-to-end pipeline; not needed to read pre-built cases.

Build:

```bash
cd fortran
make                       # RRTMG single-column (default): cfram_rrtmg_1col
make fu                    # Fu single-column: cfram_fu_1col
make TOOLCHAIN=gnu         # Mac / gfortran-only hosts (conda LAPACK)
make TOOLCHAIN=intel       # HPC ifort+MKL (default; e.g. hqlx220)
cd ..
```

Both binaries infer `nlev` at runtime from `data_prep/plev.dat` size, so a
single binary handles any vertical grid (17 / 19 / 30 / 37 / ...) without
recompile.

---

## 4. End-to-end run (EH22 as example)

```bash
python3 scripts/download_era5_flux.py          # first time only, ~5 GB
python3 scripts/download_merra2_aerosol.py     # first time only, ~35 GB
python3 scripts/build_case_input.py --case eh22    # ERA5+MERRA-2 → NetCDF
python3 run_case.py eh22                           # extract → run → plot
```

`run_case.py` chains `build → extract → run → plot`; use `--step <name>`
to run a single stage. A full `--step run` on a 256-core host with 80
workers finishes one EH case in about 10 minutes (~20 pts/s).

---

## 5. Input format

Per case, under `cases/<case>/input/`:

| File | Contents |
|---|---|
| `base_pres.nc` / `perturbed_pres.nc` | (time, lev, lat, lon) 3D fields: `ta_lay, q, o3, camt, cliq, cice, co2, bc, ocphi, ocpho, sulf, ss, dust` |
| `base_surf.nc` / `perturbed_surf.nc` | (time, lat, lon) 2D fields: `ts, ps, solar, albedo` |
| `nonrad_forcing.nc` (optional) | 3D ΔF: `lhflx, shflx` — `sfcdyn/ocndyn/atmdyn` are computed internally |

`lev` is surface → TOA (e.g. 1000, 975, …, 1 hPa). pyCFRAM flips internally
to RRTMG's bottom-to-top order. Full reference: [`input_spec.md`](input_spec.md).

---

## 6. Output variables in `cfram_result.nc`

Shape `(lev, lat, lon)` with surface at `lev[-1]`. Fill value `-999`.

| Group | Variables | Notes |
|---|---|---|
| Radiative (bulk) | `frc_/dT_{co2, q, ts, o3, solar, albedo, cloud, aerosol, warm}` | 9 terms |
| Cloud LW/SW split | `frc_/dT_{cloud_lw, cloud_sw}` | Exactly additive: `cloud = lw + sw` to float64 precision |
| Aerosol species | `frc_/dT_{bc, ocphi, ocpho, sulf, ss, dust}` | Sum ≈ bulk; small residual = non-linear coupling |
| Non-radiative | `dT_{lhflx, shflx}` | From `nonrad_forcing.nc`; forcing placed on the **surface row only** |
| Derived | `dT_{atmdyn, sfcdyn, ocndyn, dry, observed}` | All from `frc_full` (warm-state ΔR) via the Planck inverse. `dry = atmdyn + sfcdyn`; `sfcdyn = ocndyn + lhflx + shflx` (both exact to float precision) |

---

## 7. Key algorithmic choices

**Per-species aerosol via RRTMG (path A), not external pre-computed forcings (path B).** For each species the Fortran engine builds a mixed state
where only that species is swapped base → warm, using AOD-weighted
combination for SSA and *g*:

```
tau_mix      = tau_bulk_base − tau_base_spc(i) + tau_warm_spc(i)
tau·ssa_mix  = (Σ tau·ssa)_base − (tau·ssa)_base_spc(i) + (tau·ssa)_warm_spc(i)
tau·ssa·g_mix = same with g
ssa_mix = tau·ssa_mix / tau_mix
g_mix   = tau·ssa·g_mix / tau·ssa_mix
```

This retains aerosol–aerosol and aerosol–cloud non-linear coupling that a
linear (path B) Planck-inverse split would discard. The price is a 2–5 %
additivity residual between `Σspecies` and the bulk aerosol forcing.

**Cloud LW/SW split — no extra RRTMG calls.** The LW and SW flux
components of the cloud perturbation are snapshotted from the shared
workspace immediately after the cloud `rad_driver` call. The resulting
`dT_cloud_lw + dT_cloud_sw` matches `dT_cloud` to float64 precision
(residual ~10⁻¹⁵).

**Surface `dT` in Python.** Unlike the reference Fortran code, which zeroed
surface `dT` by solving only the atmospheric sub-block, pyCFRAM inverts the
full `(nlayer+1) × (nlayer+1)` Planck matrix in Python and produces a
physically meaningful surface response.

**Dynamics terms from `frc_full`, and the surface-flux identity.** `frc_full`
is the radiative imbalance evaluated in the *warm* atmosphere (`rad_1d_full`,
not a base-state linearisation). Four dynamics terms are derived from it with
the same Planck inverse: `atmdyn` (forcing on atmosphere rows only), `sfcdyn`
(surface row only), `dry = atmdyn + sfcdyn` (full column), and `ocndyn`
(surface row = `−frc_full[sfc] − Δlh − Δsh`). Because `lhflx`/`shflx` forcings
are placed on the **surface row only**, the identity

```
dT_sfcdyn ≡ dT_ocndyn + dT_lhflx + dT_shflx
```

holds to machine precision at every level. (A bug fixed 2026-06-01 had copied
the lhflx/shflx vertical profile onto atmosphere rows, so `dT_lhflx/dT_shflx`
carried a spurious atmospheric Planck response absent from `ocndyn`, breaking
the identity with a structured ~2 K residual at all levels — the discrepancy
reported in `raw/SFCDYN-问题.docx`. The terms are **not** independent CFRAM
processes; `sfcdyn` is by construction the sum.)

**ERA5 6-hourly accumulation ×6.** CDS returns 1-hour increments per
6-hour step for `ssrd/ssr/tisr/slhf/sshf`. `data/era5_source.py` multiplies
by 6 before dividing by 86400. Missing this factor yields a 6× undercount
of `lhflx/shflx` and flips the sign of `dT_cloud`.

**Fu radiation — dual MC sub-column overlap.** When the Fu engine is
selected (`run.executable: cfram_fu_1col`), the radiation engine maintains
two separate Monte Carlo sub-column overlap patterns:

- `base_no_cloud` — sampled (via `ran3`) from `cc_base`, used by cases 0,
  2, 3, 4, 5, 6, 7 + the Planck-matrix `drdt` perturbations.
- `warm_no_cloud` — sampled from `cc_warm`, used by cases 1 (warm),
  8 (cloud), 9 (full).

This mirrors the OLD Fortran CFRAM (`raw/CFRAM.zip`: GW-base.f writes
`base_no_cloud_out`, GW-warm.f writes `warm_no_cloud_out`, GW-cloud.f
reads warm). Using a single MC pattern across both base and warm cloud
states (the pre-2026-05-10 pyCFRAM behaviour) misaligns the realised
sub-column cloud fraction with the `cc_warm` field actually loaded into
the radiation arrays, producing a ≈ +1.8 K mirror flip in CLDL/CLDS at
the global mean. Optional `data_prep/{base,warm}_no_cloud_seed.dat` may
preload OLD CFRAM's exact patterns for bit-perfect single-cell
reproduction.

---

## 8. Known limitations

- **Dust species correlation vs paper** is poor (|corr| < 0.25 on EH22 /
  EH13). Two factors: dust surface `dT` magnitudes are tiny (~0.01–0.1 K),
  and the GOCART single effective-size approximation loses skill in that
  small-signal regime. Future work could use the five DU001–DU005 size
  bins independently.
- **Sulfate bias is systematically high** (~+0.2 K). Single effective
  radius of 0.16 μm for SO₄ optics may not match Wu et al.'s treatment.
- **Cloud / aerosol complementary bias**: path A preserves non-linear
  optical coupling, paper's path B does not. They differ by 10–30 % on
  individual terms but the `cloud + aerosol` sum correlation stays
  above 0.95.
- **OLD CFRAM reference run uses a corrupted O3 input.** The
  collaborator's CESM2 4×CO2 dataset (`raw_data/cesm2_cmip6/CFRAM_xiaominghu/`)
  has `o3_base.dat == o3_warm.dat == hus_base.dat` (md5 `1cefb325...`,
  byte-identical). The OLD radiation routine reads water-vapor values
  into the O3 channel, ~5 orders of magnitude too high. `frc_o3 = 0`
  because base = warm by accident, but the absolute LW/SW radiation
  budget (and therefore drdt and frc_full) is perturbed by the
  mis-typed absorber. This propagates 5–15 % residual bias into
  DYN/ATM/OCH terms. pyCFRAM's `cesm2_4xco2_*` cases use physically
  correct O3, so they should not be expected to bit-match the OLD
  reference at the global field level — single-column dT differences
  ≤ 0.003 K are obtained only when both runs are fed the same broken O3.
- (historical) `core/radiation.py`, `core/decomposition.py`,
  `core/planck_matrix.py`, `core/cfram_runner.py`, `core/fortran_io.py`
  and the `tests/` package were removed as dead code — they were
  residue of an earlier pure-Python-port attempt and never used by the
  current RRTMG-driven pipeline.

---

## 9. Three things that will bite you

Distilled from `technical_notes_zh.md` §13:

1. **HPC LAPACK/BLAS not loaded → silent Fortran exit → NaN results + 1000+
   pts/s throughput.** Check `ldd fortran/cfram_rrtmg_1col`.
2. **ERA5 accumulation ×6 factor.** Diagnose via `build_case_input.py`
   log line `solar domain mean`. Correct value ≈ 420 W m⁻²; broken
   value ≈ 70 W m⁻².
3. **`ndarray.tofile()` is always C-order**, regardless of the array's
   Fortran flags. `np.asfortranarray(arr).tofile(...)` is a no-op for
   byte layout. The Fortran reader's implicit-DO loop order must match
   C-order disk bytes (innermost axis varies fastest on disk).

---

## 10. Single-column climlab validation

A clear-sky radiative-convective equilibrium (RCE) sanity check that runs in
1.5–2 s, independent of any external reference (no `paper_data/`, no OLD
CFRAM). It pairs two climlab columns (1×CO2 = 348 ppm, 4×CO2 = 1392 ppm,
Manabe RH-fixed, ConvectiveAdjustment 6.5 K/km, ΔTs ≈ +4.59 K) and writes
pyCFRAM-format NetCDFs straight into `cases/climlab_4xco2/input/` and
`cases/climlab_4xco2_fu/input/`.

### Run

```bash
# 1. Solve two RCE equilibria → write pyCFRAM input (~30 s for both)
/path/to/conda/python experiments/climlab_validation/run_rce_4xco2.py

# 2. Decompose with each engine (single-cell auto-detected: no mp.Pool spin-up)
python3 run_case.py climlab_4xco2     --step run     # RRTMG
python3 run_case.py climlab_4xco2_fu  --step run     # Fu

# 3. 4-panel vertical-profile diagnostic (Σ closure + co2 + q + dry)
python3 scripts/plot_singlecol_profile.py climlab_4xco2
python3 scripts/plot_singlecol_profile.py climlab_4xco2 climlab_4xco2_fu   # overlay both
```

### Higher-order CFRAM options (RCE only)

Single-column clear-sky cases expose the linearisation residual of 1st-order
CFRAM, which gets hidden by aerosol/cloud/dynamics in real-atmosphere
applications. pyCFRAM offers four opt-in `radiation:` flags that
progressively reduce specific Taylor-expansion terms; all default to the
standard 1st-order behaviour, and **none should be enabled for ERA5/CESM2
runs** — they assume the radiation problem is near-equilibrium with no
non-radiative coupling.

| Flag | Default | Effect |
|---|---|---|
| `drdt_eval: midstate` | `base` | Planck Jacobian at (T_base+T_warm)/2. Cancels the single-variable R_TT term (1-variable 2nd-order Taylor). +1 rad_driver call/cell. |
| `drdt_probe: centered` | `onesided` | T_j ± 0.5K centred FD inside `calc_drdt` instead of +1K one-sided. Cancels R_TT in the Jacobian construction itself. 2× rad_driver_lw calls in calc_drdt. |
| `co2_handling: midstate` | `base` | `frc_co2 = R(T_mid, q_mid, co2_warm) − R(T_mid, q_mid, co2_base)` instead of in the cold base atmosphere. Cancels the ∂²R/∂T∂C cross-term. +2 rad_driver calls/cell. |
| `q_handling: midstate` / `feedback` | `independent` | Diagnostic only. `midstate` makes upper-trop closure 3× worse on Manabe RH-fixed RCE — the base-path R_Tq cross-term aligns with the constrained-RH solution by accident, so cancellation removes a "correct" contribution. `feedback` flips the sign of R_Tq (one-sided warm path). |

### Verified closure

`dT_obs − (dT_co2 + dT_q + dT_ts + dT_dry)` RMS by altitude band, on
`climlab_4xco2` (RRTMG):

| Configuration | mid-trop (450–650 hPa) | upper-trop (200–400 hPa) | strato (<200 hPa) |
|---|:---:|:---:|:---:|
| Baseline 1st-order | 0.45 K | 0.07 K | 0.30 K |
| + `drdt_eval: midstate` | 0.45 K | 0.07 K | **0.28 K** (stratosphere fixed) |
| + `co2_handling: midstate` | **0.24 K (-47 %)** | 0.04 K | 0.28 K |
| + `drdt_probe: centered` | **0.21 K (-55 %)** | **0.02 K (-68 %)** | **0.25 K** |

The remaining ~0.2 K mid-trop residual is the physical floor of 1st-order
CFRAM on this setup: climlab convective-adjustment mass redistribution
(non-radiative, fundamentally outside the J⁻¹·ΔR framework), R_TTT
high-order in T, and Clausius–Clapeyron quadratic in Δq. Closing further
requires a 2nd-order CFRAM with an explicit Hessian — structurally a
different algorithm. See `temp.md` / the experiment series notes for the
full Taylor-expansion derivation and per-altitude diagnostics.

---

## 11. Remote operations (HKUST-specific)

Production runs happen on the `hqlx*` cluster:

- Edit locally → `rsync` Mac → `mini` (macOS jump box) → `mini` → `hqlx*`.
  **Always start the sync from the Mac**; a direct `ssh mini 'rsync ...'`
  uses mini's possibly-stale copy.
- Default production host: **`hqlx220`** (384 cores, ifort 2022.2 + MKL
  LP64 + netCDF4). Source `/home/lzhenn/.bashrc_liquor_i22wrf415` before
  `make` / `python3 run_case.py`. Backup hosts with the same ifort+MKL
  toolchain: `hqlx204` / `hqlx205` (256 cores) and `hqlx221–223` (384
  cores) — the ifort binary built on `hqlx220` runs on all of them via the
  NFS-shared `/home/lzhenn/work`.
- **`hqlx74` cannot run the Python runner**: its glibc is 2.17, but the
  conda `netCDF4` wheel needs `GLIBC_2.25` and fails to import. Use it for
  pure-Fortran or download tasks only.
- A single Fortran binary handles any vertical grid at runtime — `nlev` is
  inferred from `data_prep/plev.dat` (8 bytes per level) on launch. No
  per-grid rebuild needed.

Checked-in compiled binaries are **not** committed; a fresh clone on any
system runs `cd fortran && make` (intel by default, `TOOLCHAIN=gnu` on Mac
or any gfortran-only host).

---

## 12. Phase 3: DAMIP multi-model data source

A registered data source (`data/cmip6_damip_source.py`, `@register_source
('cmip6_damip')`) lets pyCFRAM ingest any CMIP6 DAMIP single-forcing
experiment (`hist-aer`, `hist-GHG`, ...) from any contributing model,
without touching `core/` or `fortran/` — same generic writer
(`scripts/build_case_input.py`) and same Fortran engine as ERA5/CESM2. Full
architecture, IO contract, and worked bug-fix examples are in
[`docs/m4_damip_module.md`](m4_damip_module.md); per-model support matrix
and known-issues are in
[`docs/m5_multimodel_userguide.md`](m5_multimodel_userguide.md). Short
version, since the Chinese technical notes carry the exhaustive account:

- **Model-agnostic numerics** live in `data/cmip6_common.py`: calendar-aware
  time decoding (`cftime`, correct for noleap/365_day/360_day/gregorian/
  proleptic_gregorian — the naive `days/365.0` arithmetic used before this
  module is silently wrong for anything but noleap), hybrid-sigma-pressure
  coefficient-naming detection (`p = a·p0 + b·ps` vs `p = ap + b·ps`),
  log-pressure plev re-interpolation, bilinear horizontal regrid (needed
  when a variable is published on a coarser native grid than its own
  model's siblings — confirmed for MRI-ESM2-0's `o3`), and an analytic TOA
  insolation fallback for the ~40% of candidate models that don't publish
  `rsdt`.
- **Per-model quirks** are pure YAML (`configs/damip_models.d/<model>.yaml`)
  — default variant/grid, calendar, and an explicit hybrid-coefficient
  override for models whose `formula_terms` CF attribute is missing (IPSL)
  or non-standard (CNRM's `formula_term`, singular). Adding a model that
  fits the standard cases needs **zero Python changes** — see
  `docs/m5_multimodel_userguide.md`'s NorESM2-LM worked example.
- **Missing cloud/O₃/`rsdt` is the norm, not the exception**: of 13
  candidate hist-aer models surveyed, only 3 publish all needed variables.
  A documented decision tree (build-time, not a post-hoc hook — the
  writer's `validate_states()` rejects non-finite values before any file
  is written) resolves each gap to a physically inert fallback (e.g.
  injected O₃ climatology with `frc_o3 ≡ 0` by construction) and records
  the choice in `provenance.json`, surfaced in a human-readable
  `*.summary.txt` after every run.
- **ESGF acquisition** (`data/esgf_fetch.py` + `scripts/download_damip.py`)
  is pure stdlib `urllib` — no `requests`/`intake-esgf` — and downloads
  without authentication (DKRZ/CEDA Solr search + HTTP Range requests).
- Coverage of the three new modules: 96–99% (`coverage run -m pytest
  tests/ -q && coverage report -m --include="*/data/cmip6_common.py,
  */data/cmip6_damip_source.py,*/data/esgf_fetch.py"`), well above the
  60% contract gate.
