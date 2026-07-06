# pyCFRAM

A Python-encapsulated implementation of the **Climate Feedback–Response Analysis Method (CFRAM)**. pyCFRAM aims to build a user-friendly, high-efficient, and extendable interface for CFRAM users. It uses `multiprocessing` (and future `f2py`/`mpi4py` support) to exploit the maximum potential of code parallelism without dramatically modifying the original Fortran radiation code.

**Two interchangeable radiation engines** are supported (selectable per-case via `case.yaml`):

- **RRTMG** (default) — modern broadband LW/SW with full GOCART aerosol coupling and per-species partial perturbations (`bc, ocphi, ocpho, sulf, ss, dust`).
- **Fu** — the classical scheme used by the original Fortran CFRAM (Lu/Cai). Bit-exact with the OLD CFRAM reference run when fed the same inputs (validated 2026-05-10, all dT_X ≤ 0.001 K residual at single column).

Reproduces:
- Temperature decomposition from Wu et al. (2025, *J. Climate*) for extreme heat events EH13 (Aug 2013) and EH22 (Aug 2022) over East China (RRTMG path).
- 13-panel global / polar dT decomposition for a CMIP6 CESM2 4×CO2 experiment (Fu path, apple-to-apple with the OLD CFRAM Fortran benchmark).

For the CFRAM methodology, please refer to:
- Lu, J., and M. Cai, 2009: A new framework for isolating individual feedback processes in coupled general circulation climate models. *Climate Dynamics*.
- Cai, M., and J. Lu, 2009: A new framework for isolating individual feedback processes in coupled general circulation climate models. Part II. *Climate Dynamics*.

For a related approach (CFRAM-A) using RRTMG with aerosol independently developed by Zhang et al. (2022):
- Zhang, T., Deng, Y., Chen, J., Yang, S., Gao, P. and Zhang, H., 2022. Disentangling physical and dynamical drivers of the 2016/17 record-breaking warm winter in China. *Environmental Research Letters*, 17(7), 074024.
- Wu, Q., Li, Q., Zhang, T., Sun, X., Yang, S. and Hu, X., 2025. Quantitative attribution of 2013 and 2022 extremely hot events in China: insights from a climate feedback–response perspective. *Journal of Climate*, 38(17), pp.4331–4349.

## Architecture
```
Fortran radiation (per grid point) — RRTMG OR Fu, picked by case.yaml
  ├── base + warm + 8 partial-perturbation rad_driver calls
  ├── Planck matrix (∂R/∂T) via nlayer+1 LW-only perturbations
  ├── 6 per-species aerosol perturbation calls (RRTMG only: bc/ocphi/ocpho/sulf/ss/dust)
  ├── cloud LW/SW component snapshot (zero extra calls)
  └── Fu only: dual MC sub-column overlap patterns
        ├── base_no_cloud (cc_base-consistent) used by base/co2/wv/o3/solar/albedo/ts + drdt
        └── warm_no_cloud (cc_warm-consistent) used by warm/cloud/full
  → writes 9 bulk + 6 per-species + 2 cloud_lw/sw forcings + Planck inverse

Python (decomposition + analysis)
  ├── dT_i = -(∂R/∂T)⁻¹ × frc_i           # radiative terms (all 17, incl. splits)
  ├── dT_lhflx/shflx                       # non-radiative; forcing on surface row only
  ├── dT_atmdyn/sfcdyn/ocndyn/dry          # from frc_full (warm-state ΔR), same Planck matrix
  │     dry = atmdyn + sfcdyn              #   (atm rows / surface row / full column)
  │     sfcdyn = ocndyn + lhflx + shflx    #   exact identity (machine precision)
  └── multiprocessing parallel execution   # embarrassingly parallel over grid points
```

### Selecting the radiation engine

Each case.yaml has a `radiation:` block:

```yaml
radiation:
  scheme: fu                     # fu | rrtmg (extensible — see core/config.RADIATION_SCHEMES)
  output_terms: [co2, q]         # optional: write only these dT/frc rows to NetCDF
                                 # (omit for "everything"; useful for clear-sky RCE
                                 # where cloud / aerosol / o3 rows are always zero)

  # Optional higher-order CFRAM options (all RRTMG-only, default = standard 1st-order
  # CFRAM; only useful for clear-sky single-column validation where the linearisation
  # residual is the dominant error source — for real-atmosphere cases keep defaults).
  drdt_eval:   midstate          # Planck Jacobian at (T_base+T_warm)/2; cancels R_TT
                                 # (single-variable 2nd-order). Resolves the
                                 # stratosphere residual.
  drdt_probe:  centered          # ±0.5K probe in calc_drdt instead of +1K one-sided;
                                 # cancels the R_TT FD term in the Jacobian itself.
  co2_handling: midstate         # frc_co2 evaluated in midstate atmosphere; cancels
                                 # the ∂²R/∂T∂C cross-term (dominant mid-trop source).
  # q_handling: midstate / feedback   # Diagnostic only. midstate makes upper-trop
                                 # worse on RCE; feedback flips R_Tq sign vs base.

run:
  nproc: auto                    # auto-caps to grid size (1×1 → sync, no Pool spin-up)
```

The legacy escape hatch `run.executable: <binary_name>` is still honoured. Both `cfram_rrtmg_1col` and `cfram_fu_1col` infer `nlev` at runtime from the size of `data_prep/plev.dat`, so a single binary per engine handles any vertical grid (17/19/30/37/...) without recompilation. `get_plev()` resolves levels by:

1. `case.yaml` `grid.pressure_levels` (explicit override)
2. Input NetCDF `lev` variable (auto-derive — works for climlab 30-level, CMIP6 19-level, ...)
3. `configs/defaults.yaml` (last-resort fallback)

## Prerequisites

- Python 3.8+ with: `numpy`, `netCDF4`, `matplotlib`, `cartopy`, `scipy`
- `gfortran` (for compiling RRTMG)
- LAPACK/BLAS libraries (on HPC clusters, may require `module load lapack`; see [`docs/technical_notes_en.md`](docs/technical_notes_en.md) §9 or the expanded Chinese `technical_notes_zh.md` §13.9)

## Setup

### 1. Clone and install dependencies

```bash
git clone git@github.com:lzhenn/pyCFRAM.git
cd pyCFRAM
pip install -r requirements.txt
```

### 2. Compile Fortran radiation engine

```bash
cd fortran
make            # default RRTMG: builds cfram_rrtmg_1col (single-column, 37-plev)
make fu         # Fu: builds cfram_fu_1col (single-column, runtime nlev)
make TOOLCHAIN=gnu      # Mac local: gfortran + conda LAPACK
make TOOLCHAIN=intel    # HPC: ifort + MKL (default on hqlx220/204)
cd ..
```

Single-column executables are used by the parallel Python runner. The RRTMG lookup tables (`rrtmg_lw.nc`, `rrtmg_sw.nc`) and Fu source (`fortran/Fu/cas_fu_radiation.f`, `fu_helpers.f`, `para.file`) are included.

Both toolchains tested:
- `intel` (ifort 2021.7 + MKL 2022.2 LP64): production HPC build.
- `gnu` (gfortran from conda + openblas via `liblapack`/`libblas`): Mac local debugging.

### 3. Obtain source data (ERA5 + MERRA-2)

The standard pipeline drives CFRAM directly from ERA5 (state) and MERRA-2 (aerosol). Place these under `pyCFRAM/era5_data/`:

```
era5_data/
├── daily/                                 # ERA5 daily, 2003-2022 Aug
│   ├── era5_pl_{var}_{YYYY}08.nc          # PL: t, q, o3, cc, ciwc, clwc (6-hourly, 37 lev)
│   └── era5_sl_{YYYY}08/                  # SL: subdirectory per month
│       ├── *stepType-instant.nc           # skt, sp
│       ├── *stepType-accum.nc             # ssrd, ssr, tisr, slhf, sshf
│       └── *stepType-max.nc               # mx2t
└── merra2/
    └── M2I3NVAER_{YYYY}08/*.nc4           # 3-hourly, 72 model levels, 13 aerosol species
```

Helper downloaders (CDS / NASA Earthdata credentials required):

```bash
python3 scripts/download_era5_flux.py        # ERA5 PL+SL via CDS API
python3 scripts/download_merra2_aerosol.py   # MERRA-2 M2I3NVAER via GES DISC
```

Optional: to validate against Wu et al. results, also download `paper_data/` (not required for standard runs).

## Quickstart: Reproduce Fig.3

### Step 1: Build case input from ERA5 + MERRA-2

```bash
python3 scripts/build_case_input.py --case eh13
python3 scripts/build_case_input.py --case eh22
```

Reads `cases/<case>/case.yaml` → generates `cases/<case>/input/{base,perturbed}_{pres,surf}.nc` and `nonrad_forcing.nc`.

### Step 2: Run CFRAM (one command per case)

```bash
python3 run_case.py eh13
python3 run_case.py eh22
```

This extracts input, runs parallel CFRAM decomposition (all CPUs by default, ~20 min/case on 80 cores), and plots results. Output in `cases/eh13/output/` and `cases/eh13/figures/`.

Or run individual steps:

```bash
python3 run_case.py eh13 --step build     # only build ERA5+MERRA-2 input
python3 run_case.py eh13 --step extract   # only extract to Fortran binary
python3 run_case.py eh13 --step run       # only run CFRAM
python3 run_case.py eh13 --step plot      # only plot
```

### Validation against Wu et al. paper_data

```bash
python3 scripts/validate_vs_paper.py eh22           # surface dT comparison
python3 scripts/plot_fig3_independent.py eh22       # 2-column Fig.3 (paper | indep)
```

### Adding a new case

1. Create `cases/my_case/case.yaml` (see `cases/eh13/case.yaml` for template — includes `source:` block for ERA5/MERRA-2 driving)
2. Run: `python3 scripts/build_case_input.py --case my_case` to generate inputs
3. Run: `python3 run_case.py my_case`

## CESM2 4×CO2 Workflow (CMIP6 raw → pyCFRAM)

In addition to ERA5+MERRA-2 (extreme-event attribution), pyCFRAM can decompose any CMIP6 abrupt-4×CO2 vs piControl experiment. The shipped `cesm2_4xco2_official_17p_fu` case is an apple-to-apple match against the OLD Fortran CFRAM benchmark on a 17-plev grid using the Fu radiation engine.

```bash
# 1. Build pyCFRAM input from CMIP6 raw (piControl + abrupt-4xCO2)
python3 scripts/build_cesm2_official.py --case cesm2_4xco2_official

# 2. Re-grid to 17 plev (matches OLD CFRAM exactly)
python3 scripts/subset_to_17p.py cesm2_4xco2_official cesm2_4xco2_official_17p_fu

# 3. Inject CESM 1850 climatology O3 (Phase A: held identical in base & warm)
python3 scripts/inject_cesm_o3.py --case cesm2_4xco2_official_17p_fu

# 4. Mask sub-surface layers below ps_warm (q/o3 = HOLD, cliq/cice/aer = 0)
python3 scripts/mask_subsurface_layers.py --case cesm2_4xco2_official_17p_fu

# 5. Build Fu binary + run on hqlx204 (256 cores)
ssh mini 'ssh lzhenn@hqlx204 "cd /home/lzhenn/work/ust-jumper/pyCFRAM/fortran && \
  source /home/lzhenn/.bashrc_liquor_i22wrf415 && make fu"'
ssh mini 'ssh lzhenn@hqlx204 "cd /home/lzhenn/work/ust-jumper/pyCFRAM && \
  source /home/lzhenn/.bashrc_liquor_i22wrf415 && \
  python3 -u run_case.py cesm2_4xco2_official_17p_fu --step run --nproc 200"'

# 6. Plot 13-panel global / north-polar / south-polar
python3 scripts/plot_13panel_global.py cesm2_4xco2_official_17p_fu
python3 scripts/plot_13panel_polar.py  cesm2_4xco2_official_17p_fu
```

**Validation against the Fortran CFRAM benchmark (collaborator's CESM2 run):**
single-column dual-MC validation at (159, 144) reproduces `partial_T_1.grd` for all 13 dT terms within ≤ 0.003 K when fed identical `.dat` inputs. Global-field discrepancy in DYN/ATM/OCH panels vs OLD's `north.jpg` is **not** a code bug — it is traceable to the OLD reference run using a corrupted O3 input (the collaborator's `o3_base.dat` and `o3_warm.dat` are byte-identical copies of `hus_base.dat` (md5 `1cefb325...`); see `session_log.md` for the full diagnosis).

## CliMLab Single-Column RCE Validation

Idealized clear-sky radiative-convective equilibrium decomposition for sanity-checking pyCFRAM's CFRAM math against a known-correct reference. A 1×1 column where ΔTs ≈ +4.59 K (4×CO2 vs 1×CO2 Manabe-fixed RH) can be verified to close energy balance to within ~0.2 K mid-trop with the higher-order CFRAM options enabled.

```bash
# 1. Run climlab RCE at 1×CO2 (348 ppm) and 4×CO2 (1392 ppm). Writes climlab
#    native NetCDFs + pyCFRAM standard input directly into both case dirs
#    (climlab default 30-level sigma grid, no Fortran rebuild required).
/path/to/conda/python experiments/climlab_validation/run_rce_4xco2.py

# 2. Decompose with each radiation engine (single-cell auto-detected → 1.5-2 s)
python3 run_case.py climlab_4xco2     --step run     # RRTMG
python3 run_case.py climlab_4xco2_fu  --step run     # Fu

# 3. Visualise vertical profile (4-panel: Σ closure + co2 + q + dry)
python3 scripts/plot_singlecol_profile.py climlab_4xco2
```

**Closure verification** (`dT_co2 + dT_q + dT_dry ≈ dT_observed` in mid-trop, RMS over 450–650 hPa):

| Configuration                                                               | mid-trop RMS | upper-trop RMS |
|-----------------------------------------------------------------------------|:------------:|:--------------:|
| 1st-order CFRAM (baseline)                                                  | 0.45 K       | 0.07 K         |
| + `drdt_eval: midstate` (single-variable 2nd-order)                         | 0.45 K       | 0.07 K *       |
| + `co2_handling: midstate` (cancels ∂²R/∂T∂C cross-term)                    | 0.24 K       | 0.04 K         |
| + `drdt_probe: centered` (cancels R_TT in FD Planck probe)                  | **0.21 K**   | **0.02 K**     |

\* midstate Planck alone resolves the stratosphere residual (2.3 → 0.7 K) but not mid-trop, which is dominated by the CO2-T cross-term (residue absorbed once `co2_handling: midstate` is added).

The remaining ~0.2 K mid-trop residual is the **physical floor of 1st-order CFRAM** for clear-sky RCE: climlab convective-adjustment mass redistribution + R_TTT high-order + Clausius–Clapeyron quadratic terms that cannot be cancelled within a one-step Taylor expansion.

The two engines differ by ~10-15 % on individual `dT_co2` / `dT_q` magnitudes (different RT solvers) but agree on the total ECS-equivalent response to within numerical noise — a useful cross-check independent of the OLD CFRAM reference.

## Project Structure

```
pyCFRAM/
├── core/                # Python modules: config, aerosol optics, heating profile, …
├── data/                # Source loaders: era5_source.py, merra2_aerosol.py,
│                        #   cesm2_cmip6_source.py (CMIP6 hybrid→plev)
├── fortran/             # CFRAM-RRTMG + Fu radiation engines
│   ├── Fu/              #   Fu source (cas_fu_radiation.f, fu_helpers.f, para.file)
│   ├── rrtmg_*/         #   upstream RRTMG LW/SW
│   ├── cfram_rrtmg.f90  #   RRTMG full-field driver (legacy, mostly obsolete)
│   ├── cfram_rrtmg_1col.f90  # RRTMG single-column driver (production)
│   ├── cfram_fu_1col.f90 #   Fu single-column driver (production, dual-MC)
│   └── makefile         #   intel (ifort+MKL) + gnu (gfortran+conda LAPACK) toolchains
├── scripts/             # Workflow + diagnostic scripts (see Key Scripts below)
├── plotting/            # Matplotlib/Cartopy visualization helpers
├── cases/               # Case configurations (case.yaml; large NetCDF inputs gitignored)
│   ├── eh13/  eh22/                         # ERA5+MERRA-2 extreme-heat events
│   ├── cesm2_4xco2_official_17p_fu/         # CMIP6 4×CO2, Fu, 17-plev (vs OLD CFRAM)
│   ├── cesm2_4xco2_official_fu/             # CMIP6 4×CO2, Fu, 19-plev
│   ├── cesm2_4xco2_official/                # CMIP6 4×CO2, RRTMG, 19-plev
│   ├── cesm2_4xco2_fu/  cesm2_4xco2/        # collaborator-style 37-plev variants
│   ├── climlab_4xco2/  climlab_4xco2_fu/    # CliMLab 1×→4×CO2 clear-sky RCE (single column)
│   └── climlab_solar/  climlab_solar_fu/    # CliMLab solar-constant perturbation (single column)
├── experiments/         # Standalone validation experiments (climlab RCE etc.)
├── configs/             # Default parameters (defaults.yaml)
└── docs/                # Algorithm + input specifications + technical notes (zh / en)
```

## Key Scripts

### Pipeline (run_case.py orchestrates these)

| Script | Purpose |
|--------|---------|
| `run_case.py` | Unified entry point: build → extract → run → plot |
| `scripts/build_case_input.py` | ERA5 + MERRA-2 → standard NetCDF input (driven by case.yaml) |
| `scripts/build_cesm2_official.py` | CMIP6 piControl + abrupt-4×CO2 → standard NetCDF input |
| `scripts/run_parallel_python.py` | Parallel CFRAM decomposition (multiprocessing, all radiation engines) |

### CESM2 / CMIP6 helpers

| Script | Purpose |
|--------|---------|
| `scripts/subset_to_17p.py` | Re-grid 19-plev → 17-plev (drops k=1, k=5 hPa to match OLD CFRAM exact grid) |
| `scripts/inject_cesm_o3.py` | Replace O3 with CESM 1850 climatology (Phase A: base = warm, frc_o3 = 0) |
| `scripts/mask_subsurface_layers.py` | Mask layers below ps_warm (q/o3 HOLD, cliq/cice/aer → 0) |
| `scripts/expand_nonrad_to_column.py` | Distribute column-integrated lhflx/shflx to per-layer profiles |
| `scripts/compare_cesm2_official_vs_collaborator.py` | Side-by-side 13-panel diagnostic |

### Plotting

| Script | Purpose |
|--------|---------|
| `scripts/plot_fig3_independent.py` | Fig.3 (Wu et al.): 2-column paper-vs-independent spatial maps |
| `scripts/plot_fig3_self.py` / `plot_fig3.py` / `plot_fig4.py` / `plot_fig5.py` | Other paper figures |
| `scripts/plot_13panel_global.py` | 13-panel global surface dT decomposition (matches OLD CFRAM layout) |
| `scripts/plot_13panel_polar.py` | 13-panel north + south polar dT decomposition |

### Verification / diagnostics

| Script | Purpose |
|--------|---------|
| `scripts/verify_phase1_spc.py` | Per-species optical-property additivity sanity check |
| `scripts/verify_phase4_spc.py` | Per-species aerosol forcing/dT additivity + regression |
| `scripts/verify_phase5_vs_paper.py` | Per-species spatial correlation vs Wu et al. paper_data |
| `scripts/verify_cloud_split.py` | Cloud LW/SW additivity (`bulk == lw + sw`) + regression |
| `scripts/validate_vs_paper.py` | Surface dT comparison vs Wu et al. results |
| `scripts/diag_cloud_column.py` / `diag_drdt_singlecol.py` | Single-column rad/Planck-matrix dumps for debugging |

### Lapse-Rate kernel module

| Script | Purpose |
|--------|---------|
| `scripts/compute_lr_kernel.py` | M2: native ΔR_LR/ΔR_PL from `dT_observed`, per configured kernel (default CloudSat=Kramer + GFDL) |
| `scripts/validate_lr_vs_climkern.py` | M2: cross-validate against ClimKern's `calc_T_feedbacks` (needs `pycfram-kern` env) |
| `scripts/plot_lr_comparison.py` | M2: native/ClimKern/diff triptych + Kramer-vs-GFDL difference maps |
| `scripts/compute_lr_attribution.py` | M3: per-physical-process (`dT_q`, `dT_atmdyn`, ...) lapse-rate attribution |
| `scripts/plot_lr_attribution.py` | M3: per-process attribution maps + zonal-mean profile |
| `data/kernel_source.py` | Stage CloudSat/GFDL kernel NetCDFs from the climkern install into `data/kernels/` |

## Decomposition Output

`cases/<case>/output/cfram_result.nc` contains partial temperature changes `dT_*` and the underlying radiative forcings `frc_*` (W/m²), shape `(lev, lat, lon)` with surface at `lev[-1]`:

| Group | Terms |
|-------|-------|
| Radiative (bulk) | `co2`, `q`, `ts`, `o3`, `solar`, `albedo`, `cloud`, `aerosol`, `warm` |
| Cloud LW/SW split | `cloud_lw`, `cloud_sw` — exact additive (`cloud == cloud_lw + cloud_sw`) |
| Aerosol species | `bc`, `ocphi`, `ocpho`, `sulf`, `ss`, `dust` — sum ≈ bulk (small non-linear residual) |
| Non-radiative | `lhflx`, `shflx` (forcing placed on surface row only) |
| Derived | `atmdyn`, `sfcdyn`, `ocndyn`, `dry`, `observed` — `dry = atmdyn + sfcdyn`, `sfcdyn = ocndyn + lhflx + shflx` (exact) |

## Lapse-Rate Kernel Module (M2/M3)

A radiative-kernel Lapse-Rate/Planck decomposition, cross-validated against [ClimKern](https://github.com/tyfolino/climkern)'s `calc_T_feedbacks` — see [docs/m2_kernel_module.md](docs/m2_kernel_module.md) for the full write-up and [docs/m3_methodology_comparison.md](docs/m3_methodology_comparison.md) for the CFRAM-vs-kernel-method comparison + per-process attribution.

```bash
# Native module only (numpy/scipy/netCDF4, no xesmf needed):
python scripts/compute_lr_kernel.py cesm2_4xco2_official
python scripts/compute_lr_attribution.py cesm2_4xco2_official   # M3 route i: per-process attribution

# Cross-validation against ClimKern (needs the `pycfram-kern` conda env, see docs/m2_kernel_module.md §Environment):
conda activate pycfram-kern
python scripts/validate_lr_vs_climkern.py cesm2_4xco2_official
python scripts/plot_lr_comparison.py cesm2_4xco2_official

# Or via run_case.py:
python run_case.py cesm2_4xco2_official --step lr
python run_case.py cesm2_4xco2_official --step lr-attr
```

## Phase 3: DAMIP Multi-Model Support

Beyond ERA5/MERRA-2 and the single CESM2 4×CO2 dataset above, pyCFRAM can
decompose **any CMIP6 DAMIP single-forcing experiment** (`hist-aer`,
`hist-GHG`, `hist-nat`, `hist-stratO3`, ...) from **any** contributing
model — a registered `cmip6_damip` data source
(`data/cmip6_damip_source.py`) auto-detects each model's calendar, hybrid
vertical-coordinate naming, native horizontal grid, and published pressure
levels, and applies a documented missing-variable decision tree (cloud/O₃/
solar/aerosol all have well-defined fallbacks) so that models missing
cloud or ozone data — the common case, not the exception, across CMIP6 —
still produce a physically valid decomposition. Base/perturbed states are
the experiment's first-decade (1850–1859) vs last-available-decade
climatology; a from-scratch, authentication-free ESGF download client
(`data/esgf_fetch.py` + `scripts/download_damip.py`) handles data
acquisition. As with the rest of pyCFRAM, `core/` and `fortran/` are
untouched by any of this — DAMIP support is entirely a new data-source
plugin feeding the existing generic writer and CFRAM engine.

Nine models are validated end-to-end against real downloaded ESGF data:
IPSL-CM6A-LR, MRI-ESM2-0, CESM2, CNRM-CM6-1, MIROC6, GISS-E2-1-G,
HadGEM3-GC31-LL, and CanESM5 (the M4/M5 model set), plus NorESM2-LM as a
worked "custom model accession" example showing that adding a model not in
that list requires touching only a `configs/damip_models.d/<model>.yaml`
and a `cases/<case>/case.yaml` — no Python changes. Each has a
corresponding case directory:

```
cases/damip_{ipsl,mri,cesm2,cnrm,miroc6,giss,hadgem3,canesm5,noresm2}_histaer/
```

```bash
# Build + run any DAMIP case exactly like an ERA5/CESM2 case:
python3 run_case.py damip_ipsl_histaer --step build
python3 run_case.py damip_ipsl_histaer --step run --nproc 200
# Every run writes a human-readable cases/<case>/output/<case>.summary.txt
# recording which processes were ACTIVE/SKIPPED/CLIMATOLOGY/ANALYTIC and why.
```

See [docs/plan_ph3.md](docs/plan_ph3.md) for the full Phase 3 execution
plan (DAMIP protocol, cross-model heterogeneity strategy, ESGF data
availability matrix), [docs/m4_damip_module.md](docs/m4_damip_module.md)
for the module's architecture and the concrete cross-model bugs that
shaped it, and [docs/m5_multimodel_userguide.md](docs/m5_multimodel_userguide.md)
for the per-model support matrix, known-issues log, and the worked
new-model-accession example.

## Input Data Format

See [docs/input_spec.md](docs/input_spec.md) for the standard NetCDF input format. pyCFRAM takes any two atmospheric states (base, perturbed) — it does not prescribe how they are defined.

## References

**CFRAM methodology**
- Lu, J. and Cai, M., 2009. A new framework for isolating individual feedback processes in coupled general circulation climate models. Part I: Formulation. *Climate Dynamics*, 32, 873–885.
- Cai, M. and Lu, J., 2009. A new framework for isolating individual feedback processes in coupled general circulation climate models. Part II: Method demonstrations and comparisons. *Climate Dynamics*, 32, 887–900.

**CFRAM-A: RRTMG with aerosol**
- Zhang, T., Deng, Y., Chen, J., Yang, S., Gao, P. and Zhang, H., 2022. Disentangling physical and dynamical drivers of the 2016/17 record-breaking warm winter in China. *Environmental Research Letters*, 17(7), 074024. https://doi.org/10.1088/1748-9326/ac79c1

**CFRAM-A applications (extreme-event attribution)**
- Wu, Q., Li, Q., Zhang, T., Sun, X., Yang, S. and Hu, X., 2025. Quantitative attribution of 2013 and 2022 extremely hot events in China: insights from a climate feedback–response perspective. *Journal of Climate*, 38(17), pp.4331–4349.