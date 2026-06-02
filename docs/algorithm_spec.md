# CFRAM-A Algorithm Specification

**Scope:** Fortran radiation engines in pyCFRAM (RRTMG: `fortran/cfram_rrtmg.f90`, single-column `cfram_rrtmg_1col.f90`; Fu: `fortran/cfram_fu_1col.f90` with sources in `fortran/Fu/`). The two engines are interchangeable per case via `case.yaml`'s `run.executable` field — see README §"Selecting the radiation engine". Sections 1–10 below describe the RRTMG path; §12 documents the Fu-specific dual-MC sub-column overlap protocol added 2026-05-10.

Python decomposition (Planck solve, non-radiative, residuals) is described separately in [`docs/technical_notes_zh.md`](technical_notes_zh.md).

**References:**
- Lu, J. and Cai, M. (2009). *Climate Dynamics*, 32, 873–885. [CFRAM Part I]
- Cai, M. and Lu, J. (2009). *Climate Dynamics*, 32, 887–900. [CFRAM Part II]
- Zhang, T., Deng, Y., Chen, J., Yang, S., Gao, P. and Zhang, H. (2022). *Environmental Research Letters*, 17(7), 074024. https://doi.org/10.1088/1748-9326/ac79c1 
- Wu, Q., Li, Q., Zhang, T., Sun, X., Yang, S. and Hu, X. (2025). *Journal of Climate*, 38(17), 4331–4349. 

## 1. Overview

CFRAM decomposes total temperature change ΔT into partial contributions from individual physical processes:

```
ΔT_i = −(∂R/∂T)⁻¹ × F_i      (one column per process)
ΔT_total = Σ_i ΔT_i          (linearised)
```

The Fortran engine emits the radiative forcings `F_i` (W/m²) and the Planck response matrix `∂R/∂T`; the Python runner applies the inverse. Non-radiative and residual terms (LH/SH flux, sfcdyn, ocndyn, atmdyn) are computed in Python using the same Planck matrix or energy-conservation residuals — see `scripts/run_parallel_python.py`.

Current outputs of the Fortran engine per grid point:

- 9 bulk forcings: `co2, q, ts, o3, solar, albedo, cloud, aerosol, warm`
- 2 cloud components: `cloud_lw, cloud_sw` (exactly additive: `cloud == cloud_lw + cloud_sw`)
- 6 per-species aerosol: `bc, ocphi, ocpho, sulf, ss, dust` (approximately additive; residual ≡ aerosol–aerosol / aerosol–cloud non-linear coupling)
- Planck inverse matrix `(∂R/∂T)⁻¹` of shape `(nlayer+1, nlayer+1)` (atmosphere + surface)

## 2. Call Flow

```
cfram_rrtmg.f90 (main)
├── Read binary inputs (bulk):     t, q, o3, ts, ps, solin, ssrd, ssru, cloud,
│                                  co2, aerosol bulk (aod_lw, aod_sw, ssa_sw, g_sw)
├── Read binary inputs (per-spc):  aerosol_{aod_lw,aod_sw,ssa_sw,g_sw}_{base,warm}_spc.dat
│                                  (species is innermost axis on disk)
├── Derived quantities: albedo_sw = ssru/ssrd, zenith = solin/scon
├── Determine active layers: nlayer (plev < ps)
│
├── rad_driver(base)     → rad_1d_base, lw_1d_base, fd/u_1d_base, fdl/ful_1d_base
├── rad_driver(warm)     → rad_1d_warm, fd/u_1d_warm  [ALL vars → warm]
│
├── calc_drdt(base)      → drdt(nlayer+1, nlayer+1)      [LW-only Planck matrix]
├── drdt_inv = inv(drdt(1:nlayer+1, 1:nlayer+1))         [atm + surface]
│
├── rad_driver(co2)      → rad_1d_co2
├── rad_driver(q)        → rad_1d_q
├── rad_driver(ts)       → rad_1d_ts
├── rad_driver(o3)       → rad_1d_o3
├── rad_driver(solar)    → rad_1d_solar
├── rad_driver(albedo)   → rad_1d_albedo
├── rad_driver(cloud)    → rad_1d_cloud, fd/u_1d_cloud
│    └── snapshot: lw_1d_cloud, sw_1d_cloud, fdl/ful/fds/fus_1d_cloud
├── rad_driver(aerosol)  → rad_1d_aerosol                [all aerosols → warm]
│
├── do isp = 1, nspecies (6)                             [per-species aerosol]
│   ├── build mix state: swap only species[isp] base → warm using AOD-weighted
│   │   combination for ssa and g
│   ├── rad_driver(mix) → rad_1d_spc(isp)
│   └── frc_spc(isp) = rad_1d_spc − rad_1d_base
│
├── Compute 17 layer-resolved forcings
│   ├── frc_X(1:nlayer)   = rad_perturbed_X(1:nlayer) − rad_base(1:nlayer)
│   ├── frc_X(nlayer+1)   = [fd_X(sfc) − fu_X(sfc)] − [fd_base(sfc) − fu_base(sfc)]
│   ├── Cloud split uses lw_1d_cloud / sw_1d_cloud and corresponding surface
│   │   fluxes (SW base derived from total − LW base; no extra base storage)
│
├── Write drdt_inv.dat and frc_*.dat (9 bulk + 2 cloud_lw/sw + 6 per-species)
│
└── Python runner reads frc + drdt_inv and computes dT = −drdt_inv · frc
    for every term (see `scripts/run_parallel_python.py`).
```

Total rad_driver calls per grid point: 2 (base/warm) + 8 (partial perturbations) +
(nlayer+1) for Planck + 6 (per-species) ≈ 56. Per-species adds ~13% overhead.

## 3. Coordinate Systems

| System | Index 1 | Index nlayer | Notes |
|--------|---------|-------------|-------|
| ERA5 input | TOA (1 hPa) | Surface (1000 hPa) | Top-to-bottom |
| RRTMG internal | Surface | TOA | Bottom-to-top |
| Output | TOA | Surface | Back to ERA order |

**Transformation:** `rrtmg_idx = nlayer - era_idx + 1`

**Interface pressures:**
- `pint(0) = 0.5 hPa` (TOA boundary)
- `pint(i) = (plev(i) + plev(i+1)) / 2` for i=1..nlayer-1
- `pint(nlayer) = ps` (surface)

## 4. Key Physical Computations

### 4.1 Column Density (molecules/cm²)

```
Mmoist = (1 - qv) × 28.966 + qv × 18.016    [g/mol, moist air molecular weight]
coldry = (Δp × 1e3 × 6.02e23) / (1e2 × 9.8 × Mmoist × (1 + qv))
```

where qv = specific humidity converted to volume mixing ratio = q × 28.966/18.016

### 4.2 Gas Mixing Ratios (wkl array, 7 species)

| Index | Species | Source | Conversion |
|-------|---------|--------|------------|
| 1 | H₂O | ERA5 specific humidity | q × 28.966/18.016 → volume mixing ratio |
| 2 | CO₂ | co2.txt (ppmv) | × 1e-6 |
| 3 | O₃ | ERA5 ozone mass mixing ratio | × 28.966/48.0 |
| 4 | N₂O | hardcoded 0.28 ppmv | × 1e-6 |
| 5 | CO | not used (0) | — |
| 6 | CH₄ | hardcoded 1.6 ppmv | × 1e-6 |
| 7 | O₂ | constant 0.209 | volume fraction |

**Final column amount:** `wkl_col = coldry × mixing_ratio` [molecules/cm²]

### 4.3 Cloud Water Path (g/m²)

```
coef = (1/9.8) × 100 × 1000 = 10204.08
ciwp = coef × cldiwc × Δp    [where Δp in hPa, cldiwc in kg/kg]
clwp = coef × cldlwc × Δp
```

Cloud effective radii: ice = 5 μm, liquid = 20 μm (hardcoded)

### 4.4 Surface Albedo

```
albedo_sw = swup_surf / swdn_surf    (0 if swdn=0)
albedo_lw = 0.0                       (hardcoded)
```

### 4.5 Solar Zenith Angle

```
zenith = solin / scon    (cosine of zenith angle)
```

where solin = TOA incident solar radiation (W/m²), scon = 1360.98 W/m²

## 5. Radiative Flux Definitions

For each atmospheric layer i (1 to nlayer):

```
sw_base(i) = fds(i) + fus(i+1) - fds(i+1) - fus(i)    [SW convergence]
lw_base(i) = fdl(i) + ful(i+1) - fdl(i+1) - ful(i)    [LW convergence]
rad_base(i) = lw_base(i) + sw_base(i)                    [total convergence]
```

Flux indices: 1 = TOA boundary, nlayer+1 = surface boundary (in RRTMG bottom-to-top order, but code handles the reversal)

## 6. Forcing Computation

For atmospheric layers (1:nlayer):
```
frc_X(1:nlayer) = rad_perturbed_X(1:nlayer) - rad_base(1:nlayer)
```

For surface (nlayer+1):
```
frc_X(nlayer+1) = [fd_perturbed(sfc) - fu_perturbed(sfc)] - [fd_base(sfc) - fu_base(sfc)]
```

## 7. Planck Matrix (drdt)

**Dimension:** (nlayer+1) × (nlayer+1)  
**Method:** LW-only, +1K perturbation per layer

```
For each layer k = 1..nlayer:
    T_perturbed = T_base
    T_perturbed(k) += 1.0 K
    Call rad_driver_lw → lw_1k, fdl_1k, ful_1k
    drdt(k, 1:nlayer) = lw_1k - lw_base
    drdt(k, nlayer+1) = [fdl_1k(sfc) - ful_1k(sfc)] - [fdl_base(sfc) - ful_base(sfc)]

For surface:
    Ts_perturbed = Ts_base + 1.0 K
    Call rad_driver_lw → lw_1k, fdl_1k, ful_1k
    drdt(nlayer+1, 1:nlayer) = lw_1k - lw_base
    drdt(nlayer+1, nlayer+1) = [fdl_1k(sfc) - ful_1k(sfc)] - [fdl_base(sfc) - ful_base(sfc)]
```

**Indexing convention:** drdt(row, col) where row = radiative response location, col = T perturbation location

## 8. Temperature Response

The full Planck matrix on `(1:nlayer+1, 1:nlayer+1)` (atmosphere + surface) is inverted and applied in Python:

```
drdt_inv = inv(drdt(1:nlayer+1, 1:nlayer+1))

dT_X = −drdt_inv · frc_X       (length nlayer+1 vector, atm[1:nlayer] + surface[nlayer+1])
```

Both atmospheric layer and surface responses are solved together. Surface dT is NOT zeroed out (this is a change from the original CFRAM-RRTMG reference code). Levels above the surface pressure (nlayer+1:nlev) are filled with `-999.0` in the output.

### 8.1 Dynamics / non-radiative terms (Python)

The full warm-state radiative imbalance `frc_full` (from the `rad_1d_full`
call, using `t_warm/q_warm/o3_warm/...`) drives the dynamics decomposition via
the same Planck inverse:

| Term | Forcing vector | Meaning |
|------|---------------|---------|
| `dT_atmdyn` | `−frc_full` on atmosphere rows, surface = 0 | atmospheric heat transport |
| `dT_sfcdyn` | `−frc_full` on surface row, atmosphere = 0 | surface energy-budget residual |
| `dT_dry`    | `−frc_full`, full column | `= dT_atmdyn + dT_sfcdyn` |
| `dT_ocndyn` | surface row `= −frc_full[sfc] − Δlh − Δsh` | ocean / sub-surface transport |

`lhflx`/`shflx` forcings are placed on the **surface row only**, so

```
dT_sfcdyn ≡ dT_ocndyn + dT_lhflx + dT_shflx     (exact to machine precision)
```

i.e. `sfcdyn` is by construction the sum of the ocean-transport and surface-flux
responses — not an independent process. (Placing the flux forcing on atmosphere
rows breaks this identity; that bug was fixed 2026-06-01.) For the OLD-CFRAM-style
closure *plot* only, `plot_closure_profile.py` redefines `atmdyn` as the closure
residual `dT_obs − Σrad − Σnonrad − dT_ocndyn`; the NetCDF `dT_atmdyn` variable
remains the forcing-based quantity above.

## 9. Aerosol Optical Properties

**Framework origin:** the extension of CFRAM to include an explicit aerosol
partial-radiative-perturbation term — i.e. the **CFRAM-A** framework that
pyCFRAM implements — was first introduced by Zhang et al. (2022, *ERL* 17:074024,
https://doi.org/10.1088/1748-9326/ac79c1). That paper notes "the effect of
aerosols ... has not been included in the previous CFRAM analysis" and
adds a 5-species aerosol term (BC, OC, sulfate, sea salt, dust) computed via
off-line RRTMG v5. pyCFRAM inherits this conceptual framework; relative to
the Fortran + NCL reference pipeline, the additions here are on the
orchestration side (automated ERA5 + MERRA-2 preprocessing, Python-side
full Planck-matrix solve, per-grid-point `multiprocessing` parallelism),
not on the radiative decomposition itself. The optical-property source
also differs — see below.

**Input:** MERRA-2 `M2I3NVAER` (13 mass mixing ratios on 72 model levels, kg/kg)
merged into 6 pyCFRAM species — see `configs/defaults.yaml` for the species map.
(Zhang et al. 2022 also drove their aerosols from MERRA-2 but used 5 species.)

**Lookup tables (pyCFRAM):** GOCART `opticsBands_*.RRTMG.nc` under
`fortran/data_prep/aerosol/` (bext/bsca/qext/qsca/g vs relative humidity and
size index, per RRTMG band). This lookup-table pipeline comes from the Wu et al.
(2025) reference code and is **different** from the MAM4/CAM6 optical treatment
used in Zhang et al. (2022); the two choices produce distinct species-level
forcings even given identical MERRA-2 mixing ratios.

**Computation** (implemented in `scripts/extract_full_field.py` and
`scripts/run_parallel_python.py:compute_aerosol_column`, ported from the
original Matlab `aer_opt_*.m`):

```
For each species (bc, ocphi, ocpho, sulf, ss, dust):
    1. kext(band, layer) = bext[r_idx, rh_idx_of_layer, band]  (nearest RH index)
    2. aod_species(band, layer) = mixing × dens × 1e3 × layer_thickness × kext × 1e-3
    3. ssa_species = qsca / qext
    4. g_species from lookup table
Bulk (AOD-weighted):
    aod_bulk = Σ aod_species
    ssa_bulk = Σ (aod_species × ssa_species) / aod_bulk
    g_bulk   = Σ (aod_species × ssa_species × g_species) / (aod_bulk × ssa_bulk)
```

**Bands:** SW 1–14 (jpband), LW 15–30 (nbndlw); 14 SW + 16 LW = 30 RRTMG bands.

**Per-species perturbation mix** (inline in `cfram_rrtmg.f90`, Phase 3): for each
species isp, rebuild an intermediate state where only species isp is warm and
the other five remain at base, using the same AOD-weighted combination. The
response `rad_1d_spc(isp)` is then differenced against `rad_1d_base` to give
`frc_<species>`. This preserves aerosol–aerosol and aerosol–cloud optical
non-linearity that a linear (path-B) Planck-inverse split would drop.

## 10. Decomposition Coverage vs Wu et al. (2025) Paper

| Paper (10 terms) | pyCFRAM output | Note |
|---|---|---|
| ΔS^solar | `frc_solar` | ✓ |
| Δ(S−R)^CO₂ | `frc_co2` | ✓ |
| Δ(S−R)^CH₄ | — | CH₄ hardcoded at 1.6 ppmv; no separate perturbation |
| Δ(S−R)^O₃ | `frc_o3` | ✓ |
| Δ(S−R)^WV | `frc_q` | ✓ |
| Δ(S−R)^CLD | `frc_cloud` | ✓ |
| — (paper aggregates) | `frc_cloud_lw`, `frc_cloud_sw` | Extra: exact LW/SW split |
| Δ(S−R)^AER | `frc_aerosol` (bulk) | ✓ |
| — (paper per-species via external forcing) | `frc_bc`, `frc_ocphi`, `frc_ocpho`, `frc_sulf`, `frc_ss`, `frc_dust` | Extra: in-RRTMG 6-species perturbation |
| ΔS^α | `frc_albedo` | ✓ |
| ΔQ^ATD | `dT_atmdyn` (Python) | Forcing-based from `frc_full` atm rows (see §8.1) |
| ΔQ^SRF | `dT_sfcdyn`, `dT_ocndyn`, `dT_lhflx`, `dT_shflx` (Python) | `sfcdyn = ocndyn + lhflx + shflx` exactly; all from `frc_full` surface row + input fluxes (see §8.1) |
| — | `dT_dry` | Extra: `= atmdyn + sfcdyn` (full-column `frc_full` response) |
| — | `frc_ts` | Extra: surface temperature radiative response |

## 11. Historical Note

The reference Fortran code was originally driven by NCL preprocessing scripts
(`OA/202604020900/CFRAM_RRTMG/`, removed from this repo). The current pipeline
is fully Python (ERA5 + MERRA-2 via `scripts/build_case_input.py`, see
[`technical_notes_zh.md`](technical_notes_zh.md) §5). Any residual references to
`co2.txt`, NCL scripts, or Matlab `aer_opt_*.m` in older commits are obsolete.

## 12. Fu Engine: Dual MC Sub-column Overlap

The Fu engine (`fortran/cfram_fu_1col.f90` + `fortran/Fu/cas_fu_radiation.f`) replaces RRTMG's stochastic in-cloud handling with the OLD Fortran CFRAM Monte Carlo Independent Column Approximation (MCICA): a 100 × 100 (sub-column × layer) integer mask `no_cloud(igrid, l)` (0 = cloudy, 1 = clear) determines which sub-columns are cloudy at each level. The fraction of cloudy sub-columns at level `l` matches `cc(l)`; the vertical correlation matches the max-random overlap convention.

**Two pattern buffers**, mirroring OLD CFRAM (`raw/CFRAM.zip`):

| Pattern | Sampled from | Used by Fortran cases | OLD CFRAM equivalent |
|---|---|---|---|
| `no_cloud_base_pat` | `cc_base` | 0 (base), 2 (co2), 3 (q), 4 (o3), 5 (ts), 6 (solar), 7 (albedo), and the 18 `drdt` perturbations | `base_no_cloud_out_*.dat` (written by GW-base.f, read by GW-co2/wv/o3/solar/albedo/ts/drdt) |
| `no_cloud_warm_pat` | `cc_warm` | 1 (warm), 8 (cloud), 9 (full) | `warm_no_cloud_out_*.dat` (written by GW-warm.f, read by GW-cloud.f) |

Each per-cell run:
1. case 0 generates the base pattern via `S_R_cloudy → ran3` from `cc_base`, saves into `no_cloud_base_pat`.
2. case 1 generates the warm pattern from `cc_warm`, saves into `no_cloud_warm_pat`.
3. Subsequent cases copy the appropriate pattern into the shared `no_cloud_fixed` array (in `/cloud_fix/` common block, read by `S_R_cloudy`) before the radiation call.
4. The `drdt` block explicitly resets to base pattern.

If `data_prep/{base,warm}_no_cloud_seed.dat` are present (100 × 100 INTEGER*4 records), they preempt the in-program generation, enabling bit-perfect reproduction of an external CFRAM run.

**Why this matters.** Pre-2026-05-10 pyCFRAM Fu used a single MC pattern across all cases (saved from case 0). Cases 1/8/9 then radiated the cc_warm fields through a sub-column realisation drawn from cc_base — i.e. `effective_cc ≈ cc_base` while the cloud-water fields filled into "cloudy" sub-columns came from cc_warm. The mis-alignment produced a +1.8 K global-mean LW/SW mirror flip in CLDL/CLDS. Single-column DP validation at (159, 144) on the collaborator's CESM2 dataset, after the fix, reproduces the OLD `partial_T_1.grd` for all 13 dT terms within ≤ 0.003 K.
