# M2: Radiative-Kernel Lapse-Rate Module

Implements the radiative-kernel Lapse-Rate/Planck decomposition specified in
`contract/contract_10.pdf` M2, cross-validated against
[ClimKern](https://github.com/tyfolino/climkern) v1.2's `calc_T_feedbacks`
(Janoski et al., GMD 2025). See `docs/plan.md` for the full execution plan
this module was built against.

## What it computes

Given a base and a perturbed atmospheric state, the total temperature
response `ΔT(p)` can be split into a vertically-uniform "Planck" part
(everything warms/cools like the surface) and a "Lapse-Rate" part (the
deviation from uniform warming):

```
ΔR_LR(lat,lon) = Σ_k  K_lw_t(k,lat,lon) · (ΔT(k) − ΔTS) · dp(k)/10000
ΔR_PL(lat,lon) = K_lw_ts·ΔTS + Σ_k K_lw_t(k,lat,lon) · ΔTS · dp(k)/10000
```

`K_lw_t`/`K_lw_ts` are pre-computed radiative kernels (W m⁻² K⁻¹ per 100 hPa)
that answer "how much does the TOA LW flux change per 1K warming at this
level/at the surface". `dp` is the tropospheric layer thickness (Pa),
computed on the kernel's own pressure grid with the tropopause defined as
`3e4 − 2e4·cos(lat)` Pa (100 hPa at the equator, 300 hPa at the poles) —
this both restricts the integral to the troposphere and automatically
excludes underground layers (the pressure-interface clipping makes their
`dp` exactly zero).

Two kernels are supported out of the box: `CloudSat` (the "Kramer" kernel
in the contract's terminology — an observationally-based kernel from
Kramer et al. 2019) as the primary kernel, and `GFDL` (Soden & Held 2006,
model-based) as a reference for kernel-choice sensitivity.

## Architecture

```
core/kernels.py       KernelSet: read TOA_<name>_Kerns.nc, bilinear regrid
                       (numpy/scipy only — no xesmf/climkern import)
core/lr_kernel.py      tropopause_pa / layer_dp_troposphere / delta_R_lr /
                       delta_R_planck / interp_to_kernel_plev /
                       extract_and_apply — the vertical-integration math
data/kernel_source.py  Locate/stage kernel NetCDFs (reuses climkern's own
                       downloaded data; core/ never imports climkern)
scripts/compute_lr_kernel.py        CLI: case -> lr_kernel.nc
scripts/validate_lr_vs_climkern.py  ClimKern reference + cross-val metrics
                                     (the ONLY script here that imports
                                     xesmf/climkern)
scripts/plot_lr_comparison.py       native/ClimKern/diff maps
```

**Dependency boundary** (deliberate, see `docs/plan.md` §3.3): `core/` and
`compute_lr_kernel.py` only need numpy/scipy/netCDF4 and run on any
machine (including hqlx210, which has neither xesmf nor climkern
installed). Only `validate_lr_vs_climkern.py` needs the full
xesmf/esmpy/climkern stack, isolated in a dedicated conda env.

## Environment (`pycfram-kern`)

```bash
conda create -n pycfram-kern -c conda-forge python=3.11 \
    numpy xarray netcdf4 scipy matplotlib cartopy cf_xarray \
    esmpy xesmf pooch importlib_resources pytest pyyaml -y
conda activate pycfram-kern
pip install climkern remotezip
```

**Kernel data**: climkern's own `python -m climkern download` pulls a
single 5.3 GB zip containing kernels for many other models we don't need.
Instead, `data/kernel_source.py` (once populated, see below) or a manual
`remotezip` extraction pulls only the two required files
(`data/kernels/CloudSat/TOA_CloudSat_Kerns.nc`,
`data/kernels/GFDL/TOA_GFDL_Kerns.nc`, ~130 MB each) directly via HTTP
Range requests against the same Zenodo record climkern's `download.py`
already points at (`doi:10.5281/zenodo.10223376`, which — being a Zenodo
*concept* DOI — resolves to whatever the current version record is; at
the time of writing that's record `18565513`):

```python
from remotezip import RemoteZip
with RemoteZip("https://zenodo.org/records/18565513/files/data.zip") as zf:
    zf.extract("data/kernels/CloudSat/TOA_CloudSat_Kerns.nc", path=dest)
    zf.extract("data/kernels/GFDL/TOA_GFDL_Kerns.nc", path=dest)
```

Stage the extracted files at `data/kernels/<Name>/TOA_<Name>_Kerns.nc`
(gitignored — reproducible from the DOI above) so
`data/kernel_source.get_kernel_path()` finds them without needing
climkern importable at runtime. `data/kernels/manifest.txt` records the
md5 sums.

Zenodo's download servers were observed to be heavily rate-limited from
one network path (~76 KB/s) but much faster from another (~350 KB/s,
observed via a different egress point) — this is a known Zenodo-wide
2026 issue (bot-traffic overload, see the Zenodo blog), not specific to
this record. If a download stalls, retry from a different network path.

## `case.yaml` configuration

```yaml
lapse_rate:
  kernels: [CloudSat, GFDL]     # CloudSat=Kramer(primary); GFDL=reference
  kernel_months: annual         # 'annual' = 12-month kernel mean
  sky: all-sky                  # all-sky | clear-sky
  tropopause: climkern          # make_tropo formula: 3e4 - 2e4*cos(lat) Pa
```

## Cross-validation results (`cesm2_4xco2_official`, global 192×288)

Both native and ClimKern paths are fed the exact same temperature
response (`dT_observed` from `cfram_result.nc`, which by construction
equals `perturbed.ta − base.ta` / `perturbed.ts − base.ts` — the same
quantity ClimKern's `calc_T_feedbacks` is given via the case's raw
`base_pres/base_surf`/`perturbed_pres/perturbed_surf` NetCDFs, tiled to
12 identical months since ClimKern requires `len(time) % 12 == 0`).

| Kernel | corr (area-weighted) | domain-mean rel. diff | Contract gate (0.85 / 15%) |
|---|---|---|---|
| CloudSat (Kramer) | **0.9997** | **1.32%** | PASS |
| GFDL | **0.9999** | **0.02%** | PASS |

Both far exceed the contract's 0.85/15% gate and the internal green line
(corr ≥ 0.98, rel_diff ≤ 5%) — see `cases/cesm2_4xco2_official/output/lr_validation_report.txt`
and `cases/cesm2_4xco2_official/figures/fig_lr_comparison_{CloudSat,GFDL}.png`
for the full triptych (native / ClimKern / diff).

Reproduce:
```bash
python scripts/compute_lr_kernel.py cesm2_4xco2_official         # any env
conda activate pycfram-kern
python scripts/validate_lr_vs_climkern.py cesm2_4xco2_official
python scripts/plot_lr_comparison.py cesm2_4xco2_official
```

## Known limitations

- A thin band at the extreme poles (~2° from each pole) shows saturated
  colorbar values in the native-vs-ClimKern diff maps — an edge-clipping
  artifact where the kernel's own native grid (2°×2.5°, not truly
  polar-covering) is extrapolated. It affects a negligible fraction of
  global area and does not move the area-weighted correlation below
  0.999.
- Regional (non-global) cases: the kernel is always regridded from its
  native *global* grid, so periodic longitude wrapping is always
  correct regardless of the target case's own extent (no special
  handling needed — see `docs/plan.md` §2.3 "关于 periodic").
- `KernelSet.select('annual', ...)` on a kernel level where every one of
  the 12 months is masked (e.g. a plev that never exists in the kernel's
  own underlying model's topography at that lat/lon) raises a benign
  numpy `RuntimeWarning: Mean of empty slice` and yields NaN, which is
  then correctly treated as a zero contribution by `delta_R_lr`'s
  `np.nansum` — this is expected, not a bug.

## Gotcha fixed during implementation

`netCDF4` returns masked arrays for variables with a `_FillValue`
attribute. `np.array(masked_array)` (without going through
`np.ma.filled`) silently replaces masked (e.g. underground) entries with
the *raw* fill value (~1e36) instead of NaN — this blew up `ΔR_LR` to
~1e40 W/m² before being caught (`core/kernels.py` now uses
`np.ma.filled(..., np.nan)` throughout).

## Tests

```bash
pytest tests/test_kernels.py tests/test_lr_kernel.py tests/test_e2e_smoke.py -v
# xesmf cross-check (only runs inside pycfram-kern):
conda activate pycfram-kern && pytest tests/test_kernels.py::test_regrid_vs_xesmf -v
```

`core/kernels.py`, `core/lr_kernel.py` and `core/lr_attribution.py` are
at 94–100% line coverage (`pytest --cov=core`). `tests/data/smoke/`
holds an 8×8 lat/lon crop of `cesm2_4xco2_official`'s `cfram_result.nc`
plus coarse (~0.5 MB) subsampled CloudSat/GFDL kernel files, driving the
end-to-end `compute_lr_kernel.py` CLI test in `tests/test_e2e_smoke.py`.
