# M4: CMIP6 DAMIP Multi-Model Data Source

Technical documentation of the M4 DAMIP module — the module that lets pyCFRAM
consume **any** CMIP6 DAMIP (Detection & Attribution MIP) single-forcing
experiment (`hist-aer`, `hist-GHG`, `hist-nat`, `hist-stratO3`, ...) from
**any** contributing model, without touching the CFRAM engine
(`core/`, `fortran/`). See `docs/plan_ph3.md` for the full execution plan
this module was built against (§1–§7 cover the design rationale in more
depth than repeated here); this document is the as-built reference plus the
concrete bugs that shaped it.

## What it does

Given `model` + `experiment` + a raw-data directory of downloaded CMIP6
NetCDF files, the module:

1. Discovers the right files (variant, grid label, time-chunked or
   single-file) for a fixed list of Amon variables.
2. Builds a day-weighted monthly climatology for the experiment's first
   decade (`base_years`, default 1850–1859) and last available decade
   (`warm_years`, per-model default), correctly under **any** CF calendar.
3. Normalizes every field onto one target grid and pressure-level list
   (from `case.yaml`'s `grid.pressure_levels`), regardless of the model's
   native vertical coordinate, native horizontal grid, or published plev
   list.
4. Applies pyCFRAM's existing missing-variable decision tree (cloud/O₃/
   solar/aerosol/flux/albedo) so any of the 13 candidate DAMIP models —
   most of which are missing at least one of cloud/O₃/`rsdt` — still
   produces a physically valid `(base_state, perturbed_state, nonrad)`
   triple.
5. Writes `provenance.json` recording exactly what happened to each
   variable family, consumed later by `scripts/write_run_summary.py`.

The output feeds the **existing, unmodified** generic writer
(`scripts/build_case_input.py`) and the **existing, unmodified** CFRAM
runner (`scripts/run_parallel_python.py` → `fortran/cfram_rrtmg_1col`).
Phase 3's acceptance criterion (`docs/plan_ph3.md` §11.2) is that `core/`
and `fortran/` are touched **zero** times by this work — verified
mechanically by `tests/test_damip_userguide_example.py::
test_phase3_diff_never_touches_core_or_fortran`, which diffs the whole
Phase-3 branch against its start commit.

## Architecture

```
case.yaml (source.type: cmip6_damip)
    │
    ▼
run_case.py --step build
    │  dispatch: src_type not in {cesm2_cmip6, None} → generic path
    ▼
scripts/build_case_input.py
    │  SOURCE_MODULES['cmip6_damip'] = 'data.cmip6_damip_source'
    │  importlib.import_module(...) triggers @register_source('cmip6_damip')
    │  source = get_source(cfg)          # data/source_base.py factory
    ▼
data/cmip6_damip_source.py :: CMIP6DamipSource.build_states()
    │  10-step flow (see below), using:
    ├──► data/cmip6_common.py            (model-agnostic numerics/detection)
    ├──► configs/damip_models.d/<model>.yaml   (per-model quirks, optional)
    └──► configs/damip_experiments.yaml  (per-experiment semantics, optional)
    │
    │  returns (base_state, perturbed_state, nonrad_forcing) dicts
    │  + writes cases/<case>/input/provenance.json
    ▼
scripts/build_case_input.py :: validate_states() + write_pres_nc/write_surf_nc/write_nonrad_nc
    │  rejects any non-finite value BEFORE writing (no build-after hook possible)
    ▼
cases/<case>/input/{base,perturbed}_{pres,surf}.nc + nonrad_forcing.nc + provenance.json
    │
    ▼
run_case.py --step run
    │  scripts/run_parallel_python.py → fortran/cfram_rrtmg_1col (UNCHANGED)
    ▼
cases/<case>/output/cfram_result.nc
    │
    ▼
scripts/write_run_summary.py   (run_case.py calls this at the end of --step run)
    ▼
cases/<case>/output/<case>.summary.txt
```

Data source for the ESGF acquisition side (upstream of the diagram above,
independent of it — `raw_dir` in `case.yaml` just needs to already be
populated):

```
scripts/download_damip.py  (CLI: --model/--experiment/--dry-run/--all-m4/--all-m5)
    │
    ▼
data/esgf_fetch.py   (pure-stdlib Solr search + HTTP download, no third-party deps)
    │  search_datasets() → list_files() → filename_time_overlap() filter → fetch()
    ▼
raw_data/cmip6_damip/<model>/<experiment>/*.nc + manifest.json
```

## IO contract

### `build_states()` return value

Per the `DataSource` abstract contract (`data/source_base.py`), `build_states()`
returns a 3-tuple:

```python
(base_state, perturbed_state, nonrad_forcing)
```

- `base_state`, `perturbed_state`: `dict[str, np.ndarray]`. 3D variables are
  `(nlev, nlat, nlon)`, **surface→TOA** level order (matching
  `docs/input_spec.md`); 2D variables are `(nlat, nlon)`; plus 1D
  coordinates `lat`, `lon`, `lev` (hPa, surface→TOA). Required keys: the
  full `docs/input_spec.md` variable list (`ta_lay, q, o3, camt, cliq,
  cice, co2, bc, ocphi, ocpho, sulf, ss, dust, ts, ps, solar, albedo`),
  plus the optional `huss` when the model publishes it. Every field is
  guaranteed **finite** (no NaN/Inf) by the time `build_states()` returns —
  `fill_subsurface()` (step 8 below) is what makes that guarantee.
- `nonrad_forcing`: `dict` with `lhflx`, `shflx` — `(nlat, nlon)` W/m²,
  **surface downward-flux convention** (see the sign-law subsection below).

### `provenance.json`

Written by `CMIP6DamipSource._write_provenance()` directly to
`cases/<case>/input/provenance.json` (there is no provenance-writing hook in
`build_case_input.py`, so the source writes it itself). Shape:

```jsonc
{
  "source_type": "cmip6_damip",
  "model": "IPSL-CM6A-LR", "experiment": "hist-aer",
  "variant": "r1i1p1f1", "grid_label": "gr",
  "raw_dir": "raw_data/cmip6_damip",
  "base_years": [1850, 1859], "warm_years": [2011, 2020],
  "target_plev_hpa_toa2sfc": [1, 5, 10, ..., 1000],
  "processes": {
    "albedo": {"status": "ACTIVE", "reason": "rsus/rsds both present -> compute_albedo"},
    "solar":  {"status": "ACTIVE", "reason": "model rsdt present"},
    "flux":   {"status": "ACTIVE", "reason": "hfls/hfss both present; frc_lhflx=-(warm.hfls-base.hfls), ..."},
    "cloud":  {"status": "ACTIVE", "reason": "cl/clw/cli all present and not all-NaN"},
    "o3":     {"status": "MODEL",  "reason": "model o3 (mol/mol) x 48/29 -> kg/kg"},
    "aerosol":{"status": "SKIPPED","reason": "hist-aer Amon has no 3D aerosol mixing ratio; ..."}
  },
  "sanity_checks": {
    "co2_equal": true, "o3_rel_diff": 0.0021, "solar_max_abs_diff_wm2": 0.0
  }
}
```

`scripts/write_run_summary.py` reads this file (degrading gracefully — not
failing — when it's absent, e.g. for non-DAMIP cases) plus
`cfram_result.nc` to produce the human-readable `*.summary.txt`.

## `build_states()`: the 10-step flow

(`data/cmip6_damip_source.py::CMIP6DamipSource.build_states()`, each step
commented in the source with its `docs/plan_ph3.md` §5.1 reference.)

1. **Resolve variant/grid/files.** `variant`/`grid_label` resolution order:
   `case.yaml` explicit → `configs/damip_models.d/<model>.yaml` default →
   runtime glob discovery (`cmip6_common.discover_variant` /
   `_discover_grid_label`, first entry in sorted order, with a warning).
   `cmip6_common.discover_files()` globs the CMOR filename convention
   (`<var>_<table>_<model>_<experiment>_<variant>_<grid>_<trange>.nc`); a
   missing variable is simply absent from the returned dict — not an
   error, since that's the DAMIP main path (§6.1).
2. **Time decoding + month-completeness check.** `cmip6_common.decode_time()`
   (calendar-aware via `cftime`) turns each file's raw `time` values into
   `(year, month, day_weight)`; `_check_month_completeness()` fail-fasts
   with the exact missing year/month if the raw data doesn't fully cover
   the requested `base_years`/`warm_years` window.
3. **Climatology.** `cmip6_common.annual_climo_from_monthly()`: day-weighted
   mean over the selected months, fill-value (`|x| > 1e15`) → NaN first.
4. **Cloud (`cl`/`clw`/`cli`).** All-or-nothing: if all three are present
   and not all-NaN, project hybrid→plev (mass-conserving) or interpolate
   plev→target as needed; otherwise `camt=cliq=cice=0` and the process is
   marked `SKIPPED` (§6.1's decision tree, detailed below).
5. **`ta`/`hus` plev re-interpolation.** `_maybe_interp()` calls
   `cmip6_common.interp_plev_to_target()` — but skips the call entirely
   (identity) when the source plev already equals the target plev exactly,
   which is the common case (target = the standard CMIP6 `plev19` list).
6. **Grid-dependent surface fields.** Albedo from `rsus/rsds`
   (`compute_albedo`); solar from `rsdt` if present, else
   `analytic_solar(lat)` (§6.2 below); `huss` written to the state dict
   only if the model publishes it.
7. **CO₂ + O₃.** CO₂ is a constant 3D field from `case.yaml`'s
   `source.co2.{base,perturbed}_ppmv` (via `DataSource.get_co2()`). O₃
   follows the `o3` mode (`auto`/`use_model`/`climatology`/`skip`, §6.1).
8. **`fill_subsurface()`.** Every below-ground cell (per-column
   `lev > ps`, OR `ta_lay` NaN/fill — CMIP6 masking sometimes extends
   beyond a simple `ps` comparison) is filled: `ta_lay=ts`; `q`/`o3` HOLD
   (copy the lowest real layer down — **not** zero, since RRTMG cannot
   tolerate an exact-zero H₂O/O₃ layer); clouds/aerosol → 0.
   `build_case_input.validate_states()` rejects any surviving non-finite
   value at write time, which is why this step must happen *inside*
   `build_states()` and not as a post-build hook (see "no build-after
   hook" below).
9. **Normalize grid.** `cmip6_common.normalize_grid()`: longitude wrapped
   to ascending `[0, 360)`, latitude sorted ascending S→N, applied
   consistently to every 2D+ field via `np.take`.
10. **Non-radiative forcing + single-forcing sanity checks.**
    `lhflx = -(warm.hfls - base.hfls)`, `shflx = -(warm.hfss - base.hfss)`
    (sign law below); then §5.2's soft sanity checks (CO₂ equality for
    `fixed_1850` experiments, O₃ relative difference, solar max abs diff)
    are computed and recorded in `provenance['sanity_checks']` —
    **warnings only, never a hard gate**, since a real model's data can
    legitimately violate an idealized expectation without being wrong.

### Why there is no "build-after hook"

`scripts/build_case_input.py::validate_states()` runs **before** any NetCDF
is written and rejects non-finite values outright. The pre-Phase-3 CESM2
4×CO₂ pipeline (`data/cesm2_cmip6_source.py` + `scripts/build_cesm2_official.py`)
gets away with a "write placeholder → inject O₃ → mask subsurface" sequence
of three separate scripts precisely because it bypasses `validate_states()`
(it's a bespoke "path B", not the generic writer). DAMIP uses the generic
writer ("path A"), so every fallback — O₃ injection, subsurface fill, unit
conversion — **must** complete before `build_states()` returns.

### Non-radiative forcing sign convention

CMIP6's `hfls`/`hfss` are **positive-upward** (surface losing energy).
CFRAM's non-radiative forcing convention is `frc = Δ(downward surface
energy flux)`, required for the identity
`dT_sfcdyn = dT_ocndyn + dT_lhflx + dT_shflx` to hold. So DAMIP uses
`frc_lhflx = -(warm.hfls - base.hfls)` — the same physical rule as the
pre-existing CESM2 pipeline (`build_cesm2_official.py`, positive-up →
negate), the mirror image of the ERA5 pipeline's rule (`era5_source.py`:
`slhf`/`sshf` are already positive-downward, so `frc = +(event - clim)`,
no sign flip). Both are the *same* physical statement expressed against two
different CF sign conventions.

## Cross-model heterogeneity strategy

`docs/plan_ph3.md` §1.4 surveyed 13 candidate hist-aer models and found only
3 publish all 13 needed Amon variables — **missing cloud/O₃/`rsdt` is the
main path, not an edge case**. Every axis of heterogeneity below is backed
by a real bug hit against real downloaded ESGF data this session; the
generalization exists *because* of the specific bug, not as speculative
future-proofing.

### Calendar

`cmip6_common.decode_time()` replaced `cesm2_cmip6_source.py`'s original
`days/365.0` arithmetic (silently wrong for every calendar except
noleap/365_day) with `cftime.num2date()` + explicit per-calendar day-count
via date arithmetic (`_days_in_month()`). This is required because the M5
model set alone spans five calendars: `noleap` (CESM2), `365_day`
(GISS/CanESM5), `gregorian` (IPSL/MIROC6), `proleptic_gregorian` (MRI), and
`360_day` (HadGEM3 — every month exactly 30 days).

### Vertical coordinate (hybrid sigma-pressure naming)

`cmip6_common.detect_vertical()` reads the `cl` variable's `formula_terms`
CF attribute and dispatches on which coefficient names are present:
`a`+`p0` → CESM-style `p = a·p0 + b·ps` (`hybrid_ab_p0`); `ap` (no `p0`) →
CMOR-mainstream `p = ap + b·ps` (`hybrid_ap_b`, IPSL/MRI/CNRM/HadGEM).
`normalize_hybrid_coeffs()` converts either convention to a single
`(a_eff, b, p0_eff)` triple that `hybrid_to_plev_mass_conserving()` — moved
*unmodified* from `cesm2_cmip6_source.py` since it was already validated
bit-exact — consumes directly.

**Real bug #4 (IPSL) — three separate problems in one coordinate:**
1. IPSL's `cl`/`clw`/`cli` files have **no `formula_terms` attribute at
   all** linking `ap`/`b` to the vertical coordinate — auto-detection has
   no string to probe, so `configs/damip_models.d/IPSL-CM6A-LR.yaml`
   supplies an explicit `vertical: {scheme: hybrid_ap_b, ap: ap, b: b}`
   override, which `detect_vertical(override=...)` returns verbatim,
   bypassing probing entirely.
2. IPSL publishes `ap`/`b` at **layer interfaces** (`nlev+1` values, on a
   dim like `klevp1`), not at layer **midpoints** (`nlev` values, matching
   the field data's own level dimension). `_read_hybrid_coeffs()` detects
   this via `nlev_data` and pairwise-averages adjacent interface values —
   this is *exact*, not approximate, because layer pressure
   `p = ap·p0 + b·ps` is linear in `ap`/`b` for fixed `p0`/`ps`, so
   averaging two interface pressures gives exactly the layer's midpoint
   pressure.
3. IPSL stores hybrid levels **surface→TOA** (`b` descending from ~1 to
   ~0) — the opposite of CESM2's native `hyam`/`hybm` (TOA→surface,
   `hybrid_to_plev_mass_conserving`'s only previously-tested convention).
   `_read_hybrid_coeffs()` auto-detects the direction from `b`'s monotonic
   sign (`flipped = b_eff[0] > b_eff[-1]`) and reverses the coefficient
   arrays; the caller must reverse the **field** data (`cl`/`clw`/`cli`)
   the exact same way — reversing only one side silently scrambles the
   profile into a plausible-looking but physically wrong near-zero cloud
   field (this was hit and diagnosed by sanity-checking against the real
   raw data's actual cloud-fraction magnitude, not by a crash).

**Real bug #5 (CNRM-CM6-1) — non-standard attribute name:** CNRM's `lev`
coordinate carries `formula_term` (**singular**), not the standard CF
`formula_terms` (plural) — a CERFACS CMOR-configuration quirk. Auto-probing
never finds it and silently falls through to the wrong ("already on plev")
branch. Fixed the same way as IPSL: an explicit `vertical:` override in
`configs/damip_models.d/CNRM-CM6-1.yaml`.

**Real bug #7 (HadGEM3-GC31-LL) — genuinely unsupported coordinate type:**
HadGEM3's cloud vertical coordinate is `formula_terms = "a: lev b: b
orog: orog"` — **hybrid-HEIGHT** (`z = a + b·orog`), not hybrid-sigma-
**pressure**. `detect_vertical()` correctly raises `ValueError` for a
`formula_terms` string it doesn't recognize (neither `p0`+`a` nor bare
`ap`). `cmip6_damip_source.py` catches exactly this `ValueError` around the
hybrid-projection call and **degrades to `cloud=SKIPPED`** rather than
crashing the whole build — converting height to pressure genuinely needs
the model's own temperature profile, a different physical conversion, not
a coefficient-naming variant, so this is treated as an anticipated
known-issue (`docs/plan_ph3.md` §12 R4), not a bug to force through.

### Published pressure levels

`cmip6_common.interp_plev_to_target()` — log-pressure linear interpolation
along axis 0, with an identity guarantee when target plev == source plev
exactly (the common case: target = standard CMIP6 `plev19`, verified by a
dedicated unit test). Used for `ta`/`hus`/`o3` (and `cl`/`clw`/`cli` when
already plev-based rather than hybrid) whenever a model's native published
plev list differs from the case's target grid.

### Horizontal grid

**Real bug #6 (MRI-ESM2-0) — a variable on a different native grid than
its own siblings:** MRI's `o3` is published on a coarser native grid
(64×128) than every other Amon variable in the same model/table (160×320)
— CMIP6 imposes no requirement that every variable in `Amon` share one
grid. `cmip6_damip_source.py` detects this by comparing `o3`'s own
`(lat, lon)` against the reference grid (established from `ta`) and, if
different, regrids via `cmip6_common.regrid_horizontal_bilinear()`
(bilinear, periodic in longitude, implemented independently of
`core/kernels.py`'s near-identical helper since `data/` and `core/` are
kept as separate layers per the zero-core-changes rule).

That regrid needed a **second** fix: CMIP6's standard below-ground-plev
masking (NaN at the surface end of the level axis over high terrain, e.g.
Tibetan Plateau) must be filled *before* bilinear interpolation, or a NaN
input cell propagates into its output *neighbors* too — corrupting more
cells than were originally missing. `cmip6_common.
fill_nan_hold_toward_surface()` runs first (hold the shallowest valid
level's value down through the masked run), then the regrid, then the
*real*, ps-aware `fill_subsurface()` (step 8) runs again downstream using
the **target** grid's own `ps` for the final, physically precise fill —
`fill_nan_hold_toward_surface` is only a NaN-propagation guard for the
intermediate regridding step, not the final subsurface treatment.

### Units

- `cl`: CMIP6 publishes percent (0–100) → pyCFRAM's `camt` is a 0–1
  fraction (`÷100`).
- `o3`: CMIP6 publishes mol/mol (volume mixing ratio) → pyCFRAM's `o3` is
  kg/kg (mass mixing ratio), conversion factor `VMR_TO_MMR = 48/29`
  (`M_O3`/`M_air`), the same constant used by the pre-existing
  `scripts/inject_cesm_o3.py`. **A units bug here is a ~5-order-of-
  magnitude error**, not a small one — this is exactly the class of bug
  that produced the "huss-as-O3" data-corruption incident documented in
  `session_log.md` (2026-05-10/12), where a collaborator's Fortran CFRAM
  run fed specific humidity into the O₃ channel by mistake. `data/
  cmip6_damip_source.py`'s O₃ branch is written specifically to make that
  class of error structurally hard to reintroduce (conversion applied
  exactly once, in one place, immediately after reading the raw mol/mol
  field).

## ESGF client (`data/esgf_fetch.py` + `scripts/download_damip.py`)

Pure-stdlib (`urllib` only, no `requests`/`intake-esgf`/`esgpull`) Solr
search + HTTP download client. Three real bugs shaped its retry/parsing
logic:

**Real bug #1 — unfiltered variant search inflates download volume.**
IPSL-CM6A-LR alone has 10 published `hist-aer` variants
(`r1i1p1f1`...`r10i1p1f1`); an unfiltered `--model IPSL-CM6A-LR` dry-run
matched all of them, reporting 54.8 GB instead of the correct ~20.4 GB for
one variant. Fixed via `download_damip.py`'s `MODEL_DEFAULTS` dict — a
per-model `(variant_label, grid_label)` default applied to every search —
documented in the code as "not an optional nicety — omitting it is a
correctness bug, not just a size optimization."

**Real bug #2 — cross-node replica duplication.** The same physical file
is commonly indexed at more than one ESGF data node (e.g. a CNRM file
hosted at its home institution `esg1.umr-cnrm.fr` *and* mirrored at
`esgf3.dkrz.de` with a fully independent URL). The obvious fix — filter on
Solr's `replica=false` flag — **does not work**: some models (confirmed:
MRI-ESM2-0) are indexed as `replica=true` at *every* node this client
queries, because DKRZ/CEDA aren't that model's home data node, so
`replica=false` returns **zero results even though the data is genuinely
downloadable**. Fixed instead by client-side deduplication keyed on
filename (`download_damip.search_damip_files()`), which **also** collects
every distinct URL seen for that filename into a list rather than keeping
only the first — enabling real bug #8 below.

**Real bug #3 — `checksum`/`checksum_type` are multi-valued Solr fields.**
ESGF's Solr schema returns these as **lists** (e.g.
`checksum_type=['SHA256']`), not scalars, for files that carry more than
one checksum type. This only surfaces once a real, complete file download
reaches the checksum-verification step — it does not show up in a
dry-run, which never reaches that code path. `esgf_fetch.list_files()` now
unwraps to a scalar (`checksum[0] if checksum else ''`) before returning.

**Real bug #8 (CNRM-CM6-1) — home node timeout, working mirror exists.**
CNRM's home data node (`esg1.umr-cnrm.fr`) timed out on 7 of 17 `cl`
file-time-chunks (~26 GB of CNRM's total ~35 GB) on first attempt; a fully
working mirror exists at `esgf3.dkrz.de`. `download_files()` now retries
every known mirror URL for a file (the list built by bug #2's fix) in
order before giving up, rather than failing on the first URL's error.

**Real bug #9 — `http.client.RemoteDisconnected` isn't a `URLError`
subclass.** Hit as an uncaught crash during `list_files()`'s `urlopen()`
call while dry-running NorESM2-LM (WP-M5.3 prep).
`RemoteDisconnected` is raised directly by `http.client` during
`response.begin()` — a `ConnectionResetError`/`OSError` that `urlopen()`
doesn't always wrap the way it wraps connection-time failures. Every
`except urllib.error.URLError` clause in `esgf_fetch.py` (`search_datasets`
×2, `list_files`, `fetch`) now also catches `OSError`, so a transient
mid-response disconnect surfaces as the same clean `RuntimeError` as any
other network hiccup instead of crashing the process (commit `15a2845`).

**Node fallback:** `search_with_fallback()` tries DKRZ then CEDA (LLNL's
classic Solr endpoint is dead, redirects to a bridge — not used). Both are
usable without any ESGF login for file download (Range-request `206
Partial Content` against THREDDS fileServer works unauthenticated); the
contract's "may require an ESGF account" clause turned out to be a
precaution that wasn't needed in practice.

## Missing-variable decision tree (§6.1)

All of the following happen *inside* `build_states()`, never as a
build-after hook (see above). Each row's status is one of the vocabulary
values recorded in `provenance.json['processes'][<name>]['status']`.

| Variable family | Condition | Status | Effect |
|---|---|---|---|
| Cloud (`cl`/`clw`/`cli`) | all 3 present, not all-NaN, vertical coordinate recognized | `ACTIVE` | hybrid→plev or plev→target projected normally |
| | any of the 3 missing, or all-NaN, or `detect_vertical` raises (unsupported coordinate, e.g. hybrid-height) | `SKIPPED` | `camt=cliq=cice=0`; **no partial cloud** — 2-of-3 present still counts as missing (GISS has `clw`/`cli` but not `cl`: still `SKIPPED`) |
| O₃ | `o3=use_model`, or `o3=auto` and model publishes `o3` | `MODEL` | model's own `o3` (mol/mol) × `48/29` → kg/kg |
| | `o3=auto`/`climatology` and no model `o3`, climatology source file found | `CLIMATOLOGY` | inject CESM 1850 climatology (`cmip6_common.o3_climatology()`), same field both states → `frc_o3 ≡ 0` |
| | `o3=skip`, or climatology source missing | `SKIPPED` | constant `1e-12` kg/kg (not exact zero — RRTMG cannot tolerate an exact-zero absorber layer) |
| Solar (`rsdt`) | model publishes `rsdt` | `ACTIVE` | climatological `rsdt` used directly |
| | `rsdt` not published | `ANALYTIC` | `cmip6_common.analytic_solar(lat)` fallback (§6.2 below); same field both states → `frc_solar ≡ 0` for hist-aer |
| Aerosol | always, for `hist-aer` (Amon has no 3D aerosol mixing ratio) | `SKIPPED` | `aerosol.source: zero` → all-zero 6-species arrays → runner's own `aer_max < 1e-15` auto-detection skips the aerosol decomposition entirely |
| Flux (`hfls`/`hfss`) | both present (13/13 models, in practice) | `ACTIVE` | `nonrad` dict populated, sign law above |
| Albedo (`rsus`/`rsds`) | both present (13/13 models, in practice) | `ACTIVE` | `compute_albedo()` |
| `huss` | present | (written to surf state) | optional; absent → simply not written, runner's own HOLD fallback applies |

### `analytic_solar(lat)` (§6.2)

Standard astronomical daily-mean-insolation formula (solar declination →
hour angle → daily-mean flux, averaged over 365 days), used whenever a
model doesn't publish `rsdt` (5 of 13 candidate models, including CESM2 and
CanESM5). Since `hist-aer` freezes the solar constant, base and perturbed
states get the **same** field from this function by construction, so
`frc_solar == 0` regardless of any small absolute bias in the fallback's
value — the fallback only needs to be a physically reasonable *background*
field, not an exact one. Verified against two anchors (area-weighted global
mean ≈ S₀/4 ≈ 340 W/m², equatorial annual mean ≈ 417 W/m²) and, for models
that *do* publish `rsdt`, cross-checked against the real zonal-mean curve
(near-coincident, by design of the same unit test).

## `configs/damip_models.d/` + `configs/damip_experiments.yaml`

**Design principle**: `cmip6_common.py` holds only model-agnostic
numerics/detection; every model-specific quirk (default variant/grid,
calendar, explicit vertical-coordinate override, non-default warm-year
window) lives in a small declarative YAML file, one per model, loaded by
`load_model_config()` (glob + `yaml.safe_load`, `{}` with a printed warning
if the file is absent — never a hard failure, since a model can still work
via pure runtime auto-detection). **This is the M5 "adding a model = adding
a file" contract** — see `docs/m5_multimodel_userguide.md`'s worked
NorESM2-LM example, which adds a 9th model with *zero* Python changes.

Representative entries (abbreviated; see the files themselves for full
comments):

```yaml
# configs/damip_models.d/CESM2.yaml — "skip full house" demonstrator
default_variant: r1i1p1f1
default_grid: gn
calendar: noleap
warm_years_default: [2005, 2014]      # hist-aer data ends 2014-12
vertical: {scheme: hybrid_ab_p0, a: a, b: b, p0: p0}
missing_ok: [cl, clw, cli, rsdt, o3]  # declarative: these gaps are expected

# configs/damip_models.d/IPSL-CM6A-LR.yaml — explicit override required
vertical: {scheme: hybrid_ap_b, ap: ap, b: b}   # no formula_terms attribute at all

# configs/damip_models.d/CNRM-CM6-1.yaml — non-standard attribute name
vertical: {scheme: hybrid_ap_b, ap: ap, b: b}   # attribute is `formula_term` (singular)
```

`configs/damip_experiments.yaml` records per-`experiment_id` semantics
(which forcing is time-varying, CO₂/O₃ default mode) used **only** for the
§5.2 sanity-check warnings — never as the source of truth for CO₂/O₃
values (those always come from `case.yaml`'s own `source:` block).

## `case.yaml` contract

```yaml
source:
  type: cmip6_damip
  model: IPSL-CM6A-LR
  experiment: hist-aer
  variant: r1i1p1f1            # optional — falls back to models.d default, then glob
  grid_label: gr                # optional — same fallback chain
  raw_dir: raw_data/cmip6_damip
  base_years: [1850, 1859]
  warm_years: [2011, 2020]      # per-model default; CESM2/GISS use [2005, 2014]
  co2: {source: constant, base_ppmv: 284.7, perturbed_ppmv: 284.7}
  o3: auto                      # auto | use_model | climatology | skip
  aerosol: {source: zero}

grid:
  pressure_levels: [1, 5, 10, ..., 1000]   # TOA->sfc, hPa

input:                          # REQUIRED — core/config.py's load_case parses paths from
  base_pres: input/base_pres.nc            # this block; a case.yaml missing it fails with
  base_surf: input/base_surf.nc            # KeyError at run step, not build step.
  perturbed_pres: input/perturbed_pres.nc
  perturbed_surf: input/perturbed_surf.nc
  nonrad_forcing: input/nonrad_forcing.nc

radiation:
  scheme: rrtmg
run:
  nproc: auto
```

## `*.summary.txt` (`scripts/write_run_summary.py`)

Standalone post-processing script — **not** imported by, and not coupled
to, `run_parallel_python.py` (already-accepted Phase-2 code; only its
*output* NetCDF is read). `run_case.py` invokes it via the same
`subprocess.call`-based `run_step()` helper used for build/plot, at the end
of the `run` step only. Degrades gracefully (omits the process-activity
table with an explanatory note) when `provenance.json` is absent — e.g.
ERA5, `cesm2_cmip6`-collab, or idealized climlab cases, which never
produced one.

Sections, in order:

1. **Header** — case/model/experiment/variant/grid/base_years/warm_years/
   target plev/grid dimensions/radiation engine.
2. **Process activity table** — one line per family (`albedo`, `solar`,
   `flux`, `cloud`, `o3`, `aerosol`) with its recorded status and reason,
   read straight from `provenance.json`.
3. **Additivity residual** — `dT_observed[sfc] − Σ dT_X[sfc]`, domain-mean
   and domain-max. The term set summed is computed by
   `select_additive_terms()`, which excludes known **exact or
   near-exact aggregates** to avoid double-counting: `dry` (=
   `atmdyn+sfcdyn`), `sfcdyn` (= `ocndyn+lhflx+shflx`, exact to machine
   precision per `session_log.md` 2026-06-01), `cloud` (=
   `cloud_lw+cloud_sw`, exact), and the per-species aerosol breakdown
   (≈ bulk `dT_aerosol` with a real ~2.5–4% non-linear residual — kept
   separate from the bulk sum rather than double-counted). Explicitly
   cites `docs/m3_methodology_comparison.md` for why a residual here is
   expected first-order-CFRAM non-linearity, not a bug.
4. **Acceptance gates** — the M4 end-to-end gate from `docs/plan_ph3.md`
   §8, recomputed from `cfram_result.nc` on every run rather than checked
   once by hand:
   - `sfcdyn_identity_max_abs` — max over **all levels** of
     `|dT_ocndyn + dT_lhflx + dT_shflx − dT_sfcdyn|`, gated at 1e-10 K.
     This identity holds by construction, so a non-trivial value means the
     decomposition wiring regressed; it is the cheapest tripwire in the
     pipeline, and is where a real double-counting bug first surfaced
     (`session_log.md` 2026-06-01). Measured 3.6e-15 – 1.6e-14 K across the
     nine M4/M5 cases.
   - `dT_observed[sfc] NaN frac` — gated at 0. Below-ground *atmospheric*
     levels are legitimately NaN (`plev > local ps`, ~40% of cells at
     1000 hPa), so only the surface row is gated.
   - `dT_observed[sfc] gmean` and `NH vs SH mean` — area-weighted, judged
     only for experiments carrying `net_cooling: true` in
     `configs/damip_experiments.yaml` (hist-aer expects net cooling in
     (−3, 0) K with the NH cooler, since anthropogenic aerosol is
     NH-concentrated). For other experiments the same numbers are printed
     as INFO rather than judged against the wrong expectation — hist-GHG
     warms.
5. **Single-forcing consistency** — the §5.2 sanity-check values
   (`co2_equal`, `o3_rel_diff`, `solar_max_abs_diff_wm2`) with a PASS/
   WARNING judgment looked up against `configs/damip_experiments.yaml`
   (so, e.g., `hist-nat` legitimately showing a large solar diff is not
   flagged, since solar *is* that experiment's own varying forcing).
6. **Footer** — reproduction command + md5 of every input NetCDF.

## Contract compliance

`data/cmip6_common.py`'s functions directly implement the contract's
DAMIP-specific requirements: experiment_id recognition (via `case.yaml` +
`configs/damip_experiments.yaml`), automatic base↔single-forcing pairing
(first vs last decade climatology, `_climo_pair_for_variable()`), and
"dynamic detection of each model's available pressure levels with
Python-side interpolation to a unified grid" (`interp_plev_to_target()`,
identity-preserving when source==target).

## Testing

```bash
coverage run -m pytest tests/ -q
coverage report -m --include="*/data/cmip6_common.py,*/data/cmip6_damip_source.py,*/data/esgf_fetch.py"
```

| Module | Statements | Coverage |
|---|---:|---:|
| `data/cmip6_common.py` | 272 | **98%** |
| `data/cmip6_damip_source.py` | 384 | **92%** |
| `data/esgf_fetch.py` | 171 | **99%** |
| **Total** | 827 | **96%** |

All three are well above the contract's 60% hard gate (re-measured
2026-07-06 after the regrid/fill_nan_hold/hybrid-height/multi-mirror/
OSError additions described above — no coverage regression from WP-M4.6's
original 98%/93%/98% baseline). Test files: `tests/test_cmip6_common.py`
(calendar × 3, hybrid mass-conservation, plev-interpolation identity,
grid-normalization idempotency, `analytic_solar` anchors),
`tests/test_damip_source.py` (all four decision-tree branches — full-
variable / missing-cloud / missing-rsdt / missing-o3 — on synthetic
`tests/data/damip_smoke/` fixtures spanning different calendars and
hybrid-coefficient naming conventions), `tests/test_esgf_fetch.py`
(monkeypatched Solr JSON + HTTP responses, no real network calls),
`tests/test_download_damip.py`, `tests/test_write_run_summary.py`,
`tests/test_damip_regression.py` (CESM2 4×CO₂ build-output golden
comparison, guarding the `cesm2_cmip6_source.py` refactor into a thin
shim over `cmip6_common.py`), and `tests/test_damip_userguide_example.py`
(the `core/`/`fortran/` zero-touch acceptance test).

## See also

- `docs/plan_ph3.md` — full execution plan (§1–§13), including the ESGF
  data-availability matrix behind the M4/M5 model selection.
- `docs/m5_multimodel_userguide.md` — user-facing guide: per-model support
  matrix, known-issues log, and the worked "add a new model" example.
- `docs/input_spec.md` — the NetCDF format contract this module's output
  must satisfy.
- `docs/m3_methodology_comparison.md` — why CFRAM additivity residuals are
  expected, referenced by the summary's residual section.
