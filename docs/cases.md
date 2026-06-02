# Contract M1 Case Definitions

Two CFRAM attribution case studies under deliverable **M1**. Status as of
2026-06-01:

| Case | Status | Artifacts |
|------|--------|-----------|
| ① Alaska wildfire 2022 (aerosol) | **Deferred** | definition below; no case dir yet |
| ② India–Bangladesh wet-heat 2023 | **In progress** | `cases/india_wb23/case.yaml`; ERA5 download pipeline built (`scripts/download_era5_india23{,_ncar}.py`) |

---

## ① Alaska wildfire, Jun–Jul 2022 (aerosol-driven event)

**Time**
- Main analysis period: 2022-06-10 – 2022-07-10
- Background: 1991–2020 same-period climatology
- Auxiliary windows: 2022-05-20 – 06-09, 2022-07-11 – 07-31

**Space**
- Main fire region: 58°–66°N, 170°–150°W
- Extended region: 60°–72°N, 175°–135°W
- Arctic response region: 65°–75°N, 170°–120°W

**Target diagnostics**
- 2 m air-temperature anomaly (temperature response)
- AOD / black-carbon AOD anomaly (smoke-plume / aerosol signal)
- Surface shortwave-radiation and albedo anomalies
- Spatial distribution of each CFRAM term
- Regional-mean bar charts for the three regions above

## ② India (east) – Bangladesh extreme wet-heat, Apr 2023

**Time**: April 2023
- Apr 1–16: pre-event background build-up
- Apr 17–20: extreme wet-heat core period
- Apr 21–30: post-event decay / adjustment

**Space**: 10°–35°N, 65°–110°E
- Main region: 20°–27°N, 85°–93°E
- Extended region: 15°–30°N, 75°–100°E
- Control (dry-heat) region: 24°–32°N, 68°–78°E

**Target diagnostics**: wet-heat indices anomaly + moisture transport +
CFRAM decomposition + wet-heat vs dry-heat regional contrast.

**Implementation notes** (`cases/india_wb23/case.yaml`)
- warm core days `warm_days=[16,17,18,19]` (Apr 17–20)
- CO2 base 385 / perturbed 421.3 ppmv; aerosol from MERRA-2
- ERA5 6-hourly April download via `scripts/download_era5_india23.py`
  (CDS, 37-level) or `download_era5_india23_ncar.py` (NCAR RDA mirror).
  Per-request cost cap forces 2-variable pairs; see `session_log.md`
  2026-05-31 for the download topology.
