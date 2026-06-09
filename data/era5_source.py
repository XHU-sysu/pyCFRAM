"""ERA5 daily data source for pyCFRAM input generation.

Implements the Wu et al. (2025) methodology:
  1. Load 6-hourly ERA5 data for August of each year (2003-2022)
  2. Identify warm period via regional Tmax > 90th percentile
  3. Compute base state = climatological warm-period mean
  4. Compute perturbed state = event-year warm-period mean
  5. Output standard pyCFRAM state dicts

Data layout on disk (from CDS API download):
  era5_data/daily/
    era5_pl_{var}_{YYYY}08.nc          # PL: t, q, o3, cc, ciwc, clwc
    era5_sl_{YYYY}08/
      data_stream-oper_stepType-instant.nc   # skt, sp
      data_stream-oper_stepType-accum.nc     # ssrd, ssr, tisr
      data_stream-oper_stepType-max.nc       # mx2t
"""

import os
import calendar
from datetime import date, timedelta
import numpy as np
from netCDF4 import Dataset

from .source_base import DataSource, register_source

# ERA5 PL variable short names → CFRAM variable names
PL_VAR_MAP = {
    't': 'ta_lay',
    'q': 'q',
    'o3': 'o3',
    'cc': 'camt',
    'clwc': 'cliq',
    'ciwc': 'cice',
}

# Number of 6-hourly steps per day
STEPS_PER_DAY = 4
NDAYS_AUG = 31
SECS_PER_DAY = 86400.0


def _detect_data_var(nc, skip=None):
    """Find the main data variable in a NetCDF file (skip coord vars)."""
    if skip is None:
        skip = {'time', 'valid_time', 'latitude', 'longitude', 'lat', 'lon',
                'level', 'pressure_level', 'number', 'expver'}
    for v in nc.variables:
        if v.lower() not in skip and nc.variables[v].ndim >= 3:
            return v
    return None


def _load_pl(data_dir, var_short, year, month=8):
    """Load one PL variable for one month.

    Returns:
        data: (ntime, nlev, nlat, nlon) float64
        lev, lat, lon: 1D arrays
    """
    fname = os.path.join(data_dir, f'era5_pl_{var_short}_{year}{month:02d}.nc')
    nc = Dataset(fname, 'r')
    vname = _detect_data_var(nc)
    data = np.array(nc.variables[vname][:], dtype=np.float64)
    lev = np.array(nc.variables['pressure_level'][:], dtype=np.float64)
    lat = np.array(nc.variables['latitude'][:], dtype=np.float64)
    lon = np.array(nc.variables['longitude'][:], dtype=np.float64)
    nc.close()
    return data, lev, lat, lon


def _load_sl(data_dir, year, month=8):
    """Load all SL variables for one month from subdirectory structure.

    Returns dict with keys: skt, sp, ssrd, ssr, tisr, mx2t
    Each value is (ntime, nlat, nlon) float64.
    Also returns lat, lon.
    """
    sl_dir = os.path.join(data_dir, f'era5_sl_{year}{month:02d}')
    result = {}

    # Instant variables: skt, sp
    nc = Dataset(os.path.join(sl_dir, 'data_stream-oper_stepType-instant.nc'))
    for v in ('skt', 'sp'):
        if v in nc.variables:
            result[v] = np.array(nc.variables[v][:], dtype=np.float64)
    lat = np.array(nc.variables['latitude'][:], dtype=np.float64)
    lon = np.array(nc.variables['longitude'][:], dtype=np.float64)
    nc.close()

    # Accumulated variables: ssrd, ssr, tisr
    nc = Dataset(os.path.join(sl_dir, 'data_stream-oper_stepType-accum.nc'))
    for v in ('ssrd', 'ssr', 'tisr'):
        if v in nc.variables:
            result[v] = np.array(nc.variables[v][:], dtype=np.float64)
    nc.close()

    # Surface heat fluxes: slhf, sshf (separate accum file)
    fflux = os.path.join(sl_dir, 'data_stream-oper_stepType-accum_flux.nc')
    if os.path.exists(fflux):
        nc = Dataset(fflux)
        for v in ('slhf', 'sshf'):
            if v in nc.variables:
                result[v] = np.array(nc.variables[v][:], dtype=np.float64)
        nc.close()

    # Max variables: mx2t
    fmax = os.path.join(sl_dir, 'data_stream-oper_stepType-max.nc')
    if os.path.exists(fmax):
        nc = Dataset(fmax)
        vname = _detect_data_var(nc)
        if vname:
            result['mx2t'] = np.array(nc.variables[vname][:], dtype=np.float64)
        nc.close()

    return result, lat, lon


def _sixhourly_to_daily_mean(data):
    """Convert 6-hourly data to daily means.

    Args:
        data: (..., ntime, ...) where ntime = NDAYS * STEPS_PER_DAY

    First axis is time. Returns (ndays, ...) array.
    """
    ntime = data.shape[0]
    ndays = ntime // STEPS_PER_DAY
    rest = data.shape[1:]
    reshaped = data[:ndays * STEPS_PER_DAY].reshape(ndays, STEPS_PER_DAY, *rest)
    return reshaped.mean(axis=1)


def _sixhourly_accum_to_daily_wm2(data):
    """Convert 6-hourly step accumulations (J/m²) to daily mean W/m².

    ERA5 CDS convention: each 6-hourly output step contains the accumulation
    for that SINGLE HOUR (forecast step), not the full 6-hour window.
    To extrapolate to a proper 6-hour accumulation, multiply each step by 6.
    Daily total = 6 × sum of 4 steps. Daily mean W/m² = total / 86400.

    Equivalently: daily_mean = sum(steps) / (86400 / 6) = sum(steps) / 14400.
    """
    ntime = data.shape[0]
    ndays = ntime // STEPS_PER_DAY
    rest = data.shape[1:]
    reshaped = data[:ndays * STEPS_PER_DAY].reshape(ndays, STEPS_PER_DAY, *rest)
    daily_total = reshaped.sum(axis=1) * 6  # Scale hourly accum → 6-hourly
    return daily_total / SECS_PER_DAY       # W/m²


def _sixhourly_to_daily_max(data):
    """Convert 6-hourly data to daily maximum."""
    ntime = data.shape[0]
    ndays = ntime // STEPS_PER_DAY
    rest = data.shape[1:]
    reshaped = data[:ndays * STEPS_PER_DAY].reshape(ndays, STEPS_PER_DAY, *rest)
    return reshaped.max(axis=1)


def identify_warm_period(data_dir, years, event_year, month,
                         lat, lon, detect_region,
                         threshold_pct=90, min_consecutive=3):
    """Identify warm period days using Wu et al. methodology.

    Args:
        data_dir: path to daily ERA5 data
        years: list of climatology years
        event_year: the extreme event year
        month: month number (8 for August)
        lat, lon: coordinate arrays from ERA5
        detect_region: dict with 'lat' and 'lon' ranges
        threshold_pct: percentile threshold for Tmax
        min_consecutive: minimum consecutive days

    Returns:
        warm_days: list of 0-based day indices within the month
    """
    lat_range = detect_region['lat']
    lon_range = detect_region['lon']
    lat_mask = (lat >= min(lat_range)) & (lat <= max(lat_range))
    lon_mask = (lon >= min(lon_range)) & (lon <= max(lon_range))

    # Load Tmax for all climatology years
    all_tmax = []
    for yr in years:
        try:
            sl, _, _ = _load_sl(data_dir, yr, month)
            if 'mx2t' not in sl:
                continue
            mx2t_daily = _sixhourly_to_daily_max(sl['mx2t'])[:NDAYS_AUG]
            # Regional mean
            regional = mx2t_daily[:, lat_mask, :][:, :, lon_mask]
            all_tmax.append(regional.mean(axis=(1, 2)))
        except Exception as e:
            print(f"  Warning: skip mx2t {yr}: {e}")

    if not all_tmax:
        raise RuntimeError("No mx2t data loaded for warm period detection")

    tmax_stack = np.array(all_tmax)  # (nyears, 31)
    threshold = np.percentile(tmax_stack, threshold_pct, axis=0)

    # Event year Tmax
    sl_event, _, _ = _load_sl(data_dir, event_year, month)
    mx2t_event = _sixhourly_to_daily_max(sl_event['mx2t'])[:NDAYS_AUG]
    event_regional = mx2t_event[:, lat_mask, :][:, :, lon_mask]
    event_tmax = event_regional.mean(axis=(1, 2))

    # Find consecutive days exceeding threshold
    exceed = event_tmax > threshold
    warm_days = []
    consecutive = 0
    start = -1
    for d in range(NDAYS_AUG):
        if exceed[d]:
            if consecutive == 0:
                start = d
            consecutive += 1
        else:
            if consecutive >= min_consecutive:
                warm_days.extend(range(start, start + consecutive))
            consecutive = 0
            start = -1
    if consecutive >= min_consecutive:
        warm_days.extend(range(start, start + consecutive))

    return sorted(set(warm_days))


@register_source('era5_daily')
class ERA5DailySource(DataSource):
    """Build CFRAM states from ERA5 6-hourly daily data (Wu et al. method)."""

    def build_states(self):
        tcfg = self.source_cfg['temporal']
        event_year = tcfg['event_year']
        month = tcfg.get('event_month', 8)
        clim_start, clim_end = tcfg['clim_years']
        clim_years = list(range(clim_start, clim_end + 1))

        # Resolve data directory
        data_dir = self.source_cfg['data_dir']
        if not os.path.isabs(data_dir):
            from core.config import PROJECT_ROOT
            data_dir = os.path.join(PROJECT_ROOT, data_dir)

        # ── Step 1: Identify warm period ─────────────────────────────────
        warm_detect = tcfg.get('warm_detect', {})
        manual_days = tcfg.get('warm_days', None)

        if manual_days is not None:
            warm_days = list(manual_days)
            print(f"Using manual warm days: {warm_days}")
        else:
            # Need lat/lon from a PL file to compute regional mask
            _, _, lat_raw, lon_raw = _load_pl(data_dir, 't', event_year, month)
            detect_region = warm_detect.get('detect_region',
                                            self.cfg.get('plot', {}).get('key_region', {}))
            warm_days = identify_warm_period(
                data_dir, clim_years, event_year, month,
                lat_raw, lon_raw, detect_region,
                threshold_pct=warm_detect.get('threshold_pct', 90),
                min_consecutive=warm_detect.get('min_consecutive', 3),
            )

        if not warm_days:
            raise RuntimeError("No warm period days identified. "
                               "Check detect_region and threshold settings.")

        print(f"Warm period: Aug {warm_days[0]+1}–{warm_days[-1]+1} "
              f"({len(warm_days)} days)")

        # ── Step 2: Load PL variables ────────────────────────────────────
        # First pass: get coordinates
        _, lev, lat_raw, lon_raw = _load_pl(data_dir, 't', event_year, month)

        # Flip latitude to ascending (S→N) if needed
        lat_ascending = lat_raw[0] < lat_raw[-1]
        if not lat_ascending:
            lat = lat_raw[::-1]
        else:
            lat = lat_raw
        lon = lon_raw

        nlev, nlat, nlon = len(lev), len(lat), len(lon)
        print(f"Grid: {nlev} levels, {nlat} lat, {nlon} lon")

        base_state = {'lat': lat, 'lon': lon, 'lev': lev}
        pert_state = {'lat': lat, 'lon': lon, 'lev': lev}

        for era5_var, cfram_var in PL_VAR_MAP.items():
            print(f"  Processing PL: {era5_var} → {cfram_var}")
            # Climatological mean over warm days
            clim_sum = np.zeros((nlev, nlat, nlon), dtype=np.float64)
            clim_count = 0
            for yr in clim_years:
                try:
                    data, _, _, _ = _load_pl(data_dir, era5_var, yr, month)
                    daily = _sixhourly_to_daily_mean(data)[:NDAYS_AUG]
                    warm_mean = daily[warm_days].mean(axis=0)
                    if not lat_ascending:
                        warm_mean = warm_mean[:, ::-1, :]
                    clim_sum += warm_mean
                    clim_count += 1
                except Exception as e:
                    print(f"    Warning: skip {era5_var} {yr}: {e}")
            base_state[cfram_var] = clim_sum / clim_count

            # Event year warm period mean
            data_ev, _, _, _ = _load_pl(data_dir, era5_var, event_year, month)
            daily_ev = _sixhourly_to_daily_mean(data_ev)[:NDAYS_AUG]
            pert_mean = daily_ev[warm_days].mean(axis=0)
            if not lat_ascending:
                pert_mean = pert_mean[:, ::-1, :]
            pert_state[cfram_var] = pert_mean

        # ── Step 3: Load SL variables ────────────────────────────────────
        print("  Processing SL variables...")

        # Climatological SL mean
        sl_keys = ['skt', 'sp', 'ssrd', 'ssr', 'tisr', 'slhf', 'sshf']
        sl_clim = {k: np.zeros((nlat, nlon), dtype=np.float64) for k in sl_keys}
        sl_count = 0
        for yr in clim_years:
            try:
                sl, _, _ = _load_sl(data_dir, yr, month)
                # Instant vars: 6-hourly → daily mean
                for v in ('skt', 'sp'):
                    daily = _sixhourly_to_daily_mean(sl[v])[:NDAYS_AUG]
                    warm_mean = daily[warm_days].mean(axis=0)
                    if not lat_ascending:
                        warm_mean = warm_mean[::-1, :]
                    sl_clim[v] += warm_mean
                # Accum vars: 6-hourly step accum → daily W/m²
                for v in ('ssrd', 'ssr', 'tisr', 'slhf', 'sshf'):
                    if v not in sl:
                        continue
                    daily_wm2 = _sixhourly_accum_to_daily_wm2(sl[v])[:NDAYS_AUG]
                    warm_mean = daily_wm2[warm_days].mean(axis=0)
                    if not lat_ascending:
                        warm_mean = warm_mean[::-1, :]
                    sl_clim[v] += warm_mean
                sl_count += 1
            except Exception as e:
                print(f"    Warning: skip SL {yr}: {e}")

        for k in sl_clim:
            if sl_count > 0:
                sl_clim[k] /= sl_count

        # Event year SL
        sl_event, _, _ = _load_sl(data_dir, event_year, month)
        sl_ev = {}
        for v in ('skt', 'sp'):
            daily = _sixhourly_to_daily_mean(sl_event[v])[:NDAYS_AUG]
            sl_ev[v] = daily[warm_days].mean(axis=0)
            if not lat_ascending:
                sl_ev[v] = sl_ev[v][::-1, :]
        for v in ('ssrd', 'ssr', 'tisr', 'slhf', 'sshf'):
            if v not in sl_event:
                continue
            daily_wm2 = _sixhourly_accum_to_daily_wm2(sl_event[v])[:NDAYS_AUG]
            sl_ev[v] = daily_wm2[warm_days].mean(axis=0)
            if not lat_ascending:
                sl_ev[v] = sl_ev[v][::-1, :]

        # Map to CFRAM variable names
        for label, sl_data in [('base', sl_clim), ('perturbed', sl_ev)]:
            state = base_state if label == 'base' else pert_state
            state['ts'] = sl_data['skt']
            state['ps'] = sl_data['sp']
            state['solar'] = sl_data['tisr']
            # Albedo = (ssrd - ssr) / ssrd  (upwelling fraction)
            ssrd = sl_data['ssrd']
            ssr = sl_data['ssr']
            albedo = np.where(ssrd > 0, (ssrd - ssr) / ssrd, 0.2)
            albedo = np.clip(albedo, 0.0, 1.0)
            state['albedo'] = albedo

        # ── Non-radiative forcing (surface fluxes) ──────────────────────
        # Forcing = event_year - climatology (W/m²), only at surface level
        nonrad = {}
        has_flux = ('slhf' in sl_ev and 'sshf' in sl_ev and
                    np.any(sl_clim['slhf'] != 0))
        if has_flux:
            # ERA5 convention: slhf/sshf are negative for upward (energy loss)
            # CFRAM forcing convention: positive = warming
            # lhflx forcing = Δ(slhf) = event - clim
            nonrad['lhflx'] = sl_ev['slhf'] - sl_clim['slhf']
            nonrad['shflx'] = sl_ev['sshf'] - sl_clim['sshf']
            print(f"  lhflx forcing: mean={nonrad['lhflx'].mean():.3f} W/m2")
            print(f"  shflx forcing: mean={nonrad['shflx'].mean():.3f} W/m2")
        else:
            print("  WARNING: slhf/sshf not available, no nonrad forcing")

        # ── Step 4: CO2 ──────────────────────────────────────────────────
        co2_base = self.get_co2('base')
        co2_pert = self.get_co2('perturbed')
        shape_3d = (nlev, nlat, nlon)
        base_state['co2'] = np.full(shape_3d, co2_base, dtype=np.float64)
        pert_state['co2'] = np.full(shape_3d, co2_pert, dtype=np.float64)
        print(f"  CO2: base={co2_base*1e6:.1f} ppmv, "
              f"perturbed={co2_pert*1e6:.1f} ppmv")

        # ── Step 5: Aerosol ──────────────────────────────────────────────
        aer_src = self.source_cfg.get('aerosol', {}).get('source', 'zero')
        aer_kwargs = dict(
            target_plev_hpa=lev, target_lat=lat, target_lon=lon,
            years=clim_years, month=month, warm_days=warm_days)
        print(f"  Aerosol: {aer_src}")
        aer_base = self.get_aerosol(shape_3d, period='clim', **aer_kwargs)
        aer_pert = self.get_aerosol(shape_3d, period=event_year, **aer_kwargs)
        base_state.update(aer_base)
        pert_state.update(aer_pert)

        # ── Summary ──────────────────────────────────────────────────────
        print(f"\nBase state summary:")
        print(f"  T range: [{base_state['ta_lay'].min():.1f}, "
              f"{base_state['ta_lay'].max():.1f}] K")
        print(f"  ts range: [{base_state['ts'].min():.1f}, "
              f"{base_state['ts'].max():.1f}] K")
        print(f"  ps range: [{base_state['ps'].min():.0f}, "
              f"{base_state['ps'].max():.0f}] Pa")
        print(f"  albedo range: [{base_state['albedo'].min():.3f}, "
              f"{base_state['albedo'].max():.3f}]")

        dT_sfc = pert_state['ts'] - base_state['ts']
        print(f"\nPerturbed - Base (surface T):")
        print(f"  dTs range: [{dT_sfc.min():.2f}, {dT_sfc.max():.2f}] K")
        print(f"  dTs mean:  {dT_sfc.mean():.2f} K")

        return base_state, pert_state, nonrad


def _parse_mmdd(value):
    month, day = str(value).split('-')
    return int(month), int(day)


def _date_pairs_for_year(year, start_mmdd, end_mmdd):
    sm, sd = _parse_mmdd(start_mmdd)
    em, ed = _parse_mmdd(end_mmdd)
    start = date(year, sm, sd)
    end = date(year, em, ed)
    if end < start:
        raise ValueError(f"date range crosses year boundary: {start_mmdd}..{end_mmdd}")
    out = []
    cur = start
    while cur <= end:
        out.append((cur.month, cur.day))
        cur += timedelta(days=1)
    return out


def _select_daily_mean_for_dates(loader, data_dir, var_or_year, year_or_dates,
                                 dates=None):
    """Load monthly ERA5 data and return selected daily means for dates."""
    if dates is None:
        # Surface form: loader(data_dir, year, month)
        year = var_or_year
        date_pairs = year_or_dates
        var_short = None
    else:
        # PL form: loader(data_dir, var_short, year, month)
        var_short = var_or_year
        year = year_or_dates
        date_pairs = dates

    pieces = []
    for month in sorted(set(m for m, _ in date_pairs)):
        days = [d for m, d in date_pairs if m == month]
        if var_short is None:
            sl, _, _ = loader(data_dir, year, month)
            month_daily = {}
            for key, arr in sl.items():
                if key in ('ssrd', 'ssr', 'tisr', 'slhf', 'sshf'):
                    month_daily[key] = _sixhourly_accum_to_daily_wm2(arr)
                else:
                    month_daily[key] = _sixhourly_to_daily_mean(arr)
            pieces.append((month, days, month_daily))
        else:
            arr, _, _, _ = loader(data_dir, var_short, year, month)
            daily = _sixhourly_to_daily_mean(arr)
            pieces.append((month, days, daily))
    return pieces


@register_source('era5_date_range')
class ERA5DateRangeSource(DataSource):
    """Build CFRAM states from ERA5 over an explicit MM-DD date range."""

    def build_states(self):
        tcfg = self.source_cfg['temporal']
        event_year = int(tcfg['event_year'])
        clim_start, clim_end = tcfg['clim_years']
        clim_years = list(range(clim_start, clim_end + 1))
        start_mmdd = tcfg['start_mmdd']
        end_mmdd = tcfg['end_mmdd']

        data_dir = self.source_cfg['data_dir']
        if not os.path.isabs(data_dir):
            from core.config import PROJECT_ROOT
            data_dir = os.path.join(PROJECT_ROOT, data_dir)

        date_pairs_event = _date_pairs_for_year(event_year, start_mmdd, end_mmdd)
        year_dates = {yr: _date_pairs_for_year(yr, start_mmdd, end_mmdd)
                      for yr in clim_years}
        year_dates[event_year] = date_pairs_event
        print(f"Date-range period: {start_mmdd}..{end_mmdd} ({len(date_pairs_event)} days)")

        first_month = date_pairs_event[0][0]
        _, lev, lat_raw, lon_raw = _load_pl(data_dir, 't', event_year, first_month)
        lat_ascending = lat_raw[0] < lat_raw[-1]
        lat = lat_raw if lat_ascending else lat_raw[::-1]
        lon = lon_raw

        nlev, nlat, nlon = len(lev), len(lat), len(lon)
        print(f"Grid: {nlev} levels, {nlat} lat, {nlon} lon")

        base_state = {'lat': lat, 'lon': lon, 'lev': lev}
        pert_state = {'lat': lat, 'lon': lon, 'lev': lev}

        for era5_var, cfram_var in PL_VAR_MAP.items():
            print(f"  Processing PL: {era5_var} → {cfram_var}")
            clim_sum = np.zeros((nlev, nlat, nlon), dtype=np.float64)
            clim_count = 0
            for yr in clim_years:
                selected = []
                for month, days, daily in _select_daily_mean_for_dates(
                        _load_pl, data_dir, era5_var, yr, year_dates[yr]):
                    idx = [d - 1 for d in days]
                    selected.append(daily[idx])
                yr_mean = np.concatenate(selected, axis=0).mean(axis=0)
                if not lat_ascending:
                    yr_mean = yr_mean[:, ::-1, :]
                clim_sum += yr_mean
                clim_count += 1
            base_state[cfram_var] = clim_sum / clim_count

            selected_ev = []
            for month, days, daily in _select_daily_mean_for_dates(
                    _load_pl, data_dir, era5_var, event_year, date_pairs_event):
                selected_ev.append(daily[[d - 1 for d in days]])
            pert_mean = np.concatenate(selected_ev, axis=0).mean(axis=0)
            if not lat_ascending:
                pert_mean = pert_mean[:, ::-1, :]
            pert_state[cfram_var] = pert_mean

        print("  Processing SL variables...")
        sl_keys = ['skt', 'sp', 'ssrd', 'ssr', 'tisr', 'slhf', 'sshf']
        sl_clim = {k: np.zeros((nlat, nlon), dtype=np.float64) for k in sl_keys}
        sl_count = 0
        for yr in clim_years:
            selected = {k: [] for k in sl_keys}
            for month, days, month_daily in _select_daily_mean_for_dates(
                    _load_sl, data_dir, yr, year_dates[yr]):
                idx = [d - 1 for d in days]
                for key in sl_keys:
                    if key in month_daily:
                        selected[key].append(month_daily[key][idx])
            for key in sl_keys:
                if selected[key]:
                    mean = np.concatenate(selected[key], axis=0).mean(axis=0)
                    if not lat_ascending:
                        mean = mean[::-1, :]
                    sl_clim[key] += mean
            sl_count += 1
        for key in sl_keys:
            sl_clim[key] /= sl_count

        selected_ev = {k: [] for k in sl_keys}
        for month, days, month_daily in _select_daily_mean_for_dates(
                _load_sl, data_dir, event_year, date_pairs_event):
            idx = [d - 1 for d in days]
            for key in sl_keys:
                if key in month_daily:
                    selected_ev[key].append(month_daily[key][idx])
        sl_ev = {}
        for key in sl_keys:
            if selected_ev[key]:
                mean = np.concatenate(selected_ev[key], axis=0).mean(axis=0)
                if not lat_ascending:
                    mean = mean[::-1, :]
                sl_ev[key] = mean

        for label, sl_data in [('base', sl_clim), ('perturbed', sl_ev)]:
            state = base_state if label == 'base' else pert_state
            state['ts'] = sl_data['skt']
            state['ps'] = sl_data['sp']
            state['solar'] = sl_data['tisr']
            ssrd = sl_data['ssrd']
            ssr = sl_data['ssr']
            albedo = np.where(ssrd > 0, (ssrd - ssr) / ssrd, 0.2)
            state['albedo'] = np.clip(albedo, 0.0, 1.0)

        nonrad = {}
        if 'slhf' in sl_ev and 'sshf' in sl_ev:
            nonrad['lhflx'] = np.nan_to_num(sl_ev['slhf'] - sl_clim['slhf'])
            nonrad['shflx'] = np.nan_to_num(sl_ev['sshf'] - sl_clim['sshf'])
            print(f"  lhflx forcing: mean={nonrad['lhflx'].mean():.3f} W/m2")
            print(f"  shflx forcing: mean={nonrad['shflx'].mean():.3f} W/m2")

        co2_base = self.get_co2('base')
        co2_pert = self.get_co2('perturbed')
        shape_3d = (nlev, nlat, nlon)
        base_state['co2'] = np.full(shape_3d, co2_base, dtype=np.float64)
        pert_state['co2'] = np.full(shape_3d, co2_pert, dtype=np.float64)
        print(f"  CO2: base={co2_base*1e6:.1f} ppmv, perturbed={co2_pert*1e6:.1f} ppmv")

        aer_src = self.source_cfg.get('aerosol', {}).get('source', 'zero')
        print(f"  Aerosol: {aer_src}")
        if aer_src == 'merra2':
            from .merra2_aerosol import load_merra2_aerosol_dates
            aer_cfg = self.source_cfg.get('aerosol', {})
            m2_dir = aer_cfg.get('data_dir', 'era5_data/merra2')
            if not os.path.isabs(m2_dir):
                from core.config import PROJECT_ROOT
                m2_dir = os.path.join(PROJECT_ROOT, m2_dir)
            aer_base = load_merra2_aerosol_dates(
                m2_dir, {yr: year_dates[yr] for yr in clim_years},
                lev, lat, lon, period='clim')
            aer_pert = load_merra2_aerosol_dates(
                m2_dir, {event_year: date_pairs_event},
                lev, lat, lon, period=event_year)
        else:
            aer_base = self.get_aerosol(shape_3d, period='base')
            aer_pert = self.get_aerosol(shape_3d, period='perturbed')
        base_state.update(aer_base)
        pert_state.update(aer_pert)

        dT_sfc = pert_state['ts'] - base_state['ts']
        print(f"\nPerturbed - Base (surface T):")
        print(f"  dTs range: [{dT_sfc.min():.2f}, {dT_sfc.max():.2f}] K")
        print(f"  dTs mean:  {dT_sfc.mean():.2f} K")

        return base_state, pert_state, nonrad
