#!/usr/bin/env python3
"""Embarrassingly parallel CFRAM using multiprocessing.

Each worker:
1. Extracts single column from pre-loaded global data
2. Writes temp binary files to a unique temp dir
3. Runs single-column cfram_rrtmg executable
4. Reads forcing + Planck matrix, solves dT in Python

Usage:
    python3 scripts/run_parallel_python.py --case eh13 --nproc 40
"""
import os, sys, argparse, struct, tempfile, shutil
import numpy as np
from multiprocessing import Pool
from netCDF4 import Dataset

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.config import load_case, defaults, get_plev, get_aerosol_map, get_nproc, \
    get_fortran_dir, get_lookup_dir, get_executable, get_output_terms, \
    get_drdt_eval, get_drdt_probe, get_q_handling, get_co2_handling, PROJECT_ROOT
from core.constants import NBND_LW, NBND_SW

# Default plev / nlev (used until main() overrides with case-specific values).
# Per-case override flows through init_worker() so multiprocessing workers
# (both fork on Linux and spawn on macOS) see the correct nlev.
FORTRAN_PLEV = get_plev()
NLEV = len(FORTRAN_PLEV)
SCON = defaults()['radiation']['scon']
FORTRAN_EXE_DIR = get_fortran_dir()
LOOKUP_DIR = get_lookup_dir()

# Build aerosol map from config
_aer_cfg = get_aerosol_map()
AEROSOL_MAP = {k: (v['lw'], v['sw'], v['rd']) for k, v in _aer_cfg.items()}

# Canonical species order for per-species aerosol perturbation output.
# Must match the order expected by Fortran when reading *_spc.dat files.
SPECIES_ORDER = ['bc', 'ocphi', 'ocpho', 'sulf', 'ss', 'dust']
NSPECIES = len(SPECIES_ORDER)

# Global data (set by initializer)
G = {}


def load_lookup_table(nc_path):
    nc = Dataset(nc_path, 'r')
    table = {v: np.array(nc.variables[v][:], dtype=np.float64) for v in ('bext','bsca','qext','qsca','g','rh','radius') if v in nc.variables}
    nc.close()
    return table


def compute_rh(q, t, p_hpa):
    w = q / (1.0 - q + 1e-30)
    e = w * p_hpa / (0.622 + w)
    es = 6.112 * np.exp(17.67 * (t - 273.15) / (t - 29.65 + 1e-10))
    return np.clip(e / (es + 1e-30), 0.0, 1.0)


def compute_aerosol_column(aer_data, rh, thick, dens):
    """Compute aerosol optics for single column.

    Returns:
        aod_lw, aod_sw, eff_ssa, eff_g: bulk (all-species) arrays (existing)
        per_species: dict {species: {'aod_lw','aod_sw','ssa_sw','g_sw'}} where
                     ssa_sw/g_sw are the raw (per-species, unweighted) values.
    """
    LW_BANDS = slice(14, 30)
    SW_BANDS = slice(0, 14)
    aod_lw = np.zeros((NLEV, NBND_LW))
    aod_sw = np.zeros((NLEV, NBND_SW))
    aod_ssa = np.zeros((NLEV, NBND_SW))
    aod_g = np.zeros((NLEV, NBND_SW))
    per_species = {}

    for species, (lut_lw, lut_sw, rd) in AEROSOL_MAP.items():
        if species not in aer_data:
            continue
        mixing = aer_data[species]
        sp_aod_lw = np.zeros((NLEV, NBND_LW))
        sp_aod_sw = np.zeros((NLEV, NBND_SW))
        sp_ssa_sw = np.zeros((NLEV, NBND_SW))
        sp_g_sw = np.zeros((NLEV, NBND_SW))
        for lut_file, bands, nbnd, target in [(lut_lw, LW_BANDS, NBND_LW, 'lw'), (lut_sw, SW_BANDS, NBND_SW, 'sw')]:
            table = G['lut'][lut_file]
            r_idx = int(rd) - 1 if rd >= 1 else int(np.argmin(np.abs(table['radius'] - rd)))
            bext = table['bext'][r_idx, :, bands]
            rh_idx = np.array([np.argmin(np.abs(table['rh'] - r)) for r in rh])
            kext = np.array([bext[rh_idx[i], :] for i in range(NLEV)])
            aod_col = mixing[:, None] * dens[:, None] * 1e3 * thick[:, None] * kext * 1e-3
            aod_col = np.nan_to_num(aod_col)
            if target == 'lw':
                aod_lw += aod_col
                sp_aod_lw = aod_col
            else:
                aod_sw += aod_col
                qext = table['qext'][r_idx, :, bands]
                qsca = table['qsca'][r_idx, :, bands]
                g_tab = table['g'][r_idx, :, bands]
                ssa_col = np.array([np.where(qext[rh_idx[i], :] > 0, qsca[rh_idx[i], :] / qext[rh_idx[i], :], 0.0) for i in range(NLEV)])
                g_col = np.array([g_tab[rh_idx[i], :] for i in range(NLEV)])
                aod_ssa += aod_col * ssa_col
                aod_g += aod_col * ssa_col * g_col
                sp_aod_sw = aod_col
                sp_ssa_sw = np.nan_to_num(ssa_col)
                sp_g_sw = np.nan_to_num(g_col)
        per_species[species] = {
            'aod_lw': sp_aod_lw,
            'aod_sw': sp_aod_sw,
            'ssa_sw': sp_ssa_sw,
            'g_sw': sp_g_sw,
        }

    with np.errstate(divide='ignore', invalid='ignore'):
        eff_ssa = np.where(aod_sw > 0, aod_ssa / aod_sw, 0.0)
        eff_g = np.where(aod_ssa > 0, aod_g / aod_ssa, 0.0)
    return aod_lw, aod_sw, eff_ssa, eff_g, per_species


def stack_species(per_species, field, shape):
    """Stack per-species arrays into (*shape, NSPECIES), zero-fill missing species.

    Species is the LAST axis (innermost in C-order on disk). This matches the
    Fortran reader, which uses `isp` as the innermost implied-DO loop index.
    """
    arr = np.zeros(tuple(shape) + (NSPECIES,), dtype=np.float64)
    for i, s in enumerate(SPECIES_ORDER):
        if s in per_species:
            arr[..., i] = per_species[s][field]
    return arr


def write_bin(fp, data):
    np.asarray(data, dtype=np.float64).tofile(fp)


def write_aer_bin(fp, data):
    """Write aerosol array in C-order (band index fastest).

    Note: this is NOT Fortran column-major, but matches Phase 3 convention
    and gives better correlation with Wu et al. paper results.
    """
    np.asarray(data, dtype=np.float64).tofile(fp)


def read_fortran_seq(fp):
    with open(fp, 'rb') as f:
        raw = f.read()
    n = struct.unpack('i', raw[:4])[0]
    return np.frombuffer(raw[4:4+n], dtype=np.float64).copy()


def process_column(args):
    """Process single column. Returns (ilat, ilon, dT_dict, frc_dict)."""
    ilat, ilon = args
    d = G  # global data

    # Extract column (already TOA->surface)
    t_b = d['t_base'][:, ilat, ilon]
    q_b = d['q_base'][:, ilat, ilon]
    o3_b = d['o3_base'][:, ilat, ilon]
    cc_b = d['cc_base'][:, ilat, ilon]
    clwc_b = d['clwc_base'][:, ilat, ilon]
    ciwc_b = d['ciwc_base'][:, ilat, ilon]
    t_w = d['t_warm'][:, ilat, ilon]
    q_w = d['q_warm'][:, ilat, ilon]
    o3_w = d['o3_warm'][:, ilat, ilon]
    cc_w = d['cc_warm'][:, ilat, ilon]
    clwc_w = d['clwc_warm'][:, ilat, ilon]
    ciwc_w = d['ciwc_warm'][:, ilat, ilon]

    ts_b = d['ts_base'][ilat, ilon]
    ps_b = d['ps_base'][ilat, ilon]
    solar_b = d['solar_base'][ilat, ilon]
    albedo_b = d['albedo_base'][ilat, ilon]
    huss_b = d['huss_base'][ilat, ilon]
    ts_w = d['ts_warm'][ilat, ilon]
    ps_w = d['ps_warm'][ilat, ilon]
    solar_w = d['solar_warm'][ilat, ilon]
    albedo_w = d['albedo_warm'][ilat, ilon]
    huss_w = d['huss_warm'][ilat, ilon]

    # Aerosol optics
    aer_b = {s: d['aer_base_'+s][:, ilat, ilon] for s in AEROSOL_MAP}
    aer_w = {s: d['aer_warm_'+s][:, ilat, ilon] for s in AEROSOL_MAP}

    p_hpa = FORTRAN_PLEV
    Rd = 287.05
    p_half = np.zeros(NLEV + 1)
    p_half[0] = 0.0
    for i in range(NLEV - 1):
        p_half[i + 1] = 0.5 * (p_hpa[i] + p_hpa[i + 1])
    p_half[-1] = ps_b / 100.0
    delp = np.diff(p_half) * 100.0
    dens_b = (p_hpa * 100.0) / (Rd * t_b)
    thick_b = delp / (dens_b * 9.8)
    rh_b = compute_rh(q_b, t_b, p_hpa)

    aod_lw_b, aod_sw_b, ssa_sw_b, g_sw_b, perspc_b = compute_aerosol_column(aer_b, rh_b, thick_b, dens_b)

    p_half_w = p_half.copy()
    p_half_w[-1] = ps_w / 100.0
    delp_w = np.diff(p_half_w) * 100.0
    dens_w = (p_hpa * 100.0) / (Rd * t_w)
    thick_w = delp_w / (dens_w * 9.8)
    rh_w = compute_rh(q_w, t_w, p_hpa)
    aod_lw_w, aod_sw_w, ssa_sw_w, g_sw_w, perspc_w = compute_aerosol_column(aer_w, rh_w, thick_w, dens_w)

    # Write to temp dir
    tmpdir = tempfile.mkdtemp(prefix='cfram_')
    dp = os.path.join(tmpdir, 'data_prep')
    do = os.path.join(tmpdir, 'data_output')
    os.makedirs(dp)
    os.makedirs(do)

    ssrd_b = 300.0
    ssru_b = ssrd_b * albedo_b
    ssrd_w = 300.0
    ssru_w = ssrd_w * albedo_w

    # plev.dat: vertical grid (Fortran reads this; same for both base and warm)
    write_bin(os.path.join(dp, 'plev.dat'), FORTRAN_PLEV)
    write_bin(os.path.join(dp, 't_base.dat'), t_b)
    write_bin(os.path.join(dp, 't_warm.dat'), t_w)
    write_bin(os.path.join(dp, 'hus_base.dat'), q_b)
    write_bin(os.path.join(dp, 'hus_warm.dat'), q_w)
    write_bin(os.path.join(dp, 'O3_base.dat'), o3_b)
    write_bin(os.path.join(dp, 'O3_warm.dat'), o3_w)
    write_bin(os.path.join(dp, 'cc_base.dat'), cc_b)
    write_bin(os.path.join(dp, 'cc_warm.dat'), cc_w)
    write_bin(os.path.join(dp, 'clwc_base.dat'), clwc_b)
    write_bin(os.path.join(dp, 'clwc_warm.dat'), clwc_w)
    write_bin(os.path.join(dp, 'ciwc_base.dat'), ciwc_b)
    write_bin(os.path.join(dp, 'ciwc_warm.dat'), ciwc_w)
    write_bin(os.path.join(dp, 'skt_base.dat'), np.array([ts_b]))
    write_bin(os.path.join(dp, 'skt_warm.dat'), np.array([ts_w]))
    write_bin(os.path.join(dp, 'sp_base.dat'), np.array([ps_b]))
    write_bin(os.path.join(dp, 'sp_warm.dat'), np.array([ps_w]))
    write_bin(os.path.join(dp, 'huss_base.dat'), np.array([huss_b]))
    write_bin(os.path.join(dp, 'huss_warm.dat'), np.array([huss_w]))
    write_bin(os.path.join(dp, 'solarin_base.dat'), np.array([solar_b]))
    write_bin(os.path.join(dp, 'solarin_warm.dat'), np.array([solar_w]))
    write_bin(os.path.join(dp, 'ssrd_base.dat'), np.array([ssrd_b]))
    write_bin(os.path.join(dp, 'ssrd_warm.dat'), np.array([ssrd_w]))
    write_bin(os.path.join(dp, 'ssru_base.dat'), np.array([ssru_b]))
    write_bin(os.path.join(dp, 'ssru_warm.dat'), np.array([ssru_w]))
    write_bin(os.path.join(dp, 'co2_b.dat'), np.array([d['co2_base']]))
    write_bin(os.path.join(dp, 'co2_w.dat'), np.array([d['co2_warm']]))
    write_aer_bin(os.path.join(dp, 'aerosol_aod_lw_base.dat'), aod_lw_b)
    write_aer_bin(os.path.join(dp, 'aerosol_aod_lw_warm.dat'), aod_lw_w)
    write_aer_bin(os.path.join(dp, 'aerosol_aod_sw_base.dat'), aod_sw_b)
    write_aer_bin(os.path.join(dp, 'aerosol_aod_sw_warm.dat'), aod_sw_w)
    write_aer_bin(os.path.join(dp, 'aerosol_ssa_sw_base.dat'), ssa_sw_b)
    write_aer_bin(os.path.join(dp, 'aerosol_ssa_sw_warm.dat'), ssa_sw_w)
    write_aer_bin(os.path.join(dp, 'aerosol_g_sw_base.dat'), g_sw_b)
    write_aer_bin(os.path.join(dp, 'aerosol_g_sw_warm.dat'), g_sw_w)

    # Skip flags for Fortran (Step 1: aerosol-only). Fortran reads skip_flags.txt
    # and skips the bulk + per-species aerosol rad_driver calls when set.
    if d.get('skip_aerosol', False):
        with open(os.path.join(dp, 'skip_flags.txt'), 'w') as fh:
            fh.write('skip_aerosol=1\n')

    # Optional 2nd-order CFRAM: evaluate Planck matrix at midstate
    # (T_base+T_warm)/2 instead of base state. RRTMG checks for the
    # existence of this file; Fu currently ignores it.
    if d.get('drdt_eval') == 'midstate':
        open(os.path.join(dp, 'drdt_midstate.flag'), 'w').close()
    # Optional Manabe q-feedback: compute frc_q at warm-T atmosphere,
    # adding 1 rad_driver call for the warm-T q-base reference.
    if d.get('q_handling') == 'feedback':
        open(os.path.join(dp, 'q_feedback.flag'), 'w').close()
    # Optional q midstate path: compute frc_q in midstate atmosphere
    # (cancels ∂²R/∂T∂q cross-term). +2 rad_driver calls per cell.
    if d.get('q_handling') == 'midstate':
        open(os.path.join(dp, 'q_midstate.flag'), 'w').close()
    # Optional CO2 midstate path: compute frc_co2 in midstate atmosphere
    # (T_mid, q_mid, ...) instead of base, cancelling ∂²R/∂T∂C cross-term.
    # +2 rad_driver calls per cell.
    if d.get('co2_handling') == 'midstate':
        open(os.path.join(dp, 'co2_midstate.flag'), 'w').close()
    # Optional centered-FD Planck probe: T_j ± 0.5K instead of T_j + 1K.
    # Cancels R_TT term in the FD expansion. 2× rad_driver_lw in calc_drdt.
    if d.get('drdt_probe') == 'centered':
        open(os.path.join(dp, 'drdt_centered.flag'), 'w').close()

    # Per-species optical properties: shape (NLEV, NBND, NSPECIES), C-order.
    # Species is the LAST (fastest-varying) axis, matching Fortran READ where
    # `isp` is the innermost implied-DO loop index. Species order is
    # SPECIES_ORDER = ['bc','ocphi','ocpho','sulf','ss','dust'].
    write_aer_bin(os.path.join(dp, 'aerosol_aod_lw_base_spc.dat'),
                  stack_species(perspc_b, 'aod_lw', (NLEV, NBND_LW)))
    write_aer_bin(os.path.join(dp, 'aerosol_aod_lw_warm_spc.dat'),
                  stack_species(perspc_w, 'aod_lw', (NLEV, NBND_LW)))
    write_aer_bin(os.path.join(dp, 'aerosol_aod_sw_base_spc.dat'),
                  stack_species(perspc_b, 'aod_sw', (NLEV, NBND_SW)))
    write_aer_bin(os.path.join(dp, 'aerosol_aod_sw_warm_spc.dat'),
                  stack_species(perspc_w, 'aod_sw', (NLEV, NBND_SW)))
    write_aer_bin(os.path.join(dp, 'aerosol_ssa_sw_base_spc.dat'),
                  stack_species(perspc_b, 'ssa_sw', (NLEV, NBND_SW)))
    write_aer_bin(os.path.join(dp, 'aerosol_ssa_sw_warm_spc.dat'),
                  stack_species(perspc_w, 'ssa_sw', (NLEV, NBND_SW)))
    write_aer_bin(os.path.join(dp, 'aerosol_g_sw_base_spc.dat'),
                  stack_species(perspc_b, 'g_sw', (NLEV, NBND_SW)))
    write_aer_bin(os.path.join(dp, 'aerosol_g_sw_warm_spc.dat'),
                  stack_species(perspc_w, 'g_sw', (NLEV, NBND_SW)))

    # Symlink single-column executable
    exe_1col = d['exe_1col']
    os.symlink(exe_1col, os.path.join(tmpdir, 'cfram_rrtmg'))

    # Run
    ret = os.system('cd %s && ./cfram_rrtmg > /dev/null 2>&1' % tmpdir)

    # Read output: forcing + drdt_inv, solve dT in Python
    result = {}
    rad_terms = ['co2', 'q', 'ts', 'o3', 'solar', 'albedo', 'cloud', 'aerosol', 'warm',
                 'cloud_lw', 'cloud_sw', 'full']
    nonrad_terms = ['lhflx', 'shflx']
    dyn_terms = ['atmdyn', 'sfcdyn', 'ocndyn', 'dry']
    # Legacy path-B species (from paper_data nonrad_forcing.nc), 5 terms, OC merged.
    aer_species_terms = ['bc', 'oc', 'sulf', 'seas', 'dust']
    # Fortran-produced per-species forcings (Phase 3), 6 terms with OC split.
    fortran_aer_species = SPECIES_ORDER  # ['bc','ocphi','ocpho','sulf','ss','dust']

    if ret == 0:
        # Read forcing from Fortran
        for t in rad_terms:
            try:
                result['frc_'+t] = read_fortran_seq(os.path.join(do, 'frc_%s.dat' % t))
            except:
                result['frc_'+t] = np.full(NLEV+1, np.nan)

        # Fortran per-species aerosol forcings (if binary was compiled with Phase 3).
        # Present: override any path-B value for overlapping names.
        fortran_has_species = False
        for t in fortran_aer_species:
            path = os.path.join(do, 'frc_%s.dat' % t)
            if os.path.exists(path):
                try:
                    result['frc_'+t] = read_fortran_seq(path)
                    fortran_has_species = True
                except Exception:
                    pass

        # Read drdt_inv matrix
        try:
            drdt_file = os.path.join(do, 'drdt_inv.dat')
            with open(drdt_file, 'rb') as f:
                raw = f.read()
            pos = 0
            rec1_len = struct.unpack('i', raw[pos:pos+4])[0]; pos += 4
            if rec1_len == 8:
                nlayer_val = struct.unpack('q', raw[pos:pos+8])[0]; pos += 8
            else:
                nlayer_val = struct.unpack('i', raw[pos:pos+4])[0]; pos += 4
            pos += 4  # skip footer
            rec2_len = struct.unpack('i', raw[pos:pos+4])[0]; pos += 4
            n = nlayer_val + 1
            mat = np.frombuffer(raw[pos:pos+rec2_len], dtype=np.float64).copy().reshape(n, n, order='F')
        except:
            mat = None
            nlayer_val = 0

        # Solve dT = -drdt_inv @ frc for ALL terms (radiative + non-radiative)
        if mat is not None:
            nl = nlayer_val
            n = nl + 1

            def solve_dT(frc_full):
                """Apply Planck inverse to forcing vector."""
                frc_col = np.zeros(n)
                frc_col[:nl] = frc_full[:nl]
                frc_col[nl] = frc_full[NLEV]  # surface
                dT_vec = -mat @ frc_col  # n-dimensional result
                # Map back to full NLEV+1 array
                dT_full = np.full(NLEV + 1, -999.0)
                dT_full[:nl] = dT_vec[:nl]
                dT_full[NLEV] = dT_vec[nl]  # surface
                return dT_full

            # Radiative terms: forcing from Fortran
            for t in rad_terms:
                result['dT_'+t] = solve_dT(result['frc_'+t])

            # Dynamic terms derived from energy conservation: F_dyn = -F_rad (= -frc_warm)
            frc_warm_col = result['frc_warm']

            # Atmospheric dynamics (OLD CFRAM Cai-Tung Fortran convention):
            #   F_atmdyn[atm] = -ΔR_full[atm], surface = 0
            # where ΔR_full uses T_warm in the radiation call (full-state),
            # NOT T_base. This includes the drdt·ΔT_obs Planck response, which
            # is essential to capture polar amplification's heat-transport signal
            # (without it, polar atm_dyn flips sign — see lump_legacy below).
            frc_full_col = result.get('frc_full')
            if frc_full_col is not None and np.isfinite(frc_full_col).any():
                frc_atmdyn = np.zeros(NLEV + 1)
                frc_atmdyn[:nl] = -frc_full_col[:nl]
                result['dT_atmdyn'] = solve_dT(frc_atmdyn)
                # Lu/Cai full-state Q_dyn (no sfc zeroing): dT_dry = drdt^-1·frc_full.
                # Closure-exact: Σ_X dT_X = ΔT_obs at first order.
                result['dT_dry'] = solve_dT(-frc_full_col)
            else:
                # Fall back to legacy lumping if frc_full unavailable
                frc_atmdyn = np.zeros(NLEV + 1)
                frc_atmdyn[:nl] = -frc_warm_col[:nl]
                result['dT_atmdyn'] = solve_dT(frc_atmdyn)
                result['dT_dry'] = np.full(NLEV + 1, np.nan)

            # Surface dynamics — all dyn terms use frc_full (T_warm-state ΔR) to
            # match OLD CFRAM Cai-Tung convention. Falls back to frc_warm if
            # frc_full unavailable.
            #   sfcdyn       = -ΔR_full[sfc] only (atm = 0)
            #   ocndyn       = -ΔR_full[sfc] - Δlh - Δsh  ≡ OLD dt_sfc_dyn
            #   atmdyn       = -ΔR_full[atm], sfc = 0      ≡ OLD dt_atm_dyn
            #   dry          = atmdyn + sfcdyn (full column)  ≡ OLD dt_dyn
            sfc_forcing_src = (frc_full_col if (frc_full_col is not None
                                                and np.isfinite(frc_full_col).any())
                               else frc_warm_col)
            frc_sfcdyn = np.zeros(NLEV + 1)
            frc_sfcdyn[NLEV] = -sfc_forcing_src[NLEV]
            result['dT_sfcdyn'] = solve_dT(frc_sfcdyn)

            # Ocean circulation ≡ OLD CFRAM dt_sfc_dyn:
            # F_ocndyn[sfc] = -ΔR_full[sfc] - Δlh - Δsh
            lhflx_sfc = d['frc_lhflx'][NLEV, ilat, ilon] if 'frc_lhflx' in d else 0.0
            shflx_sfc = d['frc_shflx'][NLEV, ilat, ilon] if 'frc_shflx' in d else 0.0
            frc_ocndyn = np.zeros(NLEV + 1)
            frc_ocndyn[NLEV] = frc_sfcdyn[NLEV] - lhflx_sfc - shflx_sfc
            result['dT_ocndyn'] = solve_dT(frc_ocndyn)

            # lhflx, shflx + path-B per-species aerosol: forcing from paper_data.
            # Fortran per-species (if present) already overwrote result['frc_<spc>'];
            # path-B is used only for species whose forcing is NOT from Fortran.
            for t in nonrad_terms:
                if 'frc_' + t in d:
                    result['dT_' + t] = solve_dT(d['frc_' + t][:, ilat, ilon])
            for t in aer_species_terms:
                # Skip if Fortran already provided this species name (bc/sulf/dust overlap).
                if t in result and 'frc_' + t in result:
                    continue
                if 'frc_' + t in d:
                    result['dT_' + t] = solve_dT(d['frc_' + t][:, ilat, ilon])
            # Apply Planck inverse to Fortran per-species forcings.
            for t in fortran_aer_species:
                if 'frc_' + t in result:
                    result['dT_' + t] = solve_dT(result['frc_' + t])
        else:
            for t in rad_terms + dyn_terms + nonrad_terms + aer_species_terms + fortran_aer_species:
                result['dT_'+t] = np.full(NLEV+1, np.nan)
    else:
        for t in rad_terms:
            result['frc_'+t] = np.full(NLEV+1, np.nan)
            result['dT_'+t] = np.full(NLEV+1, np.nan)
        for t in dyn_terms + nonrad_terms + aer_species_terms + fortran_aer_species:
            result['dT_'+t] = np.full(NLEV+1, np.nan)

    shutil.rmtree(tmpdir)
    return (ilat, ilon, result)


def init_worker(data_dict):
    global G, FORTRAN_PLEV, NLEV
    G = data_dict
    # Override module globals with case-specific values (works on both fork
    # and spawn multiprocessing). Workers reference these in compute_aerosol_column,
    # process_column, etc.
    if 'fortran_plev' in data_dict:
        FORTRAN_PLEV = data_dict['fortran_plev']
        NLEV = data_dict['nlev']


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--case', required=True, help='Case name (directory under cases/)')
    parser.add_argument('--nproc', type=int, default=None,
                        help='Number of workers (default: from config or all CPUs)')
    args = parser.parse_args()

    cfg = load_case(args.case)
    OUTDIR = cfg['_output_dir']
    os.makedirs(OUTDIR, exist_ok=True)
    nproc = args.nproc or get_nproc(cfg)

    # Per-case plev / nlev (default 37, override via case.yaml grid.pressure_levels).
    global FORTRAN_PLEV, NLEV
    FORTRAN_PLEV = get_plev(cfg)
    NLEV = len(FORTRAN_PLEV)

    # Optional: case.yaml radiation.output_terms restricts which dT/frc
    # variables are written to cfram_result.nc. None = write all (default).
    output_terms = get_output_terms(cfg)
    # Optional: 2nd-order CFRAM via midstate-Planck (RRTMG only). Default 'base'.
    drdt_eval = get_drdt_eval(cfg)
    # Optional: q as Manabe feedback (RRTMG only). Default 'independent'.
    q_handling = get_q_handling(cfg)
    # Optional: CO2 partial path (RRTMG only). Default 'base'.
    co2_handling = get_co2_handling(cfg)
    # Optional: Planck Jacobian probe scheme (RRTMG only). Default 'onesided'.
    drdt_probe = get_drdt_probe(cfg)

    print("=== Python parallel CFRAM: %s, %d procs (nlev=%d) ===" %
          (cfg.get('case_name', args.case), nproc, NLEV))
    if output_terms is not None:
        print("Output filter active: writing only dT/frc for %s" % output_terms)
    if drdt_eval == 'midstate':
        print("Planck matrix mode: midstate (2nd-order CFRAM, RRTMG only)")
    if q_handling == 'feedback':
        print("q partial mode: feedback (Manabe RH-fixed, RRTMG only)")
    if co2_handling == 'midstate':
        print("CO2 partial mode: midstate (∂²R/∂T∂C cross-term cancellation, RRTMG only)")
    if drdt_probe == 'centered':
        print("Planck probe: centered FD (T_j ± 0.5K, cancels R_TT term, RRTMG only)")

    # Pre-built single-column executable (nlat=1, nlon=1, nlev=NLEV).
    exe_name = get_executable(cfg)
    exe_1col = os.path.join(FORTRAN_EXE_DIR, exe_name)
    if not os.path.exists(exe_1col):
        print("ERROR: %s not found. Build with: cd fortran && make NLEV=%d" %
              (exe_1col, NLEV))
        sys.exit(1)
    print("Using executable: %s" % exe_1col)

    # Load data
    print("Loading input data...")
    nc_bp = Dataset(cfg['input']['base_pres'])
    nc_bs = Dataset(cfg['input']['base_surf'])
    nc_ap = Dataset(cfg['input']['perturbed_pres'])
    nc_as = Dataset(cfg['input']['perturbed_surf'])

    lats = np.array(nc_bp.variables['lat'][:])
    lons = np.array(nc_bp.variables['lon'][:])
    nlat, nlon = len(lats), len(lons)

    def get3d(nc, v): return np.array(nc.variables[v][0, ::-1, :, :], dtype=np.float64)
    def get2d(nc, v): return np.array(nc.variables[v][0, :, :], dtype=np.float64)

    # huss = 2m specific humidity (CMIP6 standard). Apple-to-apple OLD CFRAM
    # raw/CFRAM.zip GW-base.f L322-330 reads huss_base.dat for ph(nv1).
    # Optional: only present in surf NCs built by build_cesm2_official.py
    # (after huss was added to the var list). Legacy cases (cesm2_4xco2 with
    # collaborator's data) lack it — fall back to sentinel; Fortran detects
    # |huss| > 900 and uses ph(nv) HOLD instead.
    def get2d_opt(nc, v, fill_shape):
        if v in nc.variables:
            return get2d(nc, v)
        return np.full(fill_shape, -999.0, dtype=np.float64)

    ts_b_2d = get2d(nc_bs, 'ts')   # for shape

    data = {
        't_base': get3d(nc_bp, 'ta_lay'), 'q_base': get3d(nc_bp, 'q'),
        'o3_base': get3d(nc_bp, 'o3'), 'cc_base': get3d(nc_bp, 'camt'),
        'clwc_base': get3d(nc_bp, 'cliq'), 'ciwc_base': get3d(nc_bp, 'cice'),
        't_warm': get3d(nc_ap, 'ta_lay'), 'q_warm': get3d(nc_ap, 'q'),
        'o3_warm': get3d(nc_ap, 'o3'), 'cc_warm': get3d(nc_ap, 'camt'),
        'clwc_warm': get3d(nc_ap, 'cliq'), 'ciwc_warm': get3d(nc_ap, 'cice'),
        'ts_base': ts_b_2d, 'ps_base': get2d(nc_bs, 'ps'),
        'solar_base': get2d(nc_bs, 'solar'), 'albedo_base': get2d(nc_bs, 'albedo'),
        'huss_base': get2d_opt(nc_bs, 'huss', ts_b_2d.shape),
        'ts_warm': get2d(nc_as, 'ts'), 'ps_warm': get2d(nc_as, 'ps'),
        'solar_warm': get2d(nc_as, 'solar'), 'albedo_warm': get2d(nc_as, 'albedo'),
        'huss_warm': get2d_opt(nc_as, 'huss', ts_b_2d.shape),
        'co2_base': get3d(nc_bp, 'co2').mean() * 1e6,
        'co2_warm': get3d(nc_ap, 'co2').mean() * 1e6,
        'exe_1col': exe_1col,
        # Pass per-case plev / nlev to workers (works on macOS spawn too)
        'fortran_plev': FORTRAN_PLEV,
        'nlev': NLEV,
        'drdt_eval': drdt_eval,
        'drdt_probe': drdt_probe,
        'q_handling': q_handling,
        'co2_handling': co2_handling,
    }
    for s in AEROSOL_MAP:
        # Aerosols may be missing in CMIP6 raw cases; tolerate absence.
        if s in nc_bp.variables:
            data['aer_base_'+s] = get3d(nc_bp, s)
            data['aer_warm_'+s] = get3d(nc_ap, s)
        else:
            zeros = np.zeros((NLEV, len(lats), len(lons)), dtype=np.float64)
            data['aer_base_'+s] = zeros
            data['aer_warm_'+s] = zeros
    nc_bp.close(); nc_bs.close(); nc_ap.close(); nc_as.close()

    # Detect inactive perturbations (Step 1: aerosol only).
    # If all 6 species are zero in BOTH base and warm states, skip aerosol
    # rad_driver calls to save ~1+nspecies calls per column.
    aer_max = 0.0
    for s in AEROSOL_MAP:
        aer_max = max(aer_max,
                      float(np.nanmax(np.abs(data['aer_base_'+s]))),
                      float(np.nanmax(np.abs(data['aer_warm_'+s]))))
    data['skip_aerosol'] = aer_max < 1e-15
    print("Aerosol detection: max|mass_mixing_ratio|=%.3e -> skip_aerosol=%s" %
          (aer_max, data['skip_aerosol']))

    # Load lookup tables
    lut = {}
    for _, (lw, sw, _) in AEROSOL_MAP.items():
        for f in [lw, sw]:
            if f not in lut:
                lut[f] = load_lookup_table(os.path.join(LOOKUP_DIR, f))
    data['lut'] = lut

    # Load non-radiative forcing (optional)
    nonrad_path = cfg['input'].get('nonrad_forcing')
    if nonrad_path and os.path.exists(nonrad_path):
        print("Loading non-radiative forcing...")
        nc_pf = Dataset(nonrad_path)
        for nonrad_term in ['lhflx', 'shflx']:
            if nonrad_term not in nc_pf.variables:
                continue
            # lhflx/shflx are SURFACE energy budget terms; only the surface row
            # of the forcing vector should carry the value. Input NC packs the
            # surface flux at lev[0] (closest-to-surface in NC's surface→TOA
            # ordering); atm levels lev[1..N-1] are zero. We extract the lev[0]
            # 2D slice and place it exclusively at the surface row of frc_full.
            #
            # Do NOT copy the full 17-level NC profile into frc_full[:NLEV] —
            # that put -7.23 W/m² simultaneously at frc_full[NLEV-1] (atm
            # bottom, plev=1000hPa) and frc_full[NLEV] (surface, plev=1013hPa).
            # The Planck inverse propagated both → dT_lhflx and dT_shflx ~2×
            # what OLD CFRAM (xiaominghu) produces. Verified single-col probe
            # at lat=-0.47, lon=180:
            #   OLD: frc[sfc]=-7.55 → dT_lhflx[lev=1000]=-2.79 K
            #   Fu (with bug): same input → dT_lhflx[lev=1000]=-5.42 K (1.94×)
            sfc_2d = np.array(nc_pf.variables[nonrad_term][0, 0, :, :], dtype=np.float64)
            if sfc_2d.shape[0] != nlat or sfc_2d.shape[1] != nlon:
                print("  WARNING: %s shape %s != grid (%d,%d), skipping" % (
                    nonrad_term, sfc_2d.shape, nlat, nlon))
                continue
            frc_full = np.zeros((NLEV + 1, nlat, nlon))
            frc_full[NLEV, :, :] = sfc_2d
            frc_full = np.where(np.abs(frc_full) > 900, 0.0, frc_full)
            data['frc_' + nonrad_term] = frc_full
            print("  frc_%s: sfc mean=%.3f W/m2" % (nonrad_term, np.nanmean(frc_full[NLEV])))
        # Also load per-species aerosol forcing from same file
        aer_species_terms = ['bc', 'oc', 'sulf', 'seas', 'dust']
        for aer_term in aer_species_terms:
            if aer_term in nc_pf.variables:
                frc_3d = np.array(nc_pf.variables[aer_term][0, ::-1, :, :], dtype=np.float64)
                if frc_3d.shape[1] != nlat or frc_3d.shape[2] != nlon:
                    print("  WARNING: %s shape %s != grid (%d,%d), skipping" % (
                        aer_term, frc_3d.shape, nlat, nlon))
                    continue
                frc_full = np.zeros((NLEV + 1, nlat, nlon))
                frc_full[:NLEV, :, :] = frc_3d[:NLEV, :, :]
                frc_full[NLEV, :, :] = np.array(nc_pf.variables[aer_term][0, 0, :, :], dtype=np.float64)
                frc_full = np.where(np.abs(frc_full) > 900, 0.0, frc_full)
                data['frc_' + aer_term] = frc_full
                print("  frc_%s: sfc mean=%.3f W/m2" % (aer_term, np.nanmean(frc_full[NLEV])))
        nc_pf.close()
    else:
        print("No non-radiative forcing provided (radiative decomposition only)")

    n_cells = nlat * nlon
    print("Grid: %d x %d = %d points" % (nlat, nlon, n_cells))

    # Build task list
    tasks = [(i, j) for i in range(nlat) for j in range(nlon)]

    # Auto-cap workers at the grid size: spawning 200 workers for a 1×1
    # idealized RCE column wastes init time and process slots. For the
    # single-cell case we skip multiprocessing.Pool entirely and call
    # process_column synchronously — useful for climlab validation.
    nproc_eff = max(1, min(nproc, n_cells))
    if nproc_eff != nproc:
        print("Auto-capping nproc: requested %d → using %d (grid has only %d cells)"
              % (nproc, nproc_eff, n_cells))
    args.nproc = nproc_eff

    # Run
    import time
    t0 = time.time()

    # Result arrays
    terms = ['co2', 'q', 'ts', 'o3', 'solar', 'albedo', 'cloud', 'aerosol', 'warm',
             'cloud_lw', 'cloud_sw', 'full']
    nonrad_terms = ['lhflx', 'shflx']
    # Naming map (pyCFRAM ↔ OLD CFRAM Cai-Tung Fortran):
    #   atmdyn ≡ OLD dt_atm_dyn  (atm rows -ΔR_full, sfc=0)
    #   sfcdyn  = pyCFRAM-specific (sfc row -ΔR_full only)
    #   ocndyn ≡ OLD dt_sfc_dyn  (sfc row -ΔR_full - Δlh - Δsh)
    #   dry    ≡ OLD dt_dyn      (full-column drdt⁻¹·frc_full = atmdyn+sfcdyn)
    dyn_terms = ['atmdyn', 'sfcdyn', 'ocndyn', 'dry']
    aer_species_terms = ['bc', 'oc', 'sulf', 'seas', 'dust']           # path-B names
    fortran_aer_species = SPECIES_ORDER                                # ['bc','ocphi','ocpho','sulf','ss','dust']
    # Union of all species names (some overlap: bc/sulf/dust)
    all_species = list(dict.fromkeys(aer_species_terms + fortran_aer_species))
    all_dT_terms = terms + dyn_terms + nonrad_terms + all_species
    dT_out = {t: np.full((NLEV+1, nlat, nlon), np.nan) for t in all_dT_terms}
    # Bulk radiative + per-species forcings stored in NetCDF
    frc_out = {t: np.full((NLEV+1, nlat, nlon), np.nan) for t in terms + fortran_aer_species}

    def _absorb(ilat, ilon, result):
        for t in terms:
            dT_out[t][:, ilat, ilon] = result['dT_'+t]
            frc_out[t][:, ilat, ilon] = result['frc_'+t]
        for t in dyn_terms + nonrad_terms + all_species:
            if 'dT_'+t in result:
                dT_out[t][:, ilat, ilon] = result['dT_'+t]
        for t in fortran_aer_species:
            if 'frc_'+t in result:
                frc_out[t][:, ilat, ilon] = result['frc_'+t]

    done = 0
    if n_cells == 1:
        # Single-cell short-circuit: avoid mp.Pool spin-up. Useful for climlab
        # 1×1 RCE validation cases.
        print("Single-cell synchronous run (no multiprocessing)")
        init_worker(data)
        ilat, ilon, result = process_column(tasks[0])
        _absorb(ilat, ilon, result)
        done = 1
    else:
        print("Starting %d workers..." % args.nproc)
        with Pool(args.nproc, initializer=init_worker, initargs=(data,)) as pool:
            for ilat, ilon, result in pool.imap_unordered(process_column, tasks, chunksize=4):
                _absorb(ilat, ilon, result)
                done += 1
                if done % 200 == 0:
                    elapsed = time.time() - t0
                    rate = done / elapsed
                    eta = (len(tasks) - done) / rate
                    print("  %d/%d done (%.1f pts/s, ETA %.0fs)" % (done, len(tasks), rate, eta))

    elapsed = time.time() - t0
    print("Completed in %.1f seconds (%.1f pts/s)" % (elapsed, len(tasks) / elapsed))

    # Save as NetCDF
    outfile = os.path.join(OUTDIR, 'cfram_result.nc')
    nc = Dataset(outfile, 'w')
    nc.createDimension('lev', NLEV + 1)
    nc.createDimension('lat', nlat)
    nc.createDimension('lon', nlon)
    nc.createVariable('lat', 'f8', ('lat',))[:] = lats
    nc.createVariable('lon', 'f8', ('lon',))[:] = lons
    # NetCDF lev convention: surface→TOA + ground (1013) at end, mirroring
    # the input *_pres.nc convention. Internal arrays (dT_out/frc_out) are in
    # TOA→surface + surface row order — atm rows must be reversed before save
    # so they align with lev. (Surface row stays at index NLEV → lev=1013.)
    levs = np.concatenate([FORTRAN_PLEV[::-1], [FORTRAN_PLEV[-1] + 13]])
    nc.createVariable('lev', 'f8', ('lev',))[:] = levs

    def reorder_for_save(arr):
        """TOA→sfc atm + sfc → sfc→TOA atm + ground (matches lev order)."""
        return np.concatenate([arr[:NLEV][::-1], arr[NLEV:NLEV + 1]], axis=0)

    # Apply output_terms filter: if set, restrict which radiative dT/frc are
    # written (always keep dyn / nonrad / lhflx / shflx — those are physically
    # required for closure). Idealized cases (clear-sky climlab RCE) use this
    # to drop always-zero `cloud / aerosol / o3 / albedo / solar` rows.
    REQUIRED_DERIVED = set(dyn_terms + nonrad_terms)   # never filtered
    if output_terms is not None:
        allowed_rad = set(output_terms)
    else:
        allowed_rad = None   # write everything

    def _allowed(t):
        if t in REQUIRED_DERIVED:
            return True
        if allowed_rad is None:
            return True
        return t in allowed_rad

    for t in all_dT_terms:
        if not _allowed(t):
            continue
        v = nc.createVariable('dT_'+t, 'f8', ('lev', 'lat', 'lon'))
        v[:] = reorder_for_save(dT_out[t])
    for t in terms + fortran_aer_species:
        if not _allowed(t):
            continue
        if t in frc_out:
            v = nc.createVariable('frc_'+t, 'f8', ('lev', 'lat', 'lon'))
            v[:] = reorder_for_save(frc_out[t])

    # Compute and save observed dT
    nc_bp2 = Dataset(cfg['input']['base_pres'])
    nc_ap2 = Dataset(cfg['input']['perturbed_pres'])
    nc_bs2 = Dataset(cfg['input']['base_surf'])
    nc_as2 = Dataset(cfg['input']['perturbed_surf'])
    dT_obs = np.zeros((NLEV + 1, nlat, nlon))
    dT_obs[:NLEV] = np.array(nc_ap2.variables['ta_lay'][0, ::-1, :, :], dtype=np.float64) - \
                     np.array(nc_bp2.variables['ta_lay'][0, ::-1, :, :], dtype=np.float64)
    dT_obs[NLEV] = np.array(nc_as2.variables['ts'][0, :, :], dtype=np.float64) - \
                    np.array(nc_bs2.variables['ts'][0, :, :], dtype=np.float64)
    nc_bp2.close(); nc_ap2.close(); nc_bs2.close(); nc_as2.close()

    v = nc.createVariable('dT_observed', 'f8', ('lev', 'lat', 'lon'))
    v[:] = reorder_for_save(dT_obs)

    nc.close()
    print("Saved: %s" % outfile)


if __name__ == '__main__':
    main()
