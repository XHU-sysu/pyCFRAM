#!/usr/bin/env python3
"""Build standard CFRAM input NetCDF files from reanalysis data.

Reads case.yaml 'source' configuration and generates:
    cases/<case>/input/base_pres.nc
    cases/<case>/input/base_surf.nc
    cases/<case>/input/perturbed_pres.nc
    cases/<case>/input/perturbed_surf.nc

Usage:
    python scripts/build_case_input.py --case eh22
    python scripts/build_case_input.py --case eh22 --dry-run
"""

import os, sys, argparse, importlib
import numpy as np
from netCDF4 import Dataset

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.config import load_case

# Registry: source.type -> module that registers it via @register_source
# (data/source_base.py). Importing the module triggers the decorator, which
# is all get_source() needs. Extend this dict when adding a new registered
# DataSource plugin (e.g. a new data/<x>_source.py). If a type is missing
# here, we simply skip the dynamic import (its module may already have been
# imported elsewhere) and let get_source() raise its own clear
# "Unknown source type" error — this stays a small explicit map, not a
# generic plugin-discovery mechanism (docs/plan_ph3.md §2.3-1).
SOURCE_MODULES = {
    'era5_daily': 'data.era5_source',
    'era5_date_range': 'data.era5_source',
    'era5_merra2': 'data.era5_source',
    'cmip6_damip': 'data.cmip6_damip_source',
}

# Variable classification: which go into pres vs surf files
PRES_3D_VARS = ['ta_lay', 'q', 'o3', 'camt', 'cliq', 'cice', 'co2',
                'bc', 'ocphi', 'ocpho', 'sulf', 'ss', 'dust']
SURF_2D_VARS = ['ts', 'ps', 'solar', 'albedo']

# Optional surface variables: written to *_surf.nc only when the source's
# state dict actually provides them (state[varname] is not None). NEVER
# fabricate a default/zero value here for a variable a source doesn't
# produce -- e.g. writing huss=0 for ERA5 (which has no huss) would silently
# change ERA5/Fu-engine behavior, since the runner's HOLD-fallback logic
# treats "variable absent" differently from "variable present but ~0"
# (docs/plan_ph3.md §2.3-2). Contrast with SURF_2D_VARS above, whose missing
# case fills zeros with a warning -- that fallback predates this module's
# optional-variable support and is left as-is for those four required vars.
OPTIONAL_SURF_2D_VARS = ['huss']

# Units metadata
VAR_UNITS = {
    'ta_lay': 'K', 'q': 'kg/kg', 'o3': 'kg/kg',
    'camt': '1', 'cliq': 'kg/kg', 'cice': 'kg/kg',
    'co2': 'mol/mol',
    'bc': 'kg/kg', 'ocphi': 'kg/kg', 'ocpho': 'kg/kg',
    'sulf': 'kg/kg', 'ss': 'kg/kg', 'dust': 'kg/kg',
    'ts': 'K', 'ps': 'Pa', 'solar': 'W/m2', 'albedo': '1',
    'huss': 'kg/kg',
}

VAR_LONG_NAMES = {
    'ta_lay': 'Layer-mean temperature',
    'q': 'Specific humidity',
    'o3': 'Ozone mass mixing ratio',
    'camt': 'Cloud fraction',
    'cliq': 'Cloud liquid water content',
    'cice': 'Cloud ice water content',
    'co2': 'CO2 volume mixing ratio',
    'bc': 'Black carbon mixing ratio',
    'ocphi': 'Hydrophilic organic carbon',
    'ocpho': 'Hydrophobic organic carbon',
    'sulf': 'Sulfate mixing ratio',
    'ss': 'Sea salt mixing ratio',
    'dust': 'Dust mixing ratio',
    'ts': 'Skin temperature',
    'ps': 'Surface pressure',
    'solar': 'TOA incident solar radiation',
    'albedo': 'Surface albedo',
    'huss': 'Near-surface (2m) specific humidity',
}

AEROSOL_VARS = ['bc', 'ocphi', 'ocpho', 'sulf', 'ss', 'dust']
AEROSOL_MAX_KGKG = 1e-5


def validate_states(*states):
    """Reject corrupt state arrays before they are written or run by CFRAM."""
    for state_name, state in states:
        for varname in PRES_3D_VARS + SURF_2D_VARS:
            if varname not in state:
                continue
            data = np.asarray(state[varname])
            if not np.all(np.isfinite(data)):
                raise ValueError(f"{state_name}.{varname} contains non-finite values")
            if varname in AEROSOL_VARS:
                vmin = float(np.min(data))
                vmax = float(np.max(data))
                if vmin < 0.0 or vmax > AEROSOL_MAX_KGKG:
                    raise ValueError(
                        f"{state_name}.{varname} aerosol sanity check failed: "
                        f"range=[{vmin:.3e}, {vmax:.3e}] kg/kg"
                    )


def write_pres_nc(filepath, state):
    """Write pressure-level variables to NetCDF."""
    lat = state['lat']
    lon = state['lon']
    lev = state['lev']

    nc = Dataset(filepath, 'w', format='NETCDF4')
    nc.createDimension('time', 1)
    nc.createDimension('lev', len(lev))
    nc.createDimension('lat', len(lat))
    nc.createDimension('lon', len(lon))

    # Coordinate variables
    v = nc.createVariable('time', 'f8', ('time',))
    v[:] = [0]
    v.units = 'days since 2000-01-01'

    v = nc.createVariable('lev', 'f8', ('lev',))
    v[:] = lev
    v.units = 'hPa'
    v.long_name = 'Pressure level'
    v.positive = 'down'

    v = nc.createVariable('lat', 'f8', ('lat',))
    v[:] = lat
    v.units = 'degrees_north'

    v = nc.createVariable('lon', 'f8', ('lon',))
    v[:] = lon
    v.units = 'degrees_east'

    # Data variables
    for varname in PRES_3D_VARS:
        if varname not in state:
            print(f"  Warning: {varname} not in state, filling with zeros")
            data = np.zeros((len(lev), len(lat), len(lon)), dtype=np.float64)
        else:
            data = state[varname]
        v = nc.createVariable(varname, 'f8', ('time', 'lev', 'lat', 'lon'),
                              zlib=True, complevel=4)
        v[0, :, :, :] = data
        v.units = VAR_UNITS.get(varname, '')
        v.long_name = VAR_LONG_NAMES.get(varname, varname)

    nc.close()
    fsize = os.path.getsize(filepath) / 1e6
    print(f"  Wrote {filepath} ({fsize:.1f} MB)")


def write_surf_nc(filepath, state):
    """Write surface variables to NetCDF."""
    lat = state['lat']
    lon = state['lon']

    nc = Dataset(filepath, 'w', format='NETCDF4')
    nc.createDimension('time', 1)
    nc.createDimension('lat', len(lat))
    nc.createDimension('lon', len(lon))

    v = nc.createVariable('time', 'f8', ('time',))
    v[:] = [0]
    v.units = 'days since 2000-01-01'

    v = nc.createVariable('lat', 'f8', ('lat',))
    v[:] = lat
    v.units = 'degrees_north'

    v = nc.createVariable('lon', 'f8', ('lon',))
    v[:] = lon
    v.units = 'degrees_east'

    for varname in SURF_2D_VARS:
        if varname not in state:
            print(f"  Warning: {varname} not in state, filling with zeros")
            data = np.zeros((len(lat), len(lon)), dtype=np.float64)
        else:
            data = state[varname]
        v = nc.createVariable(varname, 'f8', ('time', 'lat', 'lon'),
                              zlib=True, complevel=4)
        v[0, :, :] = data
        v.units = VAR_UNITS.get(varname, '')
        v.long_name = VAR_LONG_NAMES.get(varname, varname)

    # Optional variables: write only if the source actually supplied them.
    # No zero-fill fallback here -- see OPTIONAL_SURF_2D_VARS docstring above.
    for varname in OPTIONAL_SURF_2D_VARS:
        data = state.get(varname)
        if data is None:
            continue
        v = nc.createVariable(varname, 'f8', ('time', 'lat', 'lon'),
                              zlib=True, complevel=4)
        v[0, :, :] = data
        v.units = VAR_UNITS.get(varname, '')
        v.long_name = VAR_LONG_NAMES.get(varname, varname)

    nc.close()
    fsize = os.path.getsize(filepath) / 1e6
    print(f"  Wrote {filepath} ({fsize:.1f} MB)")


def write_nonrad_nc(filepath, state, nonrad):
    """Write non-radiative forcing to NetCDF.

    Format matches paper_data partial_forcing.nc convention:
      dims: (time=1, lev=38, lat, lon)
      lev[0]=1013 (surface), lev[1:37]=1000→1 (atm levels, surface→TOA)
      Only surface level (lev[0]) has real values; atm levels are fill (-999).
    """
    lat = state['lat']
    lon = state['lon']
    lev = state['lev']  # surface→TOA in hPa
    # Build the 38-level axis: [1013, lev[0], lev[1], ..., lev[-1]]
    lev_out = np.concatenate([[1013.0], lev])
    nlev_out = len(lev_out)

    FILL = -999.0

    nc = Dataset(filepath, 'w', format='NETCDF4')
    nc.createDimension('time', 1)
    nc.createDimension('lev', nlev_out)
    nc.createDimension('lat', len(lat))
    nc.createDimension('lon', len(lon))

    v = nc.createVariable('time', 'f8', ('time',)); v[:] = [0]
    v = nc.createVariable('lev', 'f8', ('lev',)); v[:] = lev_out
    v.units = 'hPa'
    v = nc.createVariable('lat', 'f8', ('lat',)); v[:] = lat
    v.units = 'degrees_north'
    v = nc.createVariable('lon', 'f8', ('lon',)); v[:] = lon
    v.units = 'degrees_east'

    for varname in ('lhflx', 'shflx'):
        if varname not in nonrad:
            continue
        data_4d = np.full((1, nlev_out, len(lat), len(lon)), FILL, dtype=np.float64)
        # Surface value at lev[0]=1013
        data_4d[0, 0, :, :] = nonrad[varname]
        v = nc.createVariable(
            varname, 'f8', ('time', 'lev', 'lat', 'lon'),
            zlib=True, complevel=4, fill_value=FILL,
        )
        v[:] = data_4d
        v.missing_value = FILL
        v.units = 'W/m2'
        print(f"  {varname}: sfc mean={nonrad[varname].mean():.3f} W/m2")

    nc.close()
    fsize = os.path.getsize(filepath) / 1e6
    print(f"  Wrote {filepath} ({fsize:.1f} MB)")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--case', required=True, help='Case name (e.g. eh22)')
    parser.add_argument('--dry-run', action='store_true',
                        help='Only show warm period detection, do not write files')
    args = parser.parse_args()

    cfg = load_case(args.case)

    if 'source' not in cfg:
        print(f"ERROR: No 'source' section in cases/{args.case}/case.yaml")
        print("Add a source configuration. See docs/input_spec.md.")
        sys.exit(1)

    print(f"=== Building CFRAM input for case: {args.case} ===")
    src_type = cfg['source']['type']
    print(f"Source type: {src_type}")

    # Registry-driven import: trigger the @register_source decorator for
    # whichever plugin module owns this source.type (see SOURCE_MODULES
    # above). If src_type isn't in the registry, we don't hard-fail here --
    # its module may already be imported by something else -- and let
    # get_source() below raise its own descriptive "Unknown source type"
    # error if it truly never got registered.
    module_name = SOURCE_MODULES.get(src_type)
    if module_name:
        importlib.import_module(module_name)
    else:
        print(f"  Note: no SOURCE_MODULES entry for '{src_type}'; "
              f"assuming its module is already imported/registered.")

    from data.source_base import get_source
    source = get_source(cfg)

    if args.dry_run:
        print("\n[DRY RUN] Would build states but not write files.")
        # Still run build_states to show warm period info
        base, pert, nonrad = source.build_states()
        print("\nDry run complete. No files written.")
        return

    base_state, pert_state, nonrad = source.build_states()
    validate_states(('base', base_state), ('perturbed', pert_state))

    # Write output files — remove any existing symlinks first
    input_dir = os.path.join(cfg['_case_dir'], 'input')
    os.makedirs(input_dir, exist_ok=True)
    for fn in ['base_pres.nc', 'base_surf.nc', 'perturbed_pres.nc',
               'perturbed_surf.nc', 'nonrad_forcing.nc']:
        fpath = os.path.join(input_dir, fn)
        if os.path.islink(fpath):
            os.remove(fpath)
            print(f"  Removed symlink: {fpath}")

    print(f"\nWriting NetCDF files to {input_dir}/")
    write_pres_nc(os.path.join(input_dir, 'base_pres.nc'), base_state)
    write_surf_nc(os.path.join(input_dir, 'base_surf.nc'), base_state)
    write_pres_nc(os.path.join(input_dir, 'perturbed_pres.nc'), pert_state)
    write_surf_nc(os.path.join(input_dir, 'perturbed_surf.nc'), pert_state)

    # Write non-radiative forcing if available
    if nonrad:
        write_nonrad_nc(os.path.join(input_dir, 'nonrad_forcing.nc'),
                        base_state, nonrad)

    print(f"\n=== Done. Input files ready in cases/{args.case}/input/ ===")


if __name__ == '__main__':
    main()
