"""Tests for WP-M4.3 (docs/plan_ph3.md §2.3/§2.4/§8): generalizing the
build-step dispatch in ``run_case.py`` and the two ``scripts/build_case_input.py``
hard points (registry-driven source import, optional-``huss`` surf writer).

This is the highest-regression-risk WP in Phase 3: both files are shared by
every existing case (ERA5 eh13/eh22/india_wb23 + the bespoke CESM2 4xCO2
path), not just the new DAMIP path. The tests below are split into three
groups matching the WP's two non-negotiable regression guardrails plus the
new-path check:

1. ``run_case.py`` build dispatch (in-process, no subprocess/no real case
   dirs): the ``cesm2_cmip6`` branch must still fire its exact 3-script
   sequence; any other non-None source.type (including brand new ones, e.g.
   ``cmip6_damip``) must route through ``build_case_input.py``; an absent/
   None source.type must still skip the build step entirely.
2. ``scripts/build_case_input.py`` write_surf_nc: an ERA5-shaped state dict
   (no ``huss`` key) must produce a surf NC with NO ``huss`` variable at all
   -- the core "don't silently change ERA5 behavior" guardrail. A
   DAMIP-shaped state dict (with ``huss``) must get it written.
3. The ``SOURCE_MODULES`` registry dynamic-import fix: importing via the
   registry actually registers the source types this WP cares about.
"""
import importlib
import os
import sys

import numpy as np
import pytest
from netCDF4 import Dataset

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_ROOT)

import run_case  # noqa: E402
from scripts import build_case_input as bci  # noqa: E402


# ---------------------------------------------------------------------------
# 1. run_case.py build dispatch
# ---------------------------------------------------------------------------

def _run_build_dispatch(monkeypatch, cfg, argv_case='dispatch_test'):
    """Drive run_case.main() for just the 'build' step, recording every
    run_step(script, args_list) call instead of actually spawning a
    subprocess. Returns the list of (script, args_list) calls made.
    """
    calls = []

    def fake_run_step(script, args_list):
        calls.append((script, list(args_list)))

    monkeypatch.setattr(run_case, 'load_case', lambda case_name: cfg)
    monkeypatch.setattr(run_case, 'run_step', fake_run_step)
    monkeypatch.setattr(sys, 'argv', ['run_case.py', argv_case, '--step', 'build'])

    run_case.main()
    return calls


def test_dispatch_cesm2_cmip6_unchanged(monkeypatch):
    """The bespoke path-B 3-script sequence must fire exactly as before,
    in order, each with ['--case', <name>] -- untouched by this WP.
    """
    cfg = {'source': {'type': 'cesm2_cmip6'}, 'input': {}}
    calls = _run_build_dispatch(monkeypatch, cfg, argv_case='cesm2_4xco2_official')

    assert calls == [
        ('build_cesm2_official.py', ['--case', 'cesm2_4xco2_official']),
        ('inject_cesm_o3.py', ['--case', 'cesm2_4xco2_official']),
        ('mask_subsurface_layers.py', ['--case', 'cesm2_4xco2_official']),
    ]


@pytest.mark.parametrize('src_type', ['era5_daily', 'era5_date_range', 'era5_merra2'])
def test_dispatch_era5_variants_go_through_generic_writer(monkeypatch, src_type):
    cfg = {'source': {'type': src_type}, 'input': {}}
    calls = _run_build_dispatch(monkeypatch, cfg, argv_case='eh13')
    assert calls == [('build_case_input.py', ['--case', 'eh13'])]


def test_dispatch_cmip6_damip_goes_through_generic_writer(monkeypatch):
    """The new source type is NOT in any hardcoded whitelist -- this is the
    actual generalization under test (docs/plan_ph3.md §2.4).
    """
    cfg = {'source': {'type': 'cmip6_damip'}, 'input': {}}
    calls = _run_build_dispatch(monkeypatch, cfg, argv_case='damip_ipsl_histaer')
    assert calls == [('build_case_input.py', ['--case', 'damip_ipsl_histaer'])]


def test_dispatch_hypothetical_future_source_also_generic(monkeypatch):
    """Any non-cesm2_cmip6, non-None source.type dispatches generically --
    even one this WP has never heard of. build_case_input.py itself (not
    run_case.py) is responsible for raising "Unknown source type" if it's
    never actually registered.
    """
    cfg = {'source': {'type': 'some_future_source'}, 'input': {}}
    calls = _run_build_dispatch(monkeypatch, cfg, argv_case='future_case')
    assert calls == [('build_case_input.py', ['--case', 'future_case'])]


def test_dispatch_no_source_type_skips_build(monkeypatch, capsys):
    cfg = {'input': {}}
    calls = _run_build_dispatch(monkeypatch, cfg, argv_case='preassembled_case')
    assert calls == []
    out = capsys.readouterr().out
    assert 'skipping build step' in out


def test_dispatch_explicit_none_source_type_skips_build(monkeypatch):
    cfg = {'source': {'type': None}, 'input': {}}
    calls = _run_build_dispatch(monkeypatch, cfg, argv_case='null_type_case')
    assert calls == []


# ---------------------------------------------------------------------------
# 2. write_surf_nc optional-huss support (the ERA5 non-regression guardrail)
# ---------------------------------------------------------------------------

def _minimal_surf_state(nlat=3, nlon=4, extra=None):
    state = {
        'lat': np.linspace(-60.0, 60.0, nlat),
        'lon': np.linspace(0.0, 270.0, nlon),
        'ts': np.full((nlat, nlon), 288.0),
        'ps': np.full((nlat, nlon), 101300.0),
        'solar': np.full((nlat, nlon), 340.0),
        'albedo': np.full((nlat, nlon), 0.15),
    }
    if extra:
        state.update(extra)
    return state


def test_write_surf_nc_era5_shaped_state_has_no_huss(tmp_path):
    """ERA5's build_states() never sets a 'huss' key. The surf NC produced
    for such a state must NOT contain a huss variable at all -- writing
    huss=0 would silently change existing ERA5/Fu-engine behavior (a 0
    value doesn't trip the runner's |x|>900-missing sentinel the way an
    absent variable does).
    """
    state = _minimal_surf_state()
    out_path = os.path.join(str(tmp_path), 'era5_shaped_surf.nc')
    bci.write_surf_nc(out_path, state)

    with Dataset(out_path) as nc:
        varnames = set(nc.variables.keys())
        assert 'huss' not in varnames
        assert {'ts', 'ps', 'solar', 'albedo'} <= varnames


def test_write_surf_nc_damip_shaped_state_writes_huss(tmp_path):
    """A DAMIP-shaped state dict always carries 'huss' -- it must be
    written, with correct values/units/long_name.
    """
    huss_vals = np.full((3, 4), 0.0065)
    state = _minimal_surf_state(extra={'huss': huss_vals})
    out_path = os.path.join(str(tmp_path), 'damip_shaped_surf.nc')
    bci.write_surf_nc(out_path, state)

    with Dataset(out_path) as nc:
        assert 'huss' in nc.variables
        np.testing.assert_allclose(nc.variables['huss'][0, :, :], huss_vals)
        assert nc.variables['huss'].units == 'kg/kg'


def test_write_surf_nc_explicit_none_huss_is_not_written(tmp_path):
    """state.get('huss') is None (key present but None) must be treated the
    same as "absent" -- no phantom variable.
    """
    state = _minimal_surf_state(extra={'huss': None})
    out_path = os.path.join(str(tmp_path), 'none_huss_surf.nc')
    bci.write_surf_nc(out_path, state)

    with Dataset(out_path) as nc:
        assert 'huss' not in nc.variables


def test_surf_2d_vars_and_optional_list_unchanged_shape():
    """Pin down the exact split this WP introduces: the four required surf
    vars are untouched, huss is the only optional one added.
    """
    assert bci.SURF_2D_VARS == ['ts', 'ps', 'solar', 'albedo']
    assert bci.OPTIONAL_SURF_2D_VARS == ['huss']


def test_validate_states_unaffected_by_optional_huss_key():
    """validate_states() only inspects PRES_3D_VARS + SURF_2D_VARS, so an
    extra 'huss' key alongside otherwise-valid data must not raise --
    confirms no logic changes were needed there (docs/plan_ph3.md §8
    WP-M4.3 Done criteria item 4).
    """
    good_state = _minimal_surf_state(extra={'huss': np.full((3, 4), 0.007)})
    # validate_states() also inspects PRES_3D_VARS but tolerates their
    # absence (continue on missing varname), so a surf-only state is valid
    # input here.
    bci.validate_states(('base', good_state))


# ---------------------------------------------------------------------------
# 3. SOURCE_MODULES registry-driven import
# ---------------------------------------------------------------------------

def test_source_modules_registry_covers_known_types():
    assert bci.SOURCE_MODULES['era5_daily'] == 'data.era5_source'
    assert bci.SOURCE_MODULES['era5_date_range'] == 'data.era5_source'
    assert bci.SOURCE_MODULES['cmip6_damip'] == 'data.cmip6_damip_source'


def test_registry_import_actually_registers_sources():
    """Exercises the exact mechanism build_case_input.main() now uses:
    importlib.import_module(SOURCE_MODULES[src_type]) must make get_source()
    resolve that type, for every type in the registry.
    """
    from data.source_base import get_source

    for src_type, module_name in bci.SOURCE_MODULES.items():
        importlib.import_module(module_name)

    # era5_daily / era5_date_range / cmip6_damip must all now resolve.
    for src_type in ('era5_daily', 'era5_date_range', 'cmip6_damip'):
        source = get_source({'source': {'type': src_type}})
        assert source.__class__.__name__  # instantiated without error
