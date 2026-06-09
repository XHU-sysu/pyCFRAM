#!/usr/bin/env python3
"""pyCFRAM entry point: run CFRAM decomposition for a case.

Usage:
    python3 run_case.py eh13                    # full workflow (build → run → plot)
    python3 run_case.py eh13 --step build       # only build ERA5+MERRA-2 input
    python3 run_case.py eh13 --step run         # only run CFRAM
    python3 run_case.py eh13 --step plot        # only plot results
    python3 run_case.py eh13 --nproc 20         # override parallel workers
"""
import os
import sys
import argparse
import subprocess

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from core.config import load_case, get_nproc, PROJECT_ROOT


def run_step(script, args_list):
    """Run a script with arguments."""
    # -u: unbuffered stdout/stderr so progress lines flush in real time when
    # output is redirected to a file (e.g. via nohup).
    cmd = [sys.executable, '-u', os.path.join(PROJECT_ROOT, 'scripts', script)] + args_list
    print("\n>>> %s" % ' '.join(cmd))
    ret = subprocess.call(cmd)
    if ret != 0:
        print("ERROR: %s failed with exit code %d" % (script, ret))
        sys.exit(ret)


def main():
    parser = argparse.ArgumentParser(description='pyCFRAM: run CFRAM decomposition')
    parser.add_argument('case', help='Case name (directory under cases/)')
    parser.add_argument('--step', choices=['build', 'run', 'plot', 'all'],
                        default='all', help='Which step to run (default: all)')
    parser.add_argument('--nproc', type=int, default=None,
                        help='Number of parallel workers (default: from config or all CPUs)')
    args = parser.parse_args()

    # Validate case exists
    cfg = load_case(args.case)
    nproc = args.nproc or get_nproc(cfg)
    print("=" * 60)
    print("pyCFRAM: %s (%s)" % (cfg.get('case_name', args.case), cfg.get('description', '')))
    print("Workers: %d" % nproc)
    print("=" * 60)

    steps = ['build', 'run', 'plot'] if args.step == 'all' else [args.step]

    if 'build' in steps:
        # Dispatch on source.type: ERA5 reanalysis → build_case_input.py;
        # cesm2_cmip6 → build_cesm2_official.py; absent → skip (input pre-supplied).
        src_type = cfg.get('source', {}).get('type')
        if src_type == 'cesm2_cmip6':
            run_step('build_cesm2_official.py', ['--case', args.case])
            # Post-build: O3 inject + subsurface mask
            run_step('inject_cesm_o3.py',           ['--case', args.case])
            run_step('mask_subsurface_layers.py',   ['--case', args.case])
        elif src_type in ('era5_daily', 'era5_date_range', 'era5_merra2', None) and 'source' in cfg:
            run_step('build_case_input.py', ['--case', args.case])
        else:
            print('No source block in case.yaml — skipping build step (input must be pre-supplied)')

    if 'run' in steps:
        for key in ['base_pres', 'base_surf', 'perturbed_pres', 'perturbed_surf']:
            path = cfg['input'].get(key)
            if not path or not os.path.exists(path):
                print("ERROR: Input file missing: %s" % path)
                print("Run: python3 scripts/build_case_input.py --case %s" % args.case)
                sys.exit(1)
        run_step('run_parallel_python.py', ['--case', args.case, '--nproc', str(nproc)])

    if 'plot' in steps:
        run_step('plot_fig3_independent.py', [args.case])

    print("\n" + "=" * 60)
    print("Done! Results in: cases/%s/output/ and cases/%s/figures/" % (args.case, args.case))
    print("=" * 60)


if __name__ == '__main__':
    main()
