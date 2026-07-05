#!/usr/bin/env python3
"""Download CMIP6 DAMIP single-forcing experiment data from ESGF.

Dry-run mode shows file list and total volume without downloading.
Actual download uses multi-node fallback and resume support.

Usage:
    python3 scripts/download_damip.py --model IPSL-CM6A-LR --experiment hist-aer --dry-run
    python3 scripts/download_damip.py --model IPSL-CM6A-LR --experiment hist-aer
    python3 scripts/download_damip.py --all-m5 --dry-run
"""

import sys
import argparse
import json
from pathlib import Path
from typing import List, Tuple

# Import ESGF client (pure stdlib)
sys.path.insert(0, str(Path(__file__).parent.parent))
from data.esgf_fetch import search_with_fallback, list_files, filename_time_overlap, fetch


# M4 and M5 model lists (per plan §1.4 and Appendix D)
M4_MODELS = [
    'IPSL-CM6A-LR',
    'MRI-ESM2-0',
    'CESM2',
]

M5_MODELS = [
    'IPSL-CM6A-LR',
    'MRI-ESM2-0',
    'CNRM-CM6-1',
    'MIROC6',
    'GISS-E2-1-G',
    'HadGEM3-GC31-LL',
    'CanESM5',
    'CESM2',
]

# Core variables (§7.3)
CORE_VARIABLES = [
    'ta', 'hus', 'ts', 'ps',
    'cl', 'clw', 'cli',
    'rsdt', 'rsds', 'rsus',
    'hfls', 'hfss', 'huss',
]


def get_variables_for_model(model: str, experiment: str,
                           nodes: List[str]) -> List[str]:
    """Determine which variables to fetch for a model.

    hist-aer always fetches core variables. o3 is conditional (only if available).
    """
    variables = CORE_VARIABLES.copy()

    # Check if model has o3
    for node in nodes:
        try:
            datasets, _ = search_with_fallback(
                model=model,
                experiment=experiment,
                variables=['o3'],
                nodes=[node],
            )
            if datasets:
                variables.append('o3')
                break
        except RuntimeError:
            pass

    return variables


def search_damip_files(model: str, experiment: str, base_years: Tuple[int, int],
                       warm_years: Tuple[int, int], variant: str = None,
                       grid: str = None, nodes: List[str] = None) -> List[dict]:
    """Search for DAMIP files, filtering by time overlap.

    Args:
        model: CMIP6 source_id
        experiment: experiment_id (e.g. 'hist-aer')
        base_years: (start, end) for base climatology
        warm_years: (start, end) for warm climatology
        variant: optional variant_label filter
        grid: optional grid_label filter
        nodes: list of ESGF nodes to try

    Returns:
        List of file dicts ready for download
    """
    if nodes is None:
        nodes = [
            'https://esgf-data.dkrz.de/esg-search/search',
            'https://esgf.ceda.ac.uk/esg-search/search',
        ]

    variables = get_variables_for_model(model, experiment, nodes)

    # Search for datasets
    datasets, successful_node = search_with_fallback(
        model=model,
        experiment=experiment,
        variables=variables,
        variant=variant,
        grid=grid,
        nodes=nodes,
    )

    if not datasets:
        raise RuntimeError(f"No datasets found for {model} {experiment}")

    # For each dataset, list files and filter by time overlap
    all_files = []
    for dataset in datasets:
        dataset_id = dataset['id']
        try:
            files = list_files(dataset_id, node=successful_node)
        except RuntimeError as e:
            print(f"Warning: Could not list files for {dataset_id}: {e}")
            continue

        # Filter by time overlap
        for f in files:
            if filename_time_overlap(f['title'], base_years, warm_years):
                all_files.append(f)

    return all_files


def format_size(size_bytes: int) -> str:
    """Format bytes as human-readable size."""
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if size_bytes < 1024.0:
            return f"{size_bytes:.1f} {unit}"
        size_bytes /= 1024.0
    return f"{size_bytes:.1f} PB"


def dry_run_report(model: str, experiment: str, files: List[dict]) -> Tuple[int, float]:
    """Print dry-run report of files to be downloaded.

    Returns:
        (file_count, total_size_gb)
    """
    total_bytes = sum(f['size'] for f in files)
    total_gb = total_bytes / (1024**3)

    print(f"\n{'='*70}")
    print(f"DRY-RUN: {model} / {experiment}")
    print(f"{'='*70}")
    print(f"Found {len(files)} files, total {format_size(total_bytes)} ({total_gb:.2f} GB)")
    print()
    print(f"{'Variable':<6} {'File':<50} {'Size':>12}")
    print(f"{'-'*68}")

    by_variable = {}
    for f in files:
        # Extract variable from filename (first component before _)
        var = f['title'].split('_')[0]
        if var not in by_variable:
            by_variable[var] = []
        by_variable[var].append(f)

    for var in sorted(by_variable.keys()):
        var_files = by_variable[var]
        var_size = sum(f['size'] for f in var_files)
        for i, f in enumerate(var_files):
            var_str = var if i == 0 else ''
            print(f"{var_str:<6} {f['title']:<50} {format_size(f['size']):>12}")
        if len(var_files) > 1:
            print(f"{'':6} {'Subtotal':<44} {format_size(var_size):>18}")

    print(f"{'-'*68}")
    print(f"{'Total':<50} {format_size(total_bytes):>18}")
    print(f"{'='*70}\n")

    return len(files), total_gb


def download_files(model: str, experiment: str, files: List[dict],
                   output_dir: str = 'raw_data/cmip6_damip', verbose: bool = True):
    """Download list of files with resume support.

    Args:
        model: CMIP6 source_id (for organizing output dir)
        experiment: experiment_id (for organizing output dir)
        files: list of file dicts from search_damip_files
        output_dir: base directory for downloads
        verbose: print progress
    """
    dest_dir = Path(output_dir) / model / experiment
    dest_dir.mkdir(parents=True, exist_ok=True)

    manifest_path = dest_dir / 'manifest.json'
    existing = {}
    if manifest_path.exists():
        with open(manifest_path, 'r') as f:
            existing = json.load(f)

    manifest = []
    for i, f in enumerate(files, 1):
        title = f['title']
        url = f['url']
        checksum = f.get('checksum', '')
        checksum_type = f.get('checksum_type', 'sha256')

        # Check if already in manifest and valid
        if title in existing:
            existing_entry = existing[title]
            if existing_entry.get('checksum') == checksum:
                if verbose:
                    print(f"[{i}/{len(files)}] {title} (cached)")
                manifest.append(existing_entry)
                continue

        dest = dest_dir / title
        try:
            if verbose:
                print(f"[{i}/{len(files)}] {title}...")
            fetch(url, str(dest), checksum=checksum, checksum_type=checksum_type)
            manifest_entry = {
                'title': title,
                'url': url,
                'size': f.get('size', 0),
                'checksum': checksum,
                'checksum_type': checksum_type,
            }
            manifest.append(manifest_entry)
            if verbose:
                print(f"  ✓ {format_size(f.get('size', 0))}")
        except RuntimeError as e:
            print(f"  ✗ Error: {e}", file=sys.stderr)
            # Don't add to manifest if download failed

    # Write manifest
    with open(manifest_path, 'w') as f:
        json.dump(manifest, f, indent=2)
    if verbose:
        print(f"\nManifest written to {manifest_path}")


def main():
    parser = argparse.ArgumentParser(
        description='Download CMIP6 DAMIP data from ESGF',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s --model IPSL-CM6A-LR --experiment hist-aer --dry-run
  %(prog)s --model IPSL-CM6A-LR --experiment hist-aer
  %(prog)s --all-m5 --dry-run
  %(prog)s --all-m4 --experiment hist-aer
        """)

    parser.add_argument('--model', help='CMIP6 source_id (e.g. IPSL-CM6A-LR)')
    parser.add_argument('--experiment', default='hist-aer',
                        help='DAMIP experiment_id (default: hist-aer)')
    parser.add_argument('--variant', help='variant_label filter (optional)')
    parser.add_argument('--grid', help='grid_label filter (optional)')
    parser.add_argument('--all-m4', action='store_true',
                        help='Download all M4 models (3 models)')
    parser.add_argument('--all-m5', action='store_true',
                        help='Download all M5 models (8 models)')
    parser.add_argument('--dry-run', action='store_true',
                        help='Report files and volume without downloading')
    parser.add_argument('--output-dir', default='raw_data/cmip6_damip',
                        help='Base directory for downloads (default: raw_data/cmip6_damip)')

    args = parser.parse_args()

    # Determine models to download
    if args.all_m5:
        models = M5_MODELS
    elif args.all_m4:
        models = M4_MODELS
    elif args.model:
        models = [args.model]
    else:
        parser.print_help()
        sys.exit(1)

    # Base and warm year ranges (per plan §1.2, Appendix D)
    base_years = (1850, 1859)
    warm_years_default = (2011, 2020)
    warm_years_cesm2 = (2005, 2014)  # CESM2 data ends at 2014-12

    experiment = args.experiment

    # All-model dry-run: report combined size
    if (args.all_m4 or args.all_m5) and args.dry_run:
        total_files = 0
        total_gb = 0.0
        for model in models:
            warm_years = warm_years_cesm2 if model == 'CESM2' else warm_years_default
            try:
                files = search_damip_files(
                    model=model,
                    experiment=experiment,
                    base_years=base_years,
                    warm_years=warm_years,
                    variant=args.variant,
                    grid=args.grid,
                )
                nf, gb = dry_run_report(model, experiment, files)
                total_files += nf
                total_gb += gb
            except RuntimeError as e:
                print(f"Error for {model}: {e}", file=sys.stderr)

        print(f"\n{'='*70}")
        print(f"COMBINED TOTAL ({len(models)} models): {total_files} files, {total_gb:.2f} GB")
        print(f"{'='*70}\n")
        return

    # Single model processing
    for model in models:
        warm_years = warm_years_cesm2 if model == 'CESM2' else warm_years_default

        try:
            files = search_damip_files(
                model=model,
                experiment=experiment,
                base_years=base_years,
                warm_years=warm_years,
                variant=args.variant,
                grid=args.grid,
            )
        except RuntimeError as e:
            print(f"Error for {model}: {e}", file=sys.stderr)
            continue

        if args.dry_run:
            dry_run_report(model, experiment, files)
        else:
            print(f"Downloading {model} / {experiment}...")
            download_files(
                model=model,
                experiment=experiment,
                files=files,
                output_dir=args.output_dir,
                verbose=True,
            )


if __name__ == '__main__':
    main()
