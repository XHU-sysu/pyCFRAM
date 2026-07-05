"""Pure stdlib ESGF CMIP6 search and download client.

Uses only urllib (no requests, no third-party HTTP libs) for CMIP6 DAMIP data discovery
and download from ESGF nodes. Supports free-tier download without authentication.

Key features:
- Solr-based dataset/file search (DKRZ, CEDA nodes)
- HTTP Range request resume (206 Partial Content)
- SHA256/MD5 checksum verification
- manifest.json tracking for idempotent re-runs
- Node fallback on network errors
"""

import json
import os
import hashlib
import urllib.request
import urllib.parse
import urllib.error
import re
from pathlib import Path
from typing import Optional, List, Dict, Tuple


def search_datasets(
    model: str,
    experiment: str,
    variables: Optional[List[str]] = None,
    variant: Optional[str] = None,
    grid: Optional[str] = None,
    node: str = "https://esgf-data.dkrz.de/esg-search/search"
) -> List[Dict]:
    """Search ESGF Solr for CMIP6 DAMIP datasets.

    Args:
        model: CMIP6 source_id (e.g. "IPSL-CM6A-LR")
        experiment: experiment_id (e.g. "hist-aer")
        variables: list of variable_ids to include in search (optional)
        variant: variant_label filter (e.g. "r1i1p1f1", optional)
        grid: grid_label filter (e.g. "gn", optional)
        node: ESGF search node URL

    Returns:
        List of dataset dicts with keys: id, variant_label, grid_label, version
    """
    params = {
        'project': 'CMIP6',
        'activity_id': 'DAMIP',
        'experiment_id': experiment,
        'source_id': model,
        'table_id': 'Amon',
        'type': 'Dataset',
        'format': 'application/solr+json',
        'limit': 500,
        'latest': 'true',
    }
    # Note: deliberately NOT filtering replica=false here. Some models (e.g.
    # MRI-ESM2-0) are indexed as replica=true at every node this client
    # queries (DKRZ/CEDA aren't the model's home data node), so replica=false
    # returns zero results even though the data is genuinely downloadable.
    # Cross-data-node duplicates are instead deduplicated client-side by
    # filename in download_damip.search_damip_files().
    if variant:
        params['variant_label'] = variant
    if grid:
        params['grid_label'] = grid
    if variables:
        # Search is AND-ed over variables, so we search separately per variable
        # and union results by dataset_id
        all_datasets = {}
        for var in variables:
            params_var = params.copy()
            params_var['variable_id'] = var
            url = node + '?' + urllib.parse.urlencode(params_var)
            try:
                with urllib.request.urlopen(url, timeout=30) as response:
                    data = json.loads(response.read())
                    for doc in data.get('response', {}).get('docs', []):
                        did = doc['id']
                        if did not in all_datasets:
                            all_datasets[did] = doc
            except urllib.error.URLError as e:
                raise RuntimeError(f"ESGF search failed ({node}): {e}")
        return list(all_datasets.values())
    else:
        # No variable filter: just search once
        url = node + '?' + urllib.parse.urlencode(params)
        try:
            with urllib.request.urlopen(url, timeout=30) as response:
                data = json.loads(response.read())
                return data.get('response', {}).get('docs', [])
        except urllib.error.URLError as e:
            raise RuntimeError(f"ESGF search failed ({node}): {e}")


def list_files(
    dataset_id: str,
    node: str = "https://esgf-data.dkrz.de/esg-search/search"
) -> List[Dict]:
    """List all files in a CMIP6 DAMIP dataset.

    Args:
        dataset_id: ESGF dataset_id
        node: ESGF search node URL

    Returns:
        List of file dicts with keys: title, size, url, checksum, checksum_type
    """
    params = {
        'dataset_id': dataset_id,
        'type': 'File',
        'format': 'application/solr+json',
        'limit': 10000,
        'fields': 'title,size,url,checksum,checksum_type',
        'latest': 'true',
    }
    url = node + '?' + urllib.parse.urlencode(params)
    try:
        with urllib.request.urlopen(url, timeout=30) as response:
            data = json.loads(response.read())
            files = []
            for doc in data.get('response', {}).get('docs', []):
                # url field is pipe-separated: "http://...|HTTPServer|...|..."
                # Extract HTTPServer URL
                url_str = doc.get('url', '')
                if isinstance(url_str, list):
                    url_str = url_str[0] if url_str else ''

                # Parse pipe-separated format
                parts = url_str.split('|')
                http_url = None
                for i, part in enumerate(parts):
                    if i > 0 and parts[i - 1] == 'HTTPServer' and part.startswith('http'):
                        http_url = part
                        break
                    elif part.startswith('http') and i + 1 < len(parts) and parts[i + 1] == 'HTTPServer':
                        http_url = part
                        break

                # Fallback: just take first URL that looks like http
                if not http_url:
                    for part in parts:
                        if part.startswith('http'):
                            http_url = part
                            break

                if http_url:
                    # checksum/checksum_type are Solr multi-valued fields
                    # (returned as lists, e.g. checksum_type=['SHA256']) --
                    # unwrap to a scalar. Files can carry more than one
                    # checksum type; take the first.
                    checksum = doc.get('checksum', '')
                    if isinstance(checksum, list):
                        checksum = checksum[0] if checksum else ''
                    checksum_type = doc.get('checksum_type', 'sha256')
                    if isinstance(checksum_type, list):
                        checksum_type = checksum_type[0] if checksum_type else 'sha256'

                    files.append({
                        'title': doc.get('title', ''),
                        'size': doc.get('size', 0),
                        'url': http_url,
                        'checksum': checksum,
                        'checksum_type': checksum_type,
                    })
            return files
    except urllib.error.URLError as e:
        raise RuntimeError(f"ESGF file list failed ({node}): {e}")


def filename_time_overlap(filename: str, base_years: Tuple[int, int],
                          warm_years: Tuple[int, int]) -> bool:
    """Check if file's time range overlaps with base or warm period.

    CMIP6 filenames typically end with _YYYYMM-YYYYMM.nc
    Extract start/end years and check overlap.

    Args:
        filename: CMIP6 filename
        base_years: (start_year, end_year) tuple
        warm_years: (start_year, end_year) tuple

    Returns:
        True if file overlaps either base or warm period
    """
    # Try to extract _YYYYMM-YYYYMM pattern
    match = re.search(r'_(\d{6})-(\d{6})\.nc$', filename)
    if not match:
        # Can't parse time range, include file (better conservative)
        return True

    file_start_yyyymm = match.group(1)
    file_end_yyyymm = match.group(2)

    try:
        file_start_year = int(file_start_yyyymm[:4])
        file_end_year = int(file_end_yyyymm[:4])
    except ValueError:
        return True

    base_start, base_end = base_years
    warm_start, warm_end = warm_years

    # Check overlap: [file_start, file_end] ∩ [base_start, base_end] or [warm_start, warm_end]
    base_overlap = not (file_end_year < base_start or file_start_year > base_end)
    warm_overlap = not (file_end_year < warm_start or file_start_year > warm_end)

    return base_overlap or warm_overlap


def fetch(
    url: str,
    dest: str,
    checksum: Optional[str] = None,
    checksum_type: str = 'sha256',
    resume: bool = True,
    manifest_path: Optional[str] = None,
) -> bool:
    """Download a single file from ESGF with resume and checksum support.

    Args:
        url: HTTP URL of file
        dest: destination path (created with parent dirs if needed)
        checksum: expected checksum value (optional)
        checksum_type: 'sha256' or 'md5'
        resume: support HTTP 206 range request if file partially exists
        manifest_path: optional path to JSON manifest (unused in this function,
                       passed for API compatibility)

    Returns:
        True if successful (and checksum OK if provided)
    """
    dest_path = Path(dest)
    dest_path.parent.mkdir(parents=True, exist_ok=True)

    # Check if file already exists and is complete
    if dest_path.exists():
        if checksum:
            actual = compute_checksum(dest, checksum_type)
            if actual == checksum:
                return True
            else:
                # Checksum mismatch, redownload
                dest_path.unlink()
        else:
            # No checksum provided; assume existing file is OK
            return True

    # Download with resume support
    req = urllib.request.Request(url)

    if resume and dest_path.exists():
        # Request range starting from current file size
        current_size = dest_path.stat().st_size
        req.add_header('Range', f'bytes={current_size}-')

    try:
        with urllib.request.urlopen(req, timeout=60) as response:
            mode = 'ab' if (resume and dest_path.exists()) else 'wb'
            with open(dest, mode) as f:
                while True:
                    chunk = response.read(8192)
                    if not chunk:
                        break
                    f.write(chunk)
    except urllib.error.HTTPError as e:
        if e.code == 416:
            # Range not satisfiable (file complete on server)
            pass
        else:
            raise RuntimeError(f"Download failed ({url}): {e}")
    except urllib.error.URLError as e:
        raise RuntimeError(f"Download failed ({url}): {e}")

    # Verify checksum if provided
    if checksum:
        actual = compute_checksum(dest, checksum_type)
        if actual != checksum:
            raise RuntimeError(
                f"Checksum mismatch for {dest}: "
                f"expected {checksum}, got {actual}"
            )

    return True


def compute_checksum(filepath: str, checksum_type: str = 'sha256') -> str:
    """Compute SHA256 or MD5 checksum of a file."""
    if checksum_type.lower() == 'sha256':
        hasher = hashlib.sha256()
    elif checksum_type.lower() == 'md5':
        hasher = hashlib.md5()
    else:
        raise ValueError(f"Unknown checksum type: {checksum_type}")

    with open(filepath, 'rb') as f:
        while True:
            chunk = f.read(65536)
            if not chunk:
                break
            hasher.update(chunk)

    return hasher.hexdigest()


def write_manifest(manifest_path: str, files: List[Dict]):
    """Write download manifest (list of {url, size, checksum, mtime} per file)."""
    os.makedirs(os.path.dirname(manifest_path) or '.', exist_ok=True)
    with open(manifest_path, 'w') as f:
        json.dump(files, f, indent=2)


def read_manifest(manifest_path: str) -> Dict[str, Dict]:
    """Read download manifest and return dict keyed by filename."""
    if not os.path.exists(manifest_path):
        return {}
    with open(manifest_path, 'r') as f:
        files = json.load(f)
    return {f['title']: f for f in files if isinstance(f, dict)}


def search_with_fallback(
    model: str,
    experiment: str,
    variables: Optional[List[str]] = None,
    variant: Optional[str] = None,
    grid: Optional[str] = None,
    nodes: Optional[List[str]] = None,
) -> Tuple[List[Dict], str]:
    """Search for datasets, trying multiple nodes in order.

    Returns:
        (datasets, successful_node)
    """
    if nodes is None:
        nodes = [
            "https://esgf-data.dkrz.de/esg-search/search",
            "https://esgf.ceda.ac.uk/esg-search/search",
        ]

    last_error = None
    for node in nodes:
        try:
            datasets = search_datasets(model, experiment, variables, variant, grid, node)
            if datasets:
                return datasets, node
        except RuntimeError as e:
            last_error = e
            continue

    if last_error:
        raise last_error
    raise RuntimeError("No datasets found on any ESGF node")
