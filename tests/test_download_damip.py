"""Unit tests for scripts/download_damip.py's multi-mirror fallback in
download_files() (Phase 3, real bug/fix from the M5 download campaign).

CNRM-CM6-1's home ESGF data node (esg1.umr-cnrm.fr) timed out on 7 files
(cl/clw/cli x2 chunks + rsdt) while every other variable downloaded fine.
A live check confirmed a genuine replica existed at a different data node
(esgf3.dkrz.de, a fully independent URL/thredds path -- not just a second
index entry pointing at the same host). search_damip_files() used to keep
only the first-seen URL per filename and drop the rest; it now collects
every known mirror URL per file, and download_files() tries them in order,
falling back to the next mirror on a RuntimeError instead of giving up.

All tests here monkeypatch scripts.download_damip.fetch -- no real network
calls.
"""
import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))
from scripts import download_damip as dd


def test_download_files_falls_back_to_second_mirror_on_first_failure(tmp_path, monkeypatch):
    calls = []

    def fake_fetch(url, dest, checksum=None, checksum_type='sha256', expected_size=0):
        calls.append(url)
        if 'mirror-a' in url:
            raise RuntimeError('Download failed (timed out)')
        # mirror-b succeeds: actually create the file, like the real fetch() would.
        Path(dest).write_bytes(b'fake file content')
        return True

    monkeypatch.setattr(dd, 'fetch', fake_fetch)

    files = [{
        'title': 'cl_Amon_test.nc',
        'size': 100,
        'urls': ['http://mirror-a.example/cl_Amon_test.nc', 'http://mirror-b.example/cl_Amon_test.nc'],
        'checksum': 'abc123',
        'checksum_type': 'sha256',
    }]

    dd.download_files('TESTMODEL', 'hist-aer', files, output_dir=str(tmp_path))

    assert calls == [
        'http://mirror-a.example/cl_Amon_test.nc',
        'http://mirror-b.example/cl_Amon_test.nc',
    ]
    assert (tmp_path / 'TESTMODEL' / 'hist-aer' / 'cl_Amon_test.nc').read_bytes() == b'fake file content'

    manifest_path = tmp_path / 'TESTMODEL' / 'hist-aer' / 'manifest.json'
    manifest = json.loads(manifest_path.read_text())
    # Manifest records the mirror that actually succeeded, not the first-tried one.
    assert manifest[0]['url'] == 'http://mirror-b.example/cl_Amon_test.nc'


def test_download_files_reports_failure_only_after_all_mirrors_exhausted(tmp_path, monkeypatch):
    calls = []

    def fake_fetch(url, dest, checksum=None, checksum_type='sha256', expected_size=0):
        calls.append(url)
        raise RuntimeError('Download failed (timed out)')

    monkeypatch.setattr(dd, 'fetch', fake_fetch)

    files = [{
        'title': 'cl_Amon_test.nc',
        'size': 100,
        'urls': ['http://mirror-a.example/x.nc', 'http://mirror-b.example/x.nc'],
        'checksum': '', 'checksum_type': 'sha256',
    }]

    dd.download_files('TESTMODEL', 'hist-aer', files, output_dir=str(tmp_path))

    assert len(calls) == 2  # both mirrors attempted
    manifest_path = tmp_path / 'TESTMODEL' / 'hist-aer' / 'manifest.json'
    manifest = json.loads(manifest_path.read_text())
    assert manifest == []  # never recorded -- failed file must not appear as if it succeeded
    assert not (tmp_path / 'TESTMODEL' / 'hist-aer' / 'cl_Amon_test.nc').exists()


def test_download_files_first_mirror_success_does_not_try_second(tmp_path, monkeypatch):
    calls = []

    def fake_fetch(url, dest, checksum=None, checksum_type='sha256', expected_size=0):
        calls.append(url)
        Path(dest).write_bytes(b'ok')
        return True

    monkeypatch.setattr(dd, 'fetch', fake_fetch)

    files = [{
        'title': 'ta_Amon_test.nc',
        'size': 50,
        'urls': ['http://mirror-a.example/ta.nc', 'http://mirror-b.example/ta.nc'],
        'checksum': '', 'checksum_type': 'sha256',
    }]

    dd.download_files('TESTMODEL', 'hist-aer', files, output_dir=str(tmp_path))
    assert calls == ['http://mirror-a.example/ta.nc']  # second mirror never touched


def test_search_damip_files_collects_multiple_urls_for_same_filename(monkeypatch):
    """search_damip_files must merge URLs across dataset entries for the
    SAME filename into one `urls` list, rather than keeping only the
    first-seen dataset's URL (the real bug: CNRM-CM6-1's DKRZ mirror was
    silently dropped)."""

    def fake_search_with_fallback(model, experiment, variables=None, variant=None,
                                   grid=None, nodes=None):
        datasets = [
            {'id': 'ds-home-node'},
            {'id': 'ds-mirror-node'},
        ]
        return datasets, nodes[0] if nodes else 'https://esgf-data.dkrz.de/esg-search/search'

    def fake_list_files(dataset_id, node=None):
        if dataset_id == 'ds-home-node':
            return [{'title': 'cl_Amon_test_185001-202012.nc', 'size': 100,
                      'url': 'http://home.example/cl.nc', 'checksum': 'x', 'checksum_type': 'sha256'}]
        return [{'title': 'cl_Amon_test_185001-202012.nc', 'size': 100,
                  'url': 'http://mirror.example/cl.nc', 'checksum': 'x', 'checksum_type': 'sha256'}]

    monkeypatch.setattr(dd, 'search_with_fallback', fake_search_with_fallback)
    monkeypatch.setattr(dd, 'list_files', fake_list_files)
    monkeypatch.setattr(dd, 'get_variables_for_model', lambda model, experiment, nodes: ['cl'])

    files = dd.search_damip_files('TESTMODEL', 'hist-aer', (1850, 1859), (2011, 2020))
    assert len(files) == 1
    assert set(files[0]['urls']) == {'http://home.example/cl.nc', 'http://mirror.example/cl.nc'}
    assert 'url' not in files[0]  # replaced by the plural 'urls' list
