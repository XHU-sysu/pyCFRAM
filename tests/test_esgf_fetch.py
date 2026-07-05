"""Unit tests for data/esgf_fetch.py (Phase 3, WP-M4.6).

Pure-stdlib ESGF Solr search + HTTP download client. All network I/O is
monkeypatched (urllib.request.urlopen) to canned Solr JSON / fake HTTP
responses per docs/plan_ph3.md §10 -- no real network calls are made.

Coverage priorities (per WP-M4.6 brief):
  1. list_files' checksum/checksum_type list-vs-scalar unwrapping (a real
     production bug -- commit 2d13980 "fix(Ph3): unwrap list-typed
     checksum/checksum_type Solr fields" -- already fixed in this branch;
     this file adds the regression test that was missing).
  2. list_files' HTTPServer URL extraction out of the pipe-separated `url`
     Solr field, in both token orders plus the no-label fallback.
  3. search_with_fallback's node-fallback control flow (first-node success,
     empty-then-fallback, error-then-fallback, all-fail, all-empty).
  4. filename_time_overlap's regex parsing edge cases.

Also covers compute_checksum, write_manifest/read_manifest, and fetch()'s
control flow (existing-file short-circuits, checksum mismatch/redownload,
HTTP error handling).

Two genuine (pre-existing, NOT introduced or fixed by this test pass) bugs
were found while writing these tests. Per the WP-M4.6 brief they are
documented/reproduced here, not fixed (production code is out of scope for
this pass) -- see the two tests marked "KNOWN ISSUE" below, and the
corresponding docstrings, for details:
  - list_files(): when Solr's `url` field is multi-valued (a list, one
    pipe-string per access service -- the documented ESGF schema), only
    url_str[0] is ever inspected. If HTTPServer isn't the first service in
    the list (Solr gives no ordering guarantee), the file is silently
    dropped from the result instead of being found further down the list.
  - fetch(): the HTTP-Range "resume" path is unreachable dead code. By the
    time the download block runs, dest_path.exists() is always False
    (either the file never existed, or the checksum-mismatch branch just
    unlinked it) -- so `resume` never adds a Range header and mode is
    always 'wb', never 'ab'. A knock-on effect: a pre-existing file with NO
    checksum recorded is unconditionally trusted as "complete" (even a
    truncated partial download from an interrupted run), and an HTTP 416
    response on a fresh (no pre-existing file) attempt makes fetch()
    return True without ever creating the destination file.
"""
import hashlib
import json
import os
import sys
import urllib.error
import urllib.parse
import urllib.request

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from data import esgf_fetch as esgf


# ---------------------------------------------------------------------------
# Test helpers: fake urllib.request.urlopen responses (no network)
# ---------------------------------------------------------------------------

class _FakeHTTPResponse:
    """Minimal stand-in for the context-manager object urlopen() returns."""

    def __init__(self, payload: bytes):
        self._payload = payload
        self._pos = 0

    def read(self, size=-1):
        # Must behave like a real chunked-read stream (fetch()'s download
        # loop calls response.read(8192) until it gets back b'') -- a
        # helper that always returns the full payload regardless of `size`
        # never signals EOF and spins fetch() into an infinite loop.
        if size is None or size < 0:
            chunk = self._payload[self._pos:]
        else:
            chunk = self._payload[self._pos:self._pos + size]
        self._pos += len(chunk)
        return chunk

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False


def _solr_payload(docs):
    return json.dumps({'response': {'docs': docs}}).encode('utf-8')


# ---------------------------------------------------------------------------
# search_datasets
# ---------------------------------------------------------------------------

def test_search_datasets_no_variable_filter_single_request(monkeypatch):
    docs = [{'id': 'CMIP6.DAMIP.IPSL.hist-aer.v1'}]
    captured = {}

    def fake_urlopen(url, timeout=30):
        captured['url'] = url
        return _FakeHTTPResponse(_solr_payload(docs))

    monkeypatch.setattr(urllib.request, 'urlopen', fake_urlopen)
    result = esgf.search_datasets('IPSL-CM6A-LR', 'hist-aer',
                                   node='https://node.example/esg-search/search')
    assert result == docs
    qs = urllib.parse.parse_qs(urllib.parse.urlsplit(captured['url']).query)
    assert 'variable_id' not in qs
    assert qs['source_id'] == ['IPSL-CM6A-LR']
    assert qs['experiment_id'] == ['hist-aer']


def test_search_datasets_variant_and_grid_params_included(monkeypatch):
    captured = {}

    def fake_urlopen(url, timeout=30):
        captured['url'] = url
        return _FakeHTTPResponse(_solr_payload([]))

    monkeypatch.setattr(urllib.request, 'urlopen', fake_urlopen)
    esgf.search_datasets('M', 'E', variant='r1i1p1f1', grid='gn',
                          node='https://node.example/esg-search/search')
    qs = urllib.parse.parse_qs(urllib.parse.urlsplit(captured['url']).query)
    assert qs['variant_label'] == ['r1i1p1f1']
    assert qs['grid_label'] == ['gn']


def test_search_datasets_variables_union_dedups_by_dataset_id(monkeypatch):
    """Search is AND-ed over variables server-side, so search_datasets()
    queries once per variable and unions results by dataset id. ds2 is
    returned by both the 'ta' and 'hus' per-variable searches here and
    must appear exactly once in the final result."""
    calls = []

    def fake_urlopen(url, timeout=30):
        qs = urllib.parse.parse_qs(urllib.parse.urlsplit(url).query)
        var = qs['variable_id'][0]
        calls.append(var)
        if var == 'ta':
            docs = [{'id': 'ds1'}, {'id': 'ds2'}]
        else:
            docs = [{'id': 'ds2'}, {'id': 'ds3'}]
        return _FakeHTTPResponse(_solr_payload(docs))

    monkeypatch.setattr(urllib.request, 'urlopen', fake_urlopen)
    result = esgf.search_datasets('M', 'E', variables=['ta', 'hus'])
    assert calls == ['ta', 'hus']
    assert sorted(d['id'] for d in result) == ['ds1', 'ds2', 'ds3']


def test_search_datasets_url_error_raises_runtimeerror_no_variables(monkeypatch):
    def fake_urlopen(url, timeout=30):
        raise urllib.error.URLError('network down')

    monkeypatch.setattr(urllib.request, 'urlopen', fake_urlopen)
    with pytest.raises(RuntimeError, match='ESGF search failed'):
        esgf.search_datasets('M', 'E')


def test_search_datasets_url_error_raises_runtimeerror_with_variables(monkeypatch):
    def fake_urlopen(url, timeout=30):
        raise urllib.error.URLError('network down')

    monkeypatch.setattr(urllib.request, 'urlopen', fake_urlopen)
    with pytest.raises(RuntimeError, match='ESGF search failed'):
        esgf.search_datasets('M', 'E', variables=['ta'])


# ---------------------------------------------------------------------------
# list_files -- HTTPServer URL extraction out of the pipe-separated field
# ---------------------------------------------------------------------------

def test_list_files_httpserver_label_before_url(monkeypatch):
    """parts = [..., 'HTTPServer', 'http://...'] -- label precedes the URL."""
    docs = [{
        'title': 'ta_Amon_test.nc', 'size': 123,
        'url': 'esgf-data.example|HTTPServer|http://real.example/ta_Amon_test.nc',
        'checksum': 'abc123', 'checksum_type': 'SHA256',
    }]

    def fake_urlopen(url, timeout=30):
        return _FakeHTTPResponse(_solr_payload(docs))

    monkeypatch.setattr(urllib.request, 'urlopen', fake_urlopen)
    files = esgf.list_files('ds1')
    assert len(files) == 1
    assert files[0]['url'] == 'http://real.example/ta_Amon_test.nc'
    assert files[0]['checksum'] == 'abc123'
    assert files[0]['checksum_type'] == 'SHA256'


def test_list_files_httpserver_label_after_url(monkeypatch):
    """parts = ['http://...', 'HTTPServer', ...] -- URL precedes the label
    (the common real-world ESGF Solr ordering)."""
    docs = [{
        'title': 'ta_Amon_test.nc', 'size': 1,
        'url': 'http://real.example/ta_Amon_test.nc|HTTPServer|download',
    }]

    def fake_urlopen(url, timeout=30):
        return _FakeHTTPResponse(_solr_payload(docs))

    monkeypatch.setattr(urllib.request, 'urlopen', fake_urlopen)
    files = esgf.list_files('ds1')
    assert files[0]['url'] == 'http://real.example/ta_Amon_test.nc'
    # doc carries no checksum fields at all -> defaults apply
    assert files[0]['checksum'] == ''
    assert files[0]['checksum_type'] == 'sha256'


def test_list_files_fallback_when_no_httpserver_label_present(monkeypatch):
    """Neither token-order branch matches (no literal 'HTTPServer' in the
    pipe-separated string at all) -> fallback picks the first http-looking
    part."""
    docs = [{
        'title': 'x.nc', 'size': 1,
        'url': ('ftp://mirror.example/x.nc|GridFTP|download'
                 '|http://fallback.example/x.nc|OtherService|foo'),
    }]

    def fake_urlopen(url, timeout=30):
        return _FakeHTTPResponse(_solr_payload(docs))

    monkeypatch.setattr(urllib.request, 'urlopen', fake_urlopen)
    files = esgf.list_files('ds1')
    assert files[0]['url'] == 'http://fallback.example/x.nc'


def test_list_files_no_http_url_at_all_is_dropped(monkeypatch):
    docs = [{'title': 'y.nc', 'size': 1, 'url': 'ftp://only.example/y.nc|FTP|download'}]

    def fake_urlopen(url, timeout=30):
        return _FakeHTTPResponse(_solr_payload(docs))

    monkeypatch.setattr(urllib.request, 'urlopen', fake_urlopen)
    files = esgf.list_files('ds1')
    assert files == []


def test_list_files_url_error_raises_runtimeerror(monkeypatch):
    def fake_urlopen(url, timeout=30):
        raise urllib.error.URLError('down')

    monkeypatch.setattr(urllib.request, 'urlopen', fake_urlopen)
    with pytest.raises(RuntimeError, match='ESGF file list failed'):
        esgf.list_files('ds1')


# ---------------------------------------------------------------------------
# list_files -- checksum/checksum_type list-vs-scalar unwrapping
# (regression test for commit 2d13980's production bug fix)
# ---------------------------------------------------------------------------

def test_list_files_checksum_and_checksum_type_scalar_passthrough(monkeypatch):
    docs = [{
        'title': 'z.nc', 'size': 1,
        'url': 'http://real.example/z.nc|HTTPServer|download',
        'checksum': 'plainscalar', 'checksum_type': 'md5',
    }]

    def fake_urlopen(url, timeout=30):
        return _FakeHTTPResponse(_solr_payload(docs))

    monkeypatch.setattr(urllib.request, 'urlopen', fake_urlopen)
    files = esgf.list_files('ds1')
    assert files[0]['checksum'] == 'plainscalar'
    assert files[0]['checksum_type'] == 'md5'


def test_list_files_checksum_and_checksum_type_list_unwrapped_to_scalar(monkeypatch):
    """Solr's File-type records return checksum/checksum_type as
    multi-valued fields (e.g. checksum_type=['SHA256']). Regression guard
    for commit 2d13980: compute_checksum() crashed on .lower() against a
    list the first time a real download ran this code path."""
    docs = [{
        'title': 'z.nc', 'size': 1,
        'url': 'http://real.example/z.nc|HTTPServer|download',
        'checksum': ['deadbeef', 'ignored-second-value'],
        'checksum_type': ['SHA256'],
    }]

    def fake_urlopen(url, timeout=30):
        return _FakeHTTPResponse(_solr_payload(docs))

    monkeypatch.setattr(urllib.request, 'urlopen', fake_urlopen)
    files = esgf.list_files('ds1')
    assert files[0]['checksum'] == 'deadbeef'
    assert files[0]['checksum_type'] == 'SHA256'
    # And the unwrapped scalar must actually work with compute_checksum's
    # .lower() call (the exact call that crashed in production pre-fix).
    assert files[0]['checksum_type'].lower() in ('sha256', 'md5')


def test_list_files_checksum_and_checksum_type_empty_list_unwrapped_to_defaults(monkeypatch):
    docs = [{
        'title': 'z.nc', 'size': 1,
        'url': 'http://real.example/z.nc|HTTPServer|download',
        'checksum': [], 'checksum_type': [],
    }]

    def fake_urlopen(url, timeout=30):
        return _FakeHTTPResponse(_solr_payload(docs))

    monkeypatch.setattr(urllib.request, 'urlopen', fake_urlopen)
    files = esgf.list_files('ds1')
    assert files[0]['checksum'] == ''
    assert files[0]['checksum_type'] == 'sha256'


def test_list_files_url_as_list_takes_first_element_only_KNOWN_ISSUE(monkeypatch):
    """KNOWN ISSUE (pre-existing, not fixed here -- reported separately).

    ESGF Solr's `url` field can itself be multi-valued: a list with one
    pipe-string per access service (HTTPServer/OPENDAP/Globus/GridFTP),
    order not guaranteed. list_files() only ever inspects url_str[0]
    (line ~127: `url_str = url_str[0] if url_str else ''`) and never scans
    the remaining list entries. If the HTTPServer entry isn't first, the
    file is silently DROPPED from the result -- even though a perfectly
    good HTTPServer download URL exists later in the same list.

    This test pins the CURRENT (buggy) behavior so a future real fix shows
    up as a deliberate, visible diff rather than a silent behavior change.
    """
    docs = [{
        'title': 'dropped.nc', 'size': 1,
        'url': [
            'globus://x.example/dropped.nc|Globus|Globus',
            'http://real.example/dropped.nc|HTTPServer|download',
        ],
    }]

    def fake_urlopen(url, timeout=30):
        return _FakeHTTPResponse(_solr_payload(docs))

    monkeypatch.setattr(urllib.request, 'urlopen', fake_urlopen)
    files = esgf.list_files('ds1')
    # Current (buggy) behavior: dropped, despite a valid HTTPServer URL
    # existing at url[1].
    assert files == []


# ---------------------------------------------------------------------------
# filename_time_overlap
# ---------------------------------------------------------------------------

def test_filename_time_overlap_no_time_range_in_filename_is_conservatively_true():
    assert esgf.filename_time_overlap('ta_Amon_model_fixed.nc',
                                       (1850, 1851), (1854, 1855)) is True


def test_filename_time_overlap_single_span_file_covers_both_periods():
    """One file spans 1850-01 through 1855-12 -- covers both the base
    (1850-1851) and warm (1854-1855) windows in a single file."""
    assert esgf.filename_time_overlap('ta_Amon_model_185001-185512.nc',
                                       (1850, 1851), (1854, 1855)) is True


def test_filename_time_overlap_non_overlapping_range_is_false():
    assert esgf.filename_time_overlap('ta_Amon_model_190001-190012.nc',
                                       (1850, 1851), (1854, 1855)) is False


def test_filename_time_overlap_overlaps_base_period_only():
    assert esgf.filename_time_overlap('ta_Amon_model_184901-185006.nc',
                                       (1850, 1851), (1900, 1901)) is True


def test_filename_time_overlap_overlaps_warm_period_only():
    assert esgf.filename_time_overlap('ta_Amon_model_185501-185512.nc',
                                       (1800, 1801), (1854, 1856)) is True


def test_filename_time_overlap_exact_boundary_year_counts_as_overlap():
    # File ends exactly at base_start's year -> inclusive boundary overlap.
    assert esgf.filename_time_overlap('ta_Amon_model_184901-185001.nc',
                                       (1850, 1851), (1900, 1901)) is True


def test_filename_time_overlap_adjacent_year_before_start_is_no_overlap():
    # File ends the year BEFORE base starts -> no overlap.
    assert esgf.filename_time_overlap('ta_Amon_model_184801-184912.nc',
                                       (1850, 1851), (1900, 1901)) is False


# ---------------------------------------------------------------------------
# search_with_fallback -- node fallback control flow
# ---------------------------------------------------------------------------

def test_search_with_fallback_first_node_succeeds_second_not_tried(monkeypatch):
    calls = []

    def fake_search(model, experiment, variables, variant, grid, node):
        calls.append(node)
        if node == 'nodeA':
            return [{'id': 'ds1'}]
        raise AssertionError('nodeB should never be queried once nodeA succeeds')

    monkeypatch.setattr(esgf, 'search_datasets', fake_search)
    datasets, node = esgf.search_with_fallback('M', 'E', nodes=['nodeA', 'nodeB'])
    assert node == 'nodeA'
    assert datasets == [{'id': 'ds1'}]
    assert calls == ['nodeA']


def test_search_with_fallback_empty_result_falls_through_to_second_node(monkeypatch):
    calls = []

    def fake_search(model, experiment, variables, variant, grid, node):
        calls.append(node)
        if node == 'nodeA':
            return []
        return [{'id': 'ds2'}]

    monkeypatch.setattr(esgf, 'search_datasets', fake_search)
    datasets, node = esgf.search_with_fallback('M', 'E', nodes=['nodeA', 'nodeB'])
    assert node == 'nodeB'
    assert datasets == [{'id': 'ds2'}]
    assert calls == ['nodeA', 'nodeB']


def test_search_with_fallback_error_on_first_node_falls_through(monkeypatch):
    def fake_search(model, experiment, variables, variant, grid, node):
        if node == 'nodeA':
            raise RuntimeError('nodeA unreachable')
        return [{'id': 'ds3'}]

    monkeypatch.setattr(esgf, 'search_datasets', fake_search)
    datasets, node = esgf.search_with_fallback('M', 'E', nodes=['nodeA', 'nodeB'])
    assert node == 'nodeB'
    assert datasets == [{'id': 'ds3'}]


def test_search_with_fallback_all_nodes_error_reraises_last_error(monkeypatch):
    def fake_search(model, experiment, variables, variant, grid, node):
        raise RuntimeError('fail:%s' % node)

    monkeypatch.setattr(esgf, 'search_datasets', fake_search)
    with pytest.raises(RuntimeError, match='fail:nodeB'):
        esgf.search_with_fallback('M', 'E', nodes=['nodeA', 'nodeB'])


def test_search_with_fallback_all_nodes_empty_raises_no_datasets_found(monkeypatch):
    def fake_search(model, experiment, variables, variant, grid, node):
        return []

    monkeypatch.setattr(esgf, 'search_datasets', fake_search)
    with pytest.raises(RuntimeError, match='No datasets found'):
        esgf.search_with_fallback('M', 'E', nodes=['nodeA', 'nodeB'])


def test_search_with_fallback_default_nodes_tries_dkrz_first(monkeypatch):
    calls = []

    def fake_search(model, experiment, variables, variant, grid, node):
        calls.append(node)
        return [{'id': 'ds'}]

    monkeypatch.setattr(esgf, 'search_datasets', fake_search)
    datasets, node = esgf.search_with_fallback('M', 'E')
    assert node == 'https://esgf-data.dkrz.de/esg-search/search'
    assert calls == ['https://esgf-data.dkrz.de/esg-search/search']


# ---------------------------------------------------------------------------
# compute_checksum
# ---------------------------------------------------------------------------

def test_compute_checksum_sha256(tmp_path):
    p = tmp_path / 'f.bin'
    p.write_bytes(b'abc')
    assert esgf.compute_checksum(str(p), 'sha256') == hashlib.sha256(b'abc').hexdigest()


def test_compute_checksum_md5_case_insensitive(tmp_path):
    p = tmp_path / 'f.bin'
    p.write_bytes(b'abc')
    assert esgf.compute_checksum(str(p), 'MD5') == hashlib.md5(b'abc').hexdigest()


def test_compute_checksum_unknown_type_raises_valueerror(tmp_path):
    p = tmp_path / 'f.bin'
    p.write_bytes(b'abc')
    with pytest.raises(ValueError, match='Unknown checksum type'):
        esgf.compute_checksum(str(p), 'crc32')


# ---------------------------------------------------------------------------
# write_manifest / read_manifest
# ---------------------------------------------------------------------------

def test_write_and_read_manifest_roundtrip(tmp_path):
    manifest_path = str(tmp_path / 'sub' / 'manifest.json')
    files = [{'title': 'a.nc', 'size': 1, 'checksum': 'x', 'checksum_type': 'sha256'}]
    esgf.write_manifest(manifest_path, files)
    assert os.path.exists(manifest_path)
    loaded = esgf.read_manifest(manifest_path)
    assert loaded == {'a.nc': files[0]}


def test_read_manifest_missing_file_returns_empty_dict(tmp_path):
    assert esgf.read_manifest(str(tmp_path / 'does_not_exist.json')) == {}


def test_read_manifest_skips_non_dict_entries(tmp_path):
    manifest_path = str(tmp_path / 'manifest.json')
    with open(manifest_path, 'w') as f:
        json.dump([{'title': 'a.nc'}, 'not-a-dict', 42], f)
    loaded = esgf.read_manifest(manifest_path)
    assert loaded == {'a.nc': {'title': 'a.nc'}}


# ---------------------------------------------------------------------------
# fetch()
# ---------------------------------------------------------------------------

def test_fetch_fresh_download_success_with_matching_checksum(tmp_path, monkeypatch):
    content = b'hello world'
    checksum = hashlib.sha256(content).hexdigest()
    dest = tmp_path / 'sub' / 'file.nc'  # parent dir does not exist yet

    def fake_urlopen(req, timeout=60):
        assert req.full_url == 'http://x/file.nc'
        assert not req.has_header('Range')
        return _FakeHTTPResponse(content)

    monkeypatch.setattr(urllib.request, 'urlopen', fake_urlopen)
    ok = esgf.fetch('http://x/file.nc', str(dest), checksum=checksum, checksum_type='sha256')
    assert ok is True
    assert dest.read_bytes() == content


def test_fetch_existing_file_checksum_match_skips_network(tmp_path, monkeypatch):
    content = b'already here'
    dest = tmp_path / 'file.nc'
    dest.write_bytes(content)
    checksum = esgf.compute_checksum(str(dest), 'sha256')

    def fake_urlopen(req, timeout=60):
        raise AssertionError('must not touch the network when checksum already matches')

    monkeypatch.setattr(urllib.request, 'urlopen', fake_urlopen)
    ok = esgf.fetch('http://x/file.nc', str(dest), checksum=checksum)
    assert ok is True
    assert dest.read_bytes() == content


def test_fetch_existing_file_checksum_mismatch_redownloads_full_overwrite(tmp_path, monkeypatch):
    """Checksum mismatch on a pre-existing file unlinks it and redownloads
    from scratch. This also demonstrates that the resume/Range-header path
    is unreachable here: dest_path.exists() is False by the time the
    request is built (just unlinked), so no Range header is ever sent and
    the write mode is always 'wb' (full overwrite) -- old bytes do not
    leak into the new file."""
    old_content = b'STALE-OLD-CONTENT-DOES-NOT-MATCH-CHECKSUM'
    new_content = b'fresh correct bytes'
    dest = tmp_path / 'file.nc'
    dest.write_bytes(old_content)
    checksum = hashlib.sha256(new_content).hexdigest()

    captured = {}

    def fake_urlopen(req, timeout=60):
        captured['has_range'] = req.has_header('Range')
        return _FakeHTTPResponse(new_content)

    monkeypatch.setattr(urllib.request, 'urlopen', fake_urlopen)
    ok = esgf.fetch('http://x/file.nc', str(dest), checksum=checksum)
    assert ok is True
    assert captured['has_range'] is False
    assert dest.read_bytes() == new_content  # not old_content + new_content


def test_fetch_http_error_non_416_raises_runtimeerror(tmp_path, monkeypatch):
    dest = tmp_path / 'file.nc'

    def fake_urlopen(req, timeout=60):
        raise urllib.error.HTTPError('http://x/file.nc', 404, 'Not Found', {}, None)

    monkeypatch.setattr(urllib.request, 'urlopen', fake_urlopen)
    with pytest.raises(RuntimeError, match='Download failed'):
        esgf.fetch('http://x/file.nc', str(dest))


def test_fetch_url_error_raises_runtimeerror(tmp_path, monkeypatch):
    dest = tmp_path / 'file.nc'

    def fake_urlopen(req, timeout=60):
        raise urllib.error.URLError('unreachable')

    monkeypatch.setattr(urllib.request, 'urlopen', fake_urlopen)
    with pytest.raises(RuntimeError, match='Download failed'):
        esgf.fetch('http://x/file.nc', str(dest))


def test_fetch_post_download_checksum_mismatch_raises(tmp_path, monkeypatch):
    dest = tmp_path / 'file.nc'

    def fake_urlopen(req, timeout=60):
        return _FakeHTTPResponse(b'wrong content')

    monkeypatch.setattr(urllib.request, 'urlopen', fake_urlopen)
    with pytest.raises(RuntimeError, match='Checksum mismatch'):
        esgf.fetch('http://x/file.nc', str(dest), checksum='0' * 64)


def test_fetch_existing_file_no_checksum_accepted_without_verification_KNOWN_ISSUE(
        tmp_path, monkeypatch):
    """KNOWN ISSUE (pre-existing, not fixed here -- reported separately).

    When no checksum is supplied, fetch() unconditionally trusts ANY
    pre-existing destination file as "complete" and returns True without
    ever contacting the network -- even if the file is only a truncated
    partial download left over from a previously interrupted run. This is
    also why the documented HTTP-Range 'resume' feature never actually
    resumes anything in this client: by the time the download block would
    run, the file has already been accepted (this branch) or deleted
    (checksum-mismatch branch, see the redownload test above), so
    dest_path.exists() is always False when the Range header would be
    attached.
    """
    partial_content = b'only-part'
    dest = tmp_path / 'file.nc'
    dest.write_bytes(partial_content)

    def fake_urlopen(req, timeout=60):
        raise AssertionError('fetch() should not touch the network here (KNOWN ISSUE)')

    monkeypatch.setattr(urllib.request, 'urlopen', fake_urlopen)
    ok = esgf.fetch('http://x/file.nc', str(dest), checksum=None)
    assert ok is True
    # File is silently accepted as-is -- still just the partial bytes.
    assert dest.read_bytes() == partial_content


def test_fetch_http_416_without_checksum_reports_success_but_creates_no_file_KNOWN_ISSUE(
        tmp_path, monkeypatch):
    """KNOWN ISSUE (pre-existing, not fixed here -- reported separately).

    If urlopen() raises HTTPError(416) on the very first attempt (no
    pre-existing destination file), fetch() takes the "Range not
    satisfiable, file already complete on server" branch (a no-op `pass`)
    and then -- since no checksum was requested to catch the discrepancy --
    returns True immediately. The destination file is never created. The
    caller is told the download "succeeded" even though nothing was ever
    written to disk.
    """
    dest = tmp_path / 'never_created.nc'

    def fake_urlopen(req, timeout=60):
        raise urllib.error.HTTPError('http://x/file.nc', 416, 'Range Not Satisfiable', {}, None)

    monkeypatch.setattr(urllib.request, 'urlopen', fake_urlopen)
    ok = esgf.fetch('http://x/file.nc', str(dest), checksum=None)
    assert ok is True
    assert not dest.exists()
