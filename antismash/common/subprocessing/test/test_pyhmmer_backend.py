# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import os
import pytest

pytest.importorskip("pyhmmer")

from antismash.common.subprocessing import _pyhmmer_backend as backend

DATA = os.path.join(os.path.dirname(__file__), "data_toy.hmm")
FASTA = ">q1\nMAGICKWAND\n"


@pytest.mark.skipif(not os.path.exists(DATA), reason="toy hmm fixture not built")
def test_adapter_exposes_searchio_shape():
    results = backend.run_hmmscan_pyhmmer(DATA, FASTA, opts=[])
    assert isinstance(results, list) and results, "expected >=1 QueryResult-like object"
    qr = results[0]
    assert hasattr(qr, "id") and hasattr(qr, "hsps")
    hsp = qr.hsps[0]
    for attr in ("bitscore", "evalue", "query_id", "hit_id",
                 "hit_description", "query_start", "query_end"):
        assert hasattr(hsp, attr), f"hsp missing {attr}"
    assert hsp.query_id == "q1"
    assert isinstance(hsp.bitscore, float)
    # bitscore is rounded to 1 decimal to match HMMER's hmmer3-text precision
    assert hsp.bitscore == round(hsp.bitscore, 1)
    assert isinstance(hsp.query_start, int)


@pytest.mark.skipif(not os.path.exists(DATA), reason="toy hmm fixture not built")
def test_load_block_caches_by_path():
    backend._BLOCK_CACHE.clear()
    first = backend.load_block(DATA)
    second = backend.load_block(DATA)
    assert first is second  # cached, not reloaded
    assert DATA in backend._BLOCK_CACHE


def test_opts_mapping_and_rejection():
    assert backend._opts_to_pipeline_kwargs(["--cut_tc"]) == {"bias_filter": False, "bit_cutoffs": "trusted"}
    assert backend._opts_to_pipeline_kwargs(["--cut_ga"]) == {"bias_filter": False, "bit_cutoffs": "gathering"}
    assert backend._opts_to_pipeline_kwargs(["-E", "1e-16"]) == {"bias_filter": False, "E": 1e-16}
    assert backend._opts_to_pipeline_kwargs([]) == {"bias_filter": False}
    with pytest.raises(ValueError):
        backend._opts_to_pipeline_kwargs(["--bogus"])


# ---------------------------------------------------------------------------
# hmmsearch adapter tests (query=PROFILE, hit=SEQUENCE)
# ---------------------------------------------------------------------------

# A tiny FASTA with multiple sequences so the search is non-trivial.
# The toy HMM (data_toy.hmm) was built from MAGICKWAND, so q1 should hit.
HSEARCH_FASTA = ">q1\nMAGICKWAND\n>q2\nAAAAAAAAAAAAAAAA\n"


@pytest.mark.skipif(not os.path.exists(DATA), reason="toy hmm fixture not built")
def test_hmmsearch_adapter_shape():
    """run_hmmsearch_pyhmmer returns objects with the SearchIO-compatible shape
    consumed by cluster_prediction.py and structures.HMMerHit.from_hsp()."""
    results = backend.run_hmmsearch_pyhmmer(DATA, HSEARCH_FASTA)
    assert isinstance(results, list), "expected a list"
    assert results, "expected >=1 _HSearchQueryResult"

    # Every result must expose .accession and .hsps (QueryResult-like)
    for qr in results:
        assert hasattr(qr, "accession"), "QueryResult missing .accession"
        assert hasattr(qr, "hsps"), "QueryResult missing .hsps"
        assert isinstance(qr.accession, str), ".accession must be str (may be empty)"

    # Find a result that actually has hits (toy profile matches q1)
    hitting_results = [qr for qr in results if qr.hsps]
    assert hitting_results, "expected at least one QueryResult with hsps"

    qr = hitting_results[0]
    assert qr.hsps, "hsps list should be non-empty"
    hsp = qr.hsps[0]

    # All required HSP attributes (contract for cluster_prediction + from_hsp)
    for attr in ("bitscore", "evalue", "query_id", "hit_id",
                 "query_start", "query_end", "hit_start", "hit_end"):
        assert hasattr(hsp, attr), f"hsp missing required attr: {attr}"

    # query_id is the profile name; hit_id is the sequence name
    assert isinstance(hsp.query_id, str) and hsp.query_id, "query_id should be non-empty str"
    assert isinstance(hsp.hit_id, str) and hsp.hit_id, "hit_id should be non-empty str"

    # The hit on q1 must have hit_id == "q1" (the sequence, not the profile)
    q1_hits = [h for qr in results for h in qr.hsps if h.hit_id == "q1"]
    assert q1_hits, "expected a hit on sequence q1 (toy profile matches MAGICKWAND)"
    h = q1_hits[0]
    assert h.query_id != "q1", "query_id should be the PROFILE name, not the sequence"

    # Numeric types and rounding
    assert isinstance(hsp.bitscore, float), "bitscore must be float"
    assert hsp.bitscore == round(hsp.bitscore, 1), "bitscore must be rounded to 1 decimal"
    assert isinstance(hsp.evalue, float), "evalue must be float"
    assert isinstance(hsp.query_start, int), "query_start must be int"
    assert isinstance(hsp.query_end, int), "query_end must be int"
    assert hsp.query_start >= 0, "query_start must be >=0 (0-based)"
    assert hsp.query_end > hsp.query_start, "query_end must be > query_start"


@pytest.mark.skipif(not os.path.exists(DATA), reason="toy hmm fixture not built")
def test_hmmsearch_empty_sequence_returns_empty():
    """run_hmmsearch_pyhmmer returns [] for empty/whitespace-only input."""
    assert backend.run_hmmsearch_pyhmmer(DATA, "") == []
    assert backend.run_hmmsearch_pyhmmer(DATA, "   \n") == []


@pytest.mark.skipif(not os.path.exists(DATA), reason="toy hmm fixture not built")
def test_hmmsearch_uses_optimized_block_cache():
    """hmmsearch queries are pre-optimised via the shared load_block cache
    (not re-optimised per call), so a hmmsearch run populates _BLOCK_CACHE."""
    backend._BLOCK_CACHE.clear()
    backend.run_hmmsearch_pyhmmer(DATA, HSEARCH_FASTA)
    assert DATA in backend._BLOCK_CACHE  # query profiles cached as an OptimizedProfileBlock


@pytest.mark.skipif(not os.path.exists(DATA), reason="toy hmm fixture not built")
def test_hmmsearch_hsps_are_hashable():
    """Detection (cluster_prediction.filter_results) stores HSPs in sets, so the
    duck objects MUST be hashable — identity-based, like Bio.SearchIO HSPs.
    Regression for `TypeError: unhashable type: '_HSearchHSP'`."""
    results = backend.run_hmmsearch_pyhmmer(DATA, HSEARCH_FASTA)
    hsps = [h for qr in results for h in qr.hsps]
    assert hsps, "expected hits to exercise hashability"
    as_set = set(hsps)  # must not raise TypeError
    assert len(as_set) == len(hsps)  # identity-based: distinct objects stay distinct
