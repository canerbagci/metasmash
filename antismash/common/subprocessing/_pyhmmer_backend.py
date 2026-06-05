# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""In-process pyhmmer backend for hmmscan and hmmsearch.

Provides the same result shape as Bio.SearchIO so callers that consume SearchIO
QueryResult / HSP objects work unchanged.

Public surface (consumed by downstream tasks):
  _BLOCK_CACHE           : dict[str, tuple[OptimizedProfileBlock, int]]
  _opts_to_pipeline_kwargs(opts) -> dict
  load_block(db_path)    -> tuple[OptimizedProfileBlock, int]  (used for both hmmscan
                           profile DBs and hmmsearch query profiles, pre-optimised once)
  run_hmmscan_pyhmmer(target_hmmfile, query_sequence, opts=None) -> list[_QueryResult]
  run_hmmsearch_pyhmmer(query_hmmfile, target_sequence, ...) -> list[_HSearchQueryResult]

"""

from __future__ import annotations

import io
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

@dataclass(eq=False)
class _HSP:
    """Single reported domain, shaped like a Bio.SearchIO HSP."""
    bitscore: float
    evalue: float
    query_id: str
    hit_id: str
    hit_description: str
    query_start: int
    query_end: int


@dataclass(eq=False)
class _QueryResult:
    """Per-query results, shaped like a Bio.SearchIO QueryResult."""
    id: str
    hsps: List[_HSP] = field(default_factory=list)


@dataclass(eq=False)
class _HSearchHSP:
    """Single reported domain from hmmsearch, shaped like a Bio.SearchIO HSP.

    hmmsearch is inverted vs hmmscan: query=PROFILE, hit=SEQUENCE.
    Attribute names match what cluster_prediction.py, structures.py, and the
    product-analysis modules (lanthipeptides/sactipeptides/thiopeptides) read:
      - query_id  : profile name  (top_hits.query.name)
      - hit_id    : sequence name (hit.name)
      - query_start/query_end : profile (HMM) coords, 0-based half-open
                                (aln.hmm_from - 1 / aln.hmm_to)
      - hit_start/hit_end     : sequence coords, 0-based half-open
                                (aln.target_from - 1 / aln.target_to) — read by
                                sactipeptides' radical-SAM position check
      - evalue    : domain.i_evalue
      - bitscore  : round(domain.score, 1)
    """
    bitscore: float
    evalue: float
    query_id: str
    hit_id: str
    query_start: int
    query_end: int
    hit_start: int
    hit_end: int


@dataclass(eq=False)
class _HSearchQueryResult:
    """Per-query (per-profile) results from hmmsearch, shaped like Bio.SearchIO QueryResult.

    Attributes:
      accession : profile accession string (empty string if absent), consumed by
                  cluster_prediction.py as ``runresult.accession.split('.')[0]``.
      hsps      : list of _HSearchHSP objects for all reported domains.
    """
    accession: str
    hsps: List[_HSearchHSP] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Option mapping
# ---------------------------------------------------------------------------

def _opts_to_pipeline_kwargs(opts: List[str]) -> dict:
    """Map hmmscan option strings to pyhmmer Pipeline keyword args.

    Handled options:
      --cut_tc   → bit_cutoffs="trusted"
      --cut_ga   → bit_cutoffs="gathering"
      -E <float> → E=float(value)

    bias_filter=False is always included (mirrors run_hmmscan --nobias).
    Raises ValueError for any unrecognised option.
    """
    kwargs: dict = {"bias_filter": False}
    i = 0
    while i < len(opts):
        opt = opts[i]
        if opt == "--cut_tc":
            kwargs["bit_cutoffs"] = "trusted"
        elif opt == "--cut_ga":
            kwargs["bit_cutoffs"] = "gathering"
        elif opt == "-E":
            i += 1
            if i >= len(opts):
                raise ValueError("-E requires a following value")
            kwargs["E"] = float(opts[i])
        else:
            raise ValueError(f"Unmapped hmmscan option: {opt!r}")
        i += 1
    return kwargs


# ---------------------------------------------------------------------------
# Profile block cache
# ---------------------------------------------------------------------------

_BLOCK_CACHE: dict = {}  # db_path -> (OptimizedProfileBlock, n_profiles)


def load_block(db_path: str) -> Tuple:
    """Return (OptimizedProfileBlock, n_profiles) for db_path, cached after first load.

    Profiles are stored as an OptimizedProfileBlock.

    When the HMM database is pressed (the .h3m/.h3f files exist), the optimised
    profiles are read DIRECTLY from the pressed files.
    Falls back to optimising from raw HMMs for unpressed files.
    """
    if db_path in _BLOCK_CACHE:
        return _BLOCK_CACHE[db_path]

    import pyhmmer.plan7
    from pyhmmer.easel import Alphabet

    abc = Alphabet.amino()
    block = pyhmmer.plan7.OptimizedProfileBlock(abc)
    with pyhmmer.plan7.HMMFile(db_path) as handle:
        if handle.is_pressed():
            for optimized in handle.optimized_profiles():
                block.append(optimized)
        else:
            background = pyhmmer.plan7.Background(abc)
            for hmm in handle:
                block.append(hmm.to_profile(background).to_optimized())

    n_profiles = len(block)
    result = (block, n_profiles)
    _BLOCK_CACHE[db_path] = result
    return result


# ---------------------------------------------------------------------------
# Main scanning entry point
# ---------------------------------------------------------------------------

def run_hmmscan_pyhmmer(
    target_hmmfile: str,
    query_sequence: str,
    opts: Optional[List[str]] = None,
) -> List[_QueryResult]:
    """Scan query_sequence against target_hmmfile using pyhmmer in-process.

    Args:
        target_hmmfile:  Path to the .hmm database (pressed or plain).
        query_sequence:  FASTA-format string (one or more sequences).
        opts:            List of hmmscan option strings (see _opts_to_pipeline_kwargs).

    Returns:
        List of _QueryResult objects in query input order, containing only
        reported hits/domains (mirrors the binary's text output).

    Raises:
        ValueError: if query_sequence is empty (mirrors run_hmmscan behaviour).
    """
    if not query_sequence or not query_sequence.strip():
        raise ValueError("run_hmmscan_pyhmmer: query_sequence is empty")

    import pyhmmer.hmmer
    from pyhmmer.easel import Alphabet, SequenceFile

    if opts is None:
        opts = []

    abc = Alphabet.amino()
    pipeline_kwargs = _opts_to_pipeline_kwargs(opts)
    block, n_profiles = load_block(target_hmmfile)

    # Parse query sequences from FASTA text
    fasta_bytes = query_sequence.encode() if isinstance(query_sequence, str) else query_sequence
    with SequenceFile(io.BytesIO(fasta_bytes), format="fasta",
                      digital=True, alphabet=abc) as sf:
        queries = sf.read_block()

    query_results: List[_QueryResult] = []

    for top_hits in pyhmmer.hmmer.hmmscan(
        queries, block, cpus=1, Z=n_profiles, **pipeline_kwargs
    ):
        raw_qname = top_hits.query.name
        query_name: str = raw_qname.decode() if isinstance(raw_qname, bytes) else raw_qname

        qr = _QueryResult(id=query_name)

        for hit in top_hits:
            if not hit.reported:
                continue

            raw_pname = hit.name
            hit_id: str = raw_pname.decode() if isinstance(raw_pname, bytes) else raw_pname

            raw_desc = hit.description
            if raw_desc is None:
                hit_description = ""
            elif isinstance(raw_desc, bytes):
                hit_description = raw_desc.decode()
            else:
                hit_description = raw_desc

            for domain in hit.domains:
                if not domain.reported:
                    continue
                aln = domain.alignment
                hsp = _HSP(
                    # round to 1 decimal to match HMMER's hmmer3-text output precision
                    bitscore=round(float(domain.score), 1),
                    # 2 sig figs to match HMMER's text/domtab evalue precision
                    evalue=float(f"{float(domain.i_evalue):.2g}"),
                    query_id=query_name,
                    hit_id=hit_id,
                    hit_description=hit_description,
                    query_start=int(aln.target_from) - 1,  # 1-based → 0-based
                    query_end=int(aln.target_to),
                )
                qr.hsps.append(hsp)

        query_results.append(qr)

    return query_results


# ---------------------------------------------------------------------------
# hmmsearch entry point (query=PROFILE, hit=SEQUENCE — inverted vs hmmscan)
# ---------------------------------------------------------------------------

def run_hmmsearch_pyhmmer(
    query_hmmfile: str,
    target_sequence: str,
    use_tempfile: bool = False,
    force_single_core: bool = False,
) -> List[_HSearchQueryResult]:
    """Search query_hmmfile profiles against target_sequence using pyhmmer in-process.

    Mirrors the subprocess path in run_hmmsearch: profiles are queries, sequences
    are targets.

    Args:
        query_hmmfile:   Path to the .hmm file (profiles are the queries).
        target_sequence: FASTA-format string of target sequences to search.
        use_tempfile:    Accepted for API parity with the subprocess path; ignored
        force_single_core: Accepted for API parity; ignored (pyhmmer runs with
                         cpus=1 unconditionally to stay thread-safe in workers).

    Returns:
        List of _HSearchQueryResult objects, one per profile, ordered by profile
        input order.  Only reported hits/domains are included.
    """
    if not target_sequence or not target_sequence.strip():
        return []

    import pyhmmer.hmmer
    from pyhmmer.easel import Alphabet, SequenceFile

    abc = Alphabet.amino()

    # Parse target sequences from FASTA text
    fasta_bytes = (target_sequence.encode()
                   if isinstance(target_sequence, str) else target_sequence)
    with SequenceFile(io.BytesIO(fasta_bytes), format="fasta",
                      digital=True, alphabet=abc) as sf:
        sequences = sf.read_block()

    n_seqs = len(sequences)
    if n_seqs == 0:
        return []

    # Load (and cache) the query profiles, PRE-OPTIMISED once via the shared
    # OptimizedProfileBlock cache (same machinery + config as hmmscan).
    block, _n_profiles = load_block(query_hmmfile)

    query_results: List[_HSearchQueryResult] = []

    for top_hits in pyhmmer.hmmer.hmmsearch(block, sequences, cpus=1, Z=n_seqs):
        # top_hits.query is the HMM (profile) that was searched
        raw_pname = top_hits.query.name
        profile_name: str = (raw_pname.decode()
                              if isinstance(raw_pname, bytes) else raw_pname)

        raw_acc = top_hits.query.accession
        if raw_acc is None:
            accession = ""
        elif isinstance(raw_acc, bytes):
            accession = raw_acc.decode()
        else:
            accession = raw_acc

        qr = _HSearchQueryResult(accession=accession)

        for hit in top_hits:
            if not hit.reported:
                continue

            raw_sname = hit.name
            seq_name: str = (raw_sname.decode()
                              if isinstance(raw_sname, bytes) else raw_sname)

            for domain in hit.domains:
                if not domain.reported:
                    continue
                aln = domain.alignment
                # Coordinate mapping: SearchIO hmmsearch3-domtab (hmm_as_hit=False)
                # exposes HMM coords as query_start/query_end after its hmm<-->ali
                # swap.  pyhmmer: aln.hmm_from/hmm_to are 1-based profile coords.
                hsp = _HSearchHSP(
                    bitscore=round(float(domain.score), 1),
                    # 2 sig figs to match HMMER's text/domtab evalue precision
                    evalue=float(f"{float(domain.i_evalue):.2g}"),
                    query_id=profile_name,
                    hit_id=seq_name,
                    query_start=int(aln.hmm_from) - 1,  # profile coords, 1-based → 0-based
                    query_end=int(aln.hmm_to),
                    hit_start=int(aln.target_from) - 1,  # sequence coords, 1-based → 0-based
                    hit_end=int(aln.target_to),
                )
                qr.hsps.append(hsp)

        query_results.append(qr)

    return query_results
