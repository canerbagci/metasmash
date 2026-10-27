# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Contig-to-bin TSV file parsing and per-bin summaries for the Genome Bins page."""

from collections import Counter
import logging
from typing import Any, Callable, Dict, Iterable, List, Optional

UNBINNED_LABEL = "Unbinned"


def parse_bins_file(filepath: str) -> Dict[str, str]:
    """Parse a two-column TSV file mapping contig IDs to genome bin names.

    Format: contig_id<TAB>bin_name
    Lines starting with '#' are treated as comments and skipped.

    Arguments:
        filepath: path to the TSV file

    Returns:
        a dictionary mapping contig_id to bin name
    """
    mapping: Dict[str, str] = {}
    duplicates = 0
    with open(filepath, encoding="utf-8") as handle:
        for line_num, line in enumerate(handle, 1):
            line = line.rstrip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) != 2:
                logging.warning("Bins file line %d: expected 2 tab-separated columns, got %d; skipping",
                                line_num, len(parts))
                continue
            contig_id, bin_name = parts[0].strip(), parts[1].strip()
            if not contig_id or not bin_name:
                logging.warning("Bins file line %d: empty contig ID or bin name; skipping", line_num)
                continue
            if contig_id in mapping:
                duplicates += 1
                continue
            mapping[contig_id] = bin_name
    if duplicates:
        logging.warning("Bins file: %d contig(s) listed more than once; kept the first bin for each",
                        duplicates)
    return mapping


def get_bin(bins_mapping: Dict[str, str], record_id: str, original_id: Optional[str]) -> str:
    """Returns the bin for a record, trying its ID first and then the ID it
    had before antiSMASH renamed it, or an empty string if it is unbinned.
    """
    if record_id in bins_mapping:
        return bins_mapping[record_id]
    if original_id:
        return bins_mapping.get(original_id, "")
    return ""


def count_unmatched_contigs(bins_mapping: Dict[str, str], known_ids: Iterable[str],
                            normalise: Optional[Callable[[str], str]] = None) -> int:
    """Returns the number of contigs in the bin mapping that are not in known_ids.

    Arguments:
        bins_mapping: mapping of contig_id to bin name
        known_ids: record IDs of the input
        normalise: optional function giving the form a contig ID takes in
                   known_ids, e.g. with invalid characters stripped

    Returns:
        the number of unmatched contigs
    """
    known = set(known_ids)
    unmatched = 0
    for contig_id in bins_mapping:
        if contig_id in known or (normalise is not None and normalise(contig_id) in known):
            continue
        unmatched += 1
    return unmatched


def load_bins_mapping(filepath: str, known_ids: Iterable[str],
                      normalise: Optional[Callable[[str], str]] = None) -> Optional[Dict[str, str]]:
    """Parse the bins file, if one was given, warning about any contigs in it
    that are not among the input's record IDs.

    Arguments:
        filepath: path to the TSV file, or an empty string if none was given
        known_ids: record IDs of the input, both current and original
        normalise: optional function giving the form a contig ID takes in known_ids

    Returns:
        a dictionary mapping contig_id to bin name, or None if no file was given
    """
    if not filepath:
        return None
    mapping = parse_bins_file(filepath)
    unmatched = count_unmatched_contigs(mapping, known_ids, normalise)
    if unmatched:
        logging.warning("%d contigs in the bin file were not found in the input "
                        "(shorter than --minlength or absent)", unmatched)
    return mapping


def _bin_row(name: str, contigs: Optional[int], rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    rows = sorted(rows, key=lambda row: (row["record_index"], row["region_number"]))
    category_counts: Counter = Counter()
    for row in rows:
        category_counts.update(category for category in row["product_categories"].split(", ") if category)
    on_edge = sum(1 for row in rows if row["contig_edge"])
    return {
        "bin": name,
        "contigs": contigs,
        "bgcs": len(rows),
        "complete": len(rows) - on_edge,
        "on_edge": on_edge,
        "known_hits": sum(1 for row in rows if row["most_similar_known"]),
        "categories": [[category, count] for category, count in
                       sorted(category_counts.items(), key=lambda item: (-item[1], item[0]))],
        "anchors": [row["anchor_id"] for row in rows],
    }


def build_bins_page_data(bgc_table: List[Dict[str, Any]],
                         bins_mapping: Dict[str, str]) -> Dict[str, Any]:
    """Build the data dict for the Genome Bins page.

    Arguments:
        bgc_table: the dashboard's BGC rows, each with a "bin" value
                   (empty for BGCs on unbinned contigs)
        bins_mapping: mapping of contig_id to bin name

    Returns:
        A JSON-serializable dict with summary counts, one row per bin
        (including bins without BGCs) and a row for unbinned BGCs
    """
    contig_counts = Counter(bins_mapping.values())
    rows_by_bin: Dict[str, List[Dict[str, Any]]] = {name: [] for name in contig_counts}
    unbinned_rows: List[Dict[str, Any]] = []
    for row in bgc_table:
        if row["bin"]:
            rows_by_bin.setdefault(row["bin"], []).append(row)
        else:
            unbinned_rows.append(row)

    bins = [_bin_row(name, contig_counts.get(name, 0), rows) for name, rows in rows_by_bin.items()]
    bins.sort(key=lambda bin_row: (-bin_row["bgcs"], bin_row["bin"]))

    bgcs_in_bins = sum(bin_row["bgcs"] for bin_row in bins)
    return {
        "summary": {
            "bins_imported": len(contig_counts),
            "bins_with_bgcs": sum(1 for bin_row in bins if bin_row["bgcs"]),
            "bgcs_in_bins": bgcs_in_bins,
            "bgcs_unbinned": len(unbinned_rows),
        },
        "bins": bins,
        "unbinned": _bin_row(UNBINNED_LABEL, None, unbinned_rows),
    }
