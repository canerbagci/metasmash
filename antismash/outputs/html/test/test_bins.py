# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Tests for contig-to-bin TSV parsing and per-bin summaries."""

import logging
import os
import tempfile

from antismash.outputs.html.bins import (
    UNBINNED_LABEL,
    build_bins_page_data,
    count_unmatched_contigs,
    get_bin,
    load_bins_mapping,
    parse_bins_file,
)


def _write_temp(content: str) -> str:
    """Write content to a temp file and return the path."""
    fd, filepath = tempfile.mkstemp(suffix=".tsv")
    with os.fdopen(fd, "w", encoding="utf-8", newline="") as handle:
        handle.write(content)
    return filepath


def _parse(content: str) -> dict:
    filepath = _write_temp(content)
    try:
        return parse_bins_file(filepath)
    finally:
        os.unlink(filepath)


def _row(bin_name, record_index, region_number, products="NRPS", categories="NRPS",
         contig_edge=False, known=""):
    return {
        "bin": bin_name,
        "anchor_id": f"r{record_index}c{region_number}",
        "record_index": record_index,
        "region_number": region_number,
        "products": products,
        "product_categories": categories,
        "contig_edge": contig_edge,
        "most_similar_known": known,
    }


def test_basic_parsing():
    assert _parse("c1\tbin.1\nc2\tbin.2\n") == {"c1": "bin.1", "c2": "bin.2"}


def test_comments_and_empty_lines_skipped():
    assert _parse("# header\nc1\tbin.1\n\n\nc2\tbin.2\n") == {"c1": "bin.1", "c2": "bin.2"}


def test_malformed_line_skipped(caplog):
    with caplog.at_level(logging.WARNING):
        result = _parse("c1\tbin.1\nc3\nc2\tbin.2\n")
    assert result == {"c1": "bin.1", "c2": "bin.2"}
    assert "line 2" in caplog.text


def test_crlf_and_trailing_whitespace():
    assert _parse("c1\tbin.1 \r\nc2\t bin.2\t\r\n") == {"c1": "bin.1", "c2": "bin.2"}


def test_duplicate_keeps_first(caplog):
    with caplog.at_level(logging.WARNING):
        result = _parse("c1\tA\nc1\tB\n")
    assert result == {"c1": "A"}
    assert "1 contig" in caplog.text


def test_get_bin_prefers_record_id():
    mapping = {"new": "A", "old": "B"}
    assert get_bin(mapping, "new", "old") == "A"
    assert get_bin(mapping, "x", "old") == "B"
    assert get_bin(mapping, "x", None) == ""
    assert get_bin(mapping, "x", "") == ""


def test_count_unmatched():
    mapping = {"c1": "A", "c2": "A", "c3": "B"}
    assert count_unmatched_contigs(mapping, ["c1", "zz"]) == 2


def test_page_data_counts():
    mapping = {"c1": "A", "c2": "A", "c3": "B", "c4": "C"}
    rows = [
        _row("A", 1, 1, products="NRPS, T1PKS", categories="NRPS, PKS", contig_edge=True),
        _row("A", 1, 2, products="NRPS", known="x"),
        _row("B", 3, 1, products="terpene", categories="terpene"),
        _row("", 5, 1, products="terpene", categories="terpene"),
    ]
    data = build_bins_page_data(rows, mapping)
    assert data["summary"] == {
        "bins_imported": 3,
        "bins_with_bgcs": 2,
        "bgcs_in_bins": 3,
        "bgcs_unbinned": 1,
    }
    first = data["bins"][0]
    assert first["bin"] == "A"
    assert first["contigs"] == 2
    assert first["bgcs"] == 2
    assert first["complete"] == 1
    assert first["on_edge"] == 1
    assert first["known_hits"] == 1
    assert first["categories"] == [["NRPS", 2], ["PKS", 1]]
    assert "types" not in first
    assert data["bins"][-1]["bin"] == "C"
    assert data["bins"][-1]["bgcs"] == 0
    assert data["unbinned"]["bin"] == UNBINNED_LABEL
    assert data["unbinned"]["contigs"] is None
    assert data["unbinned"]["bgcs"] == 1


def test_page_data_anchor_order():
    mapping = {"c1": "A", "c2": "A"}
    rows = [_row("A", 2, 1), _row("A", 1, 3), _row("A", 1, 1)]
    data = build_bins_page_data(rows, mapping)
    assert data["bins"][0]["anchors"] == ["r1c1", "r1c3", "r2c1"]


def test_page_data_empty_mapping():
    data = build_bins_page_data([_row("", 1, 1)], {})
    assert data["bins"] == []
    assert data["summary"]["bins_imported"] == 0
    assert data["unbinned"]["bgcs"] == 1


def test_page_data_keeps_raw_names():
    name = 'b<1>"x'
    data = build_bins_page_data([_row(name, 1, 1)], {"c1": name})
    assert data["bins"][0]["bin"] == name


def test_load_bins_mapping_none_without_path():
    assert load_bins_mapping("", ["c1"]) is None


def test_load_bins_mapping_warns_unmatched(caplog):
    filepath = _write_temp("c1\tA\nmissing\tA\n")
    try:
        with caplog.at_level(logging.WARNING):
            assert load_bins_mapping(filepath, ["c1"]) == {"c1": "A", "missing": "A"}
    finally:
        os.unlink(filepath)
    assert "1 contigs in the bin file were not found in the input" in caplog.text


def test_load_bins_mapping_no_warning_when_all_match(caplog):
    filepath = _write_temp("c1\tA\n")
    try:
        with caplog.at_level(logging.WARNING):
            load_bins_mapping(filepath, ["c1"])
    finally:
        os.unlink(filepath)
    assert "not found" not in caplog.text


def test_extra_columns_skipped(caplog):
    with caplog.at_level(logging.WARNING):
        result = _parse("c1\tbin.1\t0.97\nc2\tbin.2\n")
    assert result == {"c2": "bin.2"}
    assert "line 1" in caplog.text


def test_count_unmatched_with_normaliser():
    mapping = {"NODE|1": "A", "gone": "B"}
    assert count_unmatched_contigs(mapping, ["NODE1"], normalise=lambda s: s.replace("|", "")) == 1
