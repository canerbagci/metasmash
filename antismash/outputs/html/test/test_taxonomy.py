# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Tests for taxonomy TSV file parsing."""

import os
import tempfile

from antismash.outputs.html.taxonomy import parse_taxonomy_file


def _write_temp(content: str) -> str:
    """Write content to a temp file and return the path."""
    fd, filepath = tempfile.mkstemp(suffix=".tsv")
    with os.fdopen(fd, "w", encoding="utf-8") as handle:
        handle.write(content)
    return filepath


def test_basic_parsing():
    filepath = _write_temp(
        "contig_001\tBacteria; Proteobacteria; Gammaproteobacteria\n"
        "contig_002\tArchaea; Euryarchaeota\n"
    )
    try:
        result = parse_taxonomy_file(filepath)
        assert result == {
            "contig_001": "Bacteria; Proteobacteria; Gammaproteobacteria",
            "contig_002": "Archaea; Euryarchaeota",
        }
    finally:
        os.unlink(filepath)


def test_comments_skipped():
    filepath = _write_temp(
        "# This is a comment\n"
        "contig_001\tBacteria\n"
        "# Another comment\n"
        "contig_002\tArchaea\n"
    )
    try:
        result = parse_taxonomy_file(filepath)
        assert len(result) == 2
        assert "contig_001" in result
        assert "contig_002" in result
    finally:
        os.unlink(filepath)


def test_empty_lines_skipped():
    filepath = _write_temp(
        "contig_001\tBacteria\n"
        "\n"
        "\n"
        "contig_002\tArchaea\n"
    )
    try:
        result = parse_taxonomy_file(filepath)
        assert len(result) == 2
    finally:
        os.unlink(filepath)


def test_missing_contigs_not_in_mapping():
    filepath = _write_temp(
        "contig_001\tBacteria\n"
    )
    try:
        result = parse_taxonomy_file(filepath)
        assert "contig_001" in result
        assert "contig_999" not in result
    finally:
        os.unlink(filepath)


def test_whitespace_stripped():
    filepath = _write_temp(
        "  contig_001  \t  Bacteria; Proteobacteria  \n"
    )
    try:
        result = parse_taxonomy_file(filepath)
        assert result == {"contig_001": "Bacteria; Proteobacteria"}
    finally:
        os.unlink(filepath)


def test_malformed_line_skipped():
    filepath = _write_temp(
        "contig_001\tBacteria\n"
        "no_tab_here\n"
        "contig_002\tArchaea\n"
    )
    try:
        result = parse_taxonomy_file(filepath)
        assert len(result) == 2
        assert "no_tab_here" not in result
    finally:
        os.unlink(filepath)


def test_empty_file():
    filepath = _write_temp("")
    try:
        result = parse_taxonomy_file(filepath)
        assert result == {}
    finally:
        os.unlink(filepath)


def test_only_comments():
    filepath = _write_temp(
        "# comment 1\n"
        "# comment 2\n"
    )
    try:
        result = parse_taxonomy_file(filepath)
        assert result == {}
    finally:
        os.unlink(filepath)


def test_tabs_in_lineage():
    filepath = _write_temp(
        "contig_001\tBacteria\textra_col\n"
    )
    try:
        result = parse_taxonomy_file(filepath)
        # maxsplit=1 means everything after first tab is the lineage
        assert result == {"contig_001": "Bacteria\textra_col"}
    finally:
        os.unlink(filepath)
