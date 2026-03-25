# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Taxonomy TSV file parsing for optional dashboard taxonomy column."""

import logging
from typing import Dict


def parse_taxonomy_file(filepath: str) -> Dict[str, str]:
    """Parse a two-column TSV file mapping contig IDs to taxonomy lineages.

    Format: contig_id<TAB>taxonomy_lineage
    Lines starting with '#' are treated as comments and skipped.
    Empty lines are skipped.

    Arguments:
        filepath: path to the TSV file

    Returns:
        a dictionary mapping contig_id to taxonomy lineage string
    """
    mapping: Dict[str, str] = {}
    with open(filepath, encoding="utf-8") as handle:
        for line_num, line in enumerate(handle, 1):
            line = line.rstrip("\n\r")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t", maxsplit=1)
            if len(parts) != 2:
                logging.warning("Taxonomy file line %d: expected 2 tab-separated columns, got %d; skipping",
                                line_num, len(parts))
                continue
            contig_id, lineage = parts
            contig_id = contig_id.strip()
            lineage = lineage.strip()
            if contig_id:
                mapping[contig_id] = lineage
    return mapping
