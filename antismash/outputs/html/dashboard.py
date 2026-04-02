# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Dashboard data computation for the metagenomic summary page."""

from collections import defaultdict
from typing import Any, Dict, List, Optional

from antismash.common.layers import RecordLayer, RegionLayer
from antismash.outputs.html.streaming_summary import StreamingRecordSummary, StreamingRegionSummary


def _get_gene_count(region: RegionLayer) -> int:
    """Returns the number of CDS features in a region.

    Uses a pre-cached count when available (set during streaming Phase 2
    before heavy record data is stripped), falling back to the live count
    for the classic (non-streaming) code path.
    """
    cached = getattr(region, '_cached_gene_count', None)
    if cached is not None:
        return cached
    return len(region.region_feature.cds_children)


def _get_most_similar_known(region: RegionLayer) -> str:
    """Returns the description of the most similar known cluster, or empty string."""
    area = region.most_related_area
    if not area.identifier:
        return ""
    similarity = area.similarity_percentage
    if similarity is not None and similarity > 15:
        return area.description
    return ""


def _get_most_similar_url(region: RegionLayer) -> str:
    """Returns the URL of the most similar known cluster, or empty string."""
    area = region.most_related_area
    if not area.identifier:
        return ""
    similarity = area.similarity_percentage
    if similarity is not None and similarity > 15:
        return area.url
    return ""


def _get_similarity_percentage(region: RegionLayer) -> int:
    """Returns the similarity percentage, or -1 if not available."""
    area = region.most_related_area
    if not area.identifier:
        return -1
    val = area.similarity_percentage
    if val is None:
        return -1
    return val


def _get_most_similar_known_from_summary(region: StreamingRegionSummary) -> str:
    """Returns the description of the most similar known cluster, or empty string."""
    area = region.most_related_area
    if not area.identifier:
        return ""
    similarity = area.similarity_percentage
    if similarity is not None and similarity > 15:
        return area.description
    return ""


def _get_most_similar_url_from_summary(region: StreamingRegionSummary) -> str:
    """Returns the URL of the most similar known cluster, or empty string."""
    area = region.most_related_area
    if not area.identifier:
        return ""
    similarity = area.similarity_percentage
    if similarity is not None and similarity > 15:
        return area.url
    return ""


def _get_similarity_percentage_from_summary(region: StreamingRegionSummary) -> int:
    """Returns the similarity percentage, or -1 if not available."""
    area = region.most_related_area
    if not area.identifier:
        return -1
    val = area.similarity_percentage
    if val is None:
        return -1
    return val


def build_dashboard_data(records_with_regions: List[RecordLayer],
                         records_without_regions: List[RecordLayer],
                         taxonomy_mapping: Optional[Dict[str, str]] = None,
                         skipped_record_count: int = 0,
                         ) -> Dict[str, Any]:
    """Build the data dict for the dashboard page.

    Must be called AFTER generate_html_sections(), which populates
    most_related_area on each RegionLayer.

    Arguments:
        records_with_regions: RecordLayer objects that have at least one region
        records_without_regions: RecordLayer objects with no regions
        taxonomy_mapping: optional mapping of contig_id to taxonomy lineage

    Returns:
        A JSON-serializable dict with summary, type/category counts, and table rows
    """
    if taxonomy_mapping is None:
        taxonomy_mapping = {}
    has_taxonomy = bool(taxonomy_mapping)
    total_records = len(records_with_regions) + len(records_without_regions) + skipped_record_count
    records_with_bgcs = len(records_with_regions)
    total_bgcs = 0
    complete_bgcs = 0

    bgc_type_counts: Dict[str, Dict[str, int]] = defaultdict(lambda: {"complete": 0, "on_edge": 0})
    bgc_category_counts: Dict[str, Dict[str, int]] = defaultdict(lambda: {"complete": 0, "on_edge": 0})
    type_to_category: Dict[str, str] = {}
    bgc_table: List[Dict[str, Any]] = []

    for record in records_with_regions:
        for region in record.regions:
            total_bgcs += 1
            region_len = len(region.location)
            gene_count = _get_gene_count(region)

            on_edge = region.contig_edge
            if not on_edge:
                complete_bgcs += 1
            edge_key = "on_edge" if on_edge else "complete"

            for proto in region.region_feature.get_unique_protoclusters():
                type_to_category[proto.product] = proto.product_category

            for product in region.products:
                bgc_type_counts[product][edge_key] += 1

            for category in region.product_categories:
                bgc_category_counts[category][edge_key] += 1

            row: Dict[str, Any] = {
                "record_id": record.seq_record.id,
                "record_index": record.seq_record.record_index,
                "region_number": region.get_region_number(),
                "anchor_id": region.anchor_id,
                "products": ", ".join(sorted(region.products)),
                "product_categories": ", ".join(sorted(region.product_categories)),
                "start": region.start + 1,
                "end": region.end,
                "length": region_len,
                "gene_count": gene_count,
                "contig_edge": on_edge,
                "most_similar_known": _get_most_similar_known(region),
                "most_similar_url": _get_most_similar_url(region),
                "similarity_percentage": _get_similarity_percentage(region),
            }
            if has_taxonomy:
                row["taxonomy"] = taxonomy_mapping.get(record.seq_record.id, "-")
            bgc_table.append(row)

    # Convert defaultdicts to regular dicts for JSON serialization
    type_counts = {k: dict(v) for k, v in sorted(bgc_type_counts.items())}
    category_counts = {k: dict(v) for k, v in sorted(bgc_category_counts.items())}

    # Compute filter metadata for the advanced filter panel
    filter_metadata: Dict[str, Any] = {
        "products_list": sorted(bgc_type_counts.keys()),
        "categories_list": sorted(bgc_category_counts.keys()),
    }
    if bgc_table:
        sim_values = [r["similarity_percentage"] for r in bgc_table
                      if r["similarity_percentage"] >= 0]
        filter_metadata["similarity_min"] = min(sim_values) if sim_values else 0
        filter_metadata["similarity_max"] = max(sim_values) if sim_values else 100
    else:
        filter_metadata.update({
            "similarity_min": 0, "similarity_max": 100,
        })

    return {
        "summary": {
            "total_records": total_records,
            "records_with_bgcs": records_with_bgcs,
            "total_bgcs": total_bgcs,
            "complete_bgcs": complete_bgcs,
        },
        "bgc_type_counts": type_counts,
        "bgc_category_counts": category_counts,
        "type_to_category": type_to_category,
        "bgc_table": bgc_table,
        "has_taxonomy": has_taxonomy,
        "filter_metadata": filter_metadata,
    }


def build_dashboard_data_from_summaries(records_with_regions: List[StreamingRecordSummary],
                                        taxonomy_mapping: Optional[Dict[str, str]] = None,
                                        skipped_record_count: int = 0,
                                        ) -> Dict[str, Any]:
    """Build dashboard data from lightweight streaming summaries."""
    if taxonomy_mapping is None:
        taxonomy_mapping = {}
    has_taxonomy = bool(taxonomy_mapping)
    total_records = len(records_with_regions) + skipped_record_count
    records_with_bgcs = len(records_with_regions)
    total_bgcs = 0
    complete_bgcs = 0

    bgc_type_counts: Dict[str, Dict[str, int]] = defaultdict(lambda: {"complete": 0, "on_edge": 0})
    bgc_category_counts: Dict[str, Dict[str, int]] = defaultdict(lambda: {"complete": 0, "on_edge": 0})
    type_to_category: Dict[str, str] = {}
    bgc_table: List[Dict[str, Any]] = []

    for record in records_with_regions:
        for region in record.regions:
            total_bgcs += 1
            on_edge = region.contig_edge
            if not on_edge:
                complete_bgcs += 1
            edge_key = "on_edge" if on_edge else "complete"

            type_to_category.update(region.product_to_category)

            for product in region.products:
                bgc_type_counts[product][edge_key] += 1

            for category in region.product_categories:
                bgc_category_counts[category][edge_key] += 1

            row: Dict[str, Any] = {
                "record_id": record.id,
                "record_index": record.record_index,
                "region_number": region.get_region_number(),
                "anchor_id": region.anchor_id,
                "products": ", ".join(sorted(region.products)),
                "product_categories": ", ".join(sorted(region.product_categories)),
                "start": region.start + 1,
                "end": region.end,
                "length": region.length,
                "gene_count": region.gene_count,
                "contig_edge": on_edge,
                "most_similar_known": _get_most_similar_known_from_summary(region),
                "most_similar_url": _get_most_similar_url_from_summary(region),
                "similarity_percentage": _get_similarity_percentage_from_summary(region),
            }
            if has_taxonomy:
                row["taxonomy"] = taxonomy_mapping.get(record.id, "-")
            bgc_table.append(row)

    type_counts = {k: dict(v) for k, v in sorted(bgc_type_counts.items())}
    category_counts = {k: dict(v) for k, v in sorted(bgc_category_counts.items())}

    filter_metadata: Dict[str, Any] = {
        "products_list": sorted(bgc_type_counts.keys()),
        "categories_list": sorted(bgc_category_counts.keys()),
    }
    if bgc_table:
        sim_values = [r["similarity_percentage"] for r in bgc_table
                      if r["similarity_percentage"] >= 0]
        filter_metadata["similarity_min"] = min(sim_values) if sim_values else 0
        filter_metadata["similarity_max"] = max(sim_values) if sim_values else 100
    else:
        filter_metadata.update({
            "similarity_min": 0, "similarity_max": 100,
        })

    return {
        "summary": {
            "total_records": total_records,
            "records_with_bgcs": records_with_bgcs,
            "total_bgcs": total_bgcs,
            "complete_bgcs": complete_bgcs,
        },
        "bgc_type_counts": type_counts,
        "bgc_category_counts": category_counts,
        "type_to_category": type_to_category,
        "bgc_table": bgc_table,
        "has_taxonomy": has_taxonomy,
        "filter_metadata": filter_metadata,
    }
