# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Lightweight summaries for streaming HTML finalization."""

from dataclasses import dataclass, field
from typing import Dict, List, Optional

from antismash.common.html_renderer import Markup
from antismash.common.layers import AbstractRelatedArea, RecordLayer, RegionLayer
from antismash.outputs.html import js


@dataclass(frozen=True)
class StreamingRelatedAreaSummary:
    """A JSON/template-friendly snapshot of a related area."""

    identifier: str = ""
    description: str = ""
    product: str = ""
    similarity_percentage: Optional[int] = None
    url: str = ""

    @classmethod
    def from_related_area(cls, area: AbstractRelatedArea) -> "StreamingRelatedAreaSummary":
        """Build a summary from a RegionLayer related-area object."""
        return cls(
            identifier=area.identifier,
            description=area.description,
            product=area.product,
            similarity_percentage=area.similarity_percentage,
            url=area.url,
        )


@dataclass(frozen=True)
class StreamingRegionSummary:
    """A lightweight representation of a region for overview/dashboard rendering."""

    anchor_id: str
    region_number: int
    start: int
    end: int
    length: int
    products: List[str]
    product_categories: List[str]
    css_class: str
    contig_edge: bool
    gene_count: int
    product_to_category: Dict[str, str]
    most_related_area: StreamingRelatedAreaSummary = field(default_factory=StreamingRelatedAreaSummary)

    @classmethod
    def from_region_layer(cls, region: RegionLayer) -> "StreamingRegionSummary":
        """Build a summary from a fully populated RegionLayer."""
        product_to_category = {
            proto.product: proto.product_category
            for proto in region.region_feature.get_unique_protoclusters()
        }
        return cls(
            anchor_id=region.anchor_id,
            region_number=region.get_region_number(),
            start=region.start,
            end=region.end,
            length=len(region.location),
            products=list(region.products),
            product_categories=sorted(region.product_categories),
            css_class=js.get_region_css(region.region_feature),
            contig_edge=region.contig_edge,
            gene_count=len(region.region_feature.cds_children),
            product_to_category=product_to_category,
            most_related_area=StreamingRelatedAreaSummary.from_related_area(region.most_related_area),
        )

    def get_region_number(self) -> int:
        """Match the RegionLayer interface used by templates."""
        return self.region_number


@dataclass(frozen=True)
class StreamingRecordSummary:
    """A lightweight representation of a record for overview rendering."""

    id: str
    record_index: int
    original_id: str
    source: str
    combined_source_count: int
    regions: List[StreamingRegionSummary]

    @classmethod
    def from_record_layer(cls, record: RecordLayer) -> "StreamingRecordSummary":
        """Build a summary from a transient RecordLayer."""
        seq_record = record.seq_record
        source = seq_record.annotations.get("source", "")
        combined_source_count = len(seq_record.get_sources()) if seq_record.has_multiple_sources() else 1
        return cls(
            id=seq_record.id,
            record_index=seq_record.record_index or 0,
            original_id=seq_record.original_id or "",
            source=source,
            combined_source_count=combined_source_count,
            regions=[StreamingRegionSummary.from_region_layer(region) for region in record.regions],
        )

    def get_name(self) -> str:
        """Match the RecordLayer interface used by templates."""
        name = self.id
        if self.combined_source_count > 1:
            name += (
                f" (combined with {self.combined_source_count - 1} "
                f"other{'s' if self.combined_source_count > 2 else ''})"
            )
        return name

    def get_from_record(self) -> Markup:
        """Match the RecordLayer interface used by templates."""
        current_id = f"<strong>{self.get_name()}</strong>"
        if self.original_id:
            if len(self.original_id) < 40:
                current_id += f" (original name was: {self.original_id})"
            else:
                current_id += f" (original name was: {self.original_id[:60]}...)"
        source = f" ({self.source})" if self.source else ""
        return Markup(f"{current_id}{source}")
