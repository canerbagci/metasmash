# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Responsible for creating the single web page results """

import base64
import gzip
import importlib
import pkgutil
import string
import os
from typing import cast, Any, Dict, IO, List, Tuple, Optional

from antismash.common import json, path
from antismash.common.html_renderer import (
    FileTemplate,
    HTMLSections,
    docs_link,
    get_antismash_js_version,
    get_antismash_js_url,
)
from antismash.common.layers import RecordLayer, RegionLayer, OptionsLayer
from antismash.common.module_results import ModuleResults
from antismash.common.secmet import Record
from antismash.config import ConfigType
from antismash.modules import tfbs_finder as tfbs, tta
from antismash.outputs.html import js
from antismash.custom_typing import AntismashModule, VisualisationModule

from .visualisers import gene_table

TEMPLATE_PATH = path.get_full_path(__file__, "templates")


def _strip_leading_whitespace(text: str) -> str:
    """Strip leading spaces from each line — faster than regex for large strings."""
    return "\n".join(line.lstrip(" ") for line in text.split("\n"))


def _get_visualisers() -> List[VisualisationModule]:
    """ Gather all the visualisation-only submodules """
    modules = []
    for module_data in pkgutil.walk_packages([path.get_full_path(__file__, "visualisers")]):
        module = importlib.import_module(f"antismash.outputs.html.visualisers.{module_data.name}")
        assert hasattr(module, "has_enough_results"), f"bad visualisation module: {module_data.name}"
        modules.append(cast(VisualisationModule, module))
    return modules


VISUALISERS = _get_visualisers()


def build_json_data(records: List[Record], results: List[Dict[str, ModuleResults]],
                    options: ConfigType, all_modules: List[AntismashModule]) -> Tuple[
                        List[Dict[str, Any]],
                        dict[str, dict[str, json.JSONCompatible]],
                    ]:
    """ Builds JSON versions of records and domains for use in drawing SVGs with
        javascript.

        Arguments:
            records: a list of Records to convert
            results: a dictionary mapping record id to a list of ModuleResults to convert
            options: antiSMASH options

        Returns:
            a tuple of
                a list of JSON-friendly dicts representing records
                a list of JSON-friendly dicts representing domains
    """

    js_records = js.convert_records(records, results, options)

    js_results: dict[str, dict[str, json.JSONCompatible]] = {}

    for i, record in enumerate(records):
        json_record = js_records[i]
        json_record['seq_id'] = "".join(char for char in json_record['seq_id'] if char in string.printable)
        for region, json_region in zip(record.get_regions(), json_record['regions']):
            handlers = find_plugins_for_cluster(all_modules, json_region)
            region_results = {}
            for handler in handlers:
                # if there's no results for the module, don't let it try
                if handler.__name__ not in results[i]:
                    continue
                if hasattr(handler, "generate_javascript_data"):
                    data = handler.generate_javascript_data(record, region, results[i][handler.__name__])
                    region_results[handler.__name__] = data

            for aggregator in VISUALISERS:
                if not hasattr(aggregator, "generate_javascript_data"):
                    continue
                if aggregator.has_enough_results(record, region, results[i]):
                    data = aggregator.generate_javascript_data(record, region, results[i])
                    region_results[aggregator.__name__] = data

            if region_results:
                js_results[RegionLayer.build_anchor_id(region)] = region_results

    return js_records, js_results


_LIGHTWEIGHT_REGION_KEYS = {"anchor", "start", "end", "idx", "type", "products",
                            "product_categories", "cssClass"}


def _make_lightweight_region(region: Dict[str, Any]) -> Dict[str, Any]:
    """ Returns a copy of a region dict containing only the keys needed for
        the index page (no orfs, clusters, sites, or sources).
    """
    return {k: v for k, v in region.items() if k in _LIGHTWEIGHT_REGION_KEYS}


def _build_regions_index(records: List[Dict[str, Any]],
                         ) -> Tuple[List[Dict[str, Any]], Dict[str, Any], Dict[str, Dict[str, Any]]]:
    """ Builds lightweight record copies, the regions index, and collects
        per-region full JSON data.

        Returns:
            a tuple of
                lightweight records (for the index file),
                regions_index dict,
                per-region full JSON keyed by anchor
    """
    lightweight_records = []
    regions_index: Dict[str, Any] = {"order": []}
    per_region_json: Dict[str, Dict[str, Any]] = {}

    for record in records:
        if not record["regions"]:
            continue  # skip records with no regions to reduce regions.js size
        light_record = {k: v for k, v in record.items() if k != "regions"}
        light_record["regions"] = []
        for region in record["regions"]:
            anchor = region["anchor"]
            light_region = _make_lightweight_region(region)
            light_record["regions"].append(light_region)
            regions_index[anchor] = light_region
            regions_index["order"].append(anchor)
            per_region_json[anchor] = region
        lightweight_records.append(light_record)

    return lightweight_records, regions_index, per_region_json


def _write_regions_index(lightweight_records: List[Dict[str, Any]],
                         regions_index: Dict[str, Any],
                         output_dir: str) -> None:
    """ Writes the lightweight regions.js index file. """
    with open(os.path.join(output_dir, "regions.js"), "w", encoding="utf-8") as handle:
        handle.write(f"var recordData = {json.dumps(lightweight_records)};\n")
        handle.write(f"var all_regions = {json.dumps(regions_index)};\n")
        handle.write("var resultsData = {};\n")
        handle.write("var regionHTML = {};\n")
        handle.write("var lazyLoading = true;\n")


class StreamingRegionDataWriter:
    """ Writes all per-region data into a single compressed JS file.

        Each region's data (region JSON, module results, rendered HTML) is
        individually gzip-compressed, base64-encoded, and stored as a
        property of a JS object. This allows the browser to decompress
        individual regions on demand without loading everything into memory.
    """
    def __init__(self, handle: IO[str]) -> None:
        self._handle = handle
        self._first = True
        handle.write("var _regionsData = {\n")

    def write_region(self, anchor: str, region_json: Dict[str, Any],
                     results: dict[str, json.JSONCompatible],
                     html: str) -> None:
        """ Write a single region's data to the compressed file. """
        payload = json.dumps({
            "region": region_json,
            "results": results,
            "html": html,
        })
        compressed = gzip.compress(payload.encode("utf-8"), compresslevel=6)
        b64 = base64.b64encode(compressed).decode("ascii")
        if not self._first:
            self._handle.write(",\n")
        self._first = False
        self._handle.write(f"  {json.dumps(anchor)}: \"{b64}\"")

    def finalize(self) -> None:
        """ Close the JS object literal. """
        self._handle.write("\n};\n")


def _write_compressed_data_file(output_dir: str,
                                per_region_json: Dict[str, Dict[str, Any]],
                                module_results: dict[str, dict[str, json.JSONCompatible]],
                                region_html: Dict[str, str]) -> None:
    """ Writes a single regions_data.js containing all region data compressed.

        Used by the classic (non-streaming) path as an alternative to
        writing thousands of individual per-region .js files.
    """
    data_path = os.path.join(output_dir, "regions_data.js")
    with open(data_path, "w", encoding="utf-8") as handle:
        writer = StreamingRegionDataWriter(handle)
        for anchor, region_json in per_region_json.items():
            results = module_results.get(anchor, {})
            html = region_html.get(anchor, "")
            writer.write_region(anchor, region_json, results, html)
        writer.finalize()


def _write_single_region_file(regions_dir: str, anchor: str,
                               region_json: Dict[str, Any],
                               module_results: dict[str, dict[str, json.JSONCompatible]],
                               fragment_template: FileTemplate,
                               record_layer: RecordLayer,
                               region_layer: RegionLayer,
                               results_by_record_id: Dict[str, Dict[str, ModuleResults]],
                               html_sections: Dict[str, Dict[int, List[HTMLSections]]],
                               options_layer: OptionsLayer,
                               svg_tooltip: str) -> None:
    """ Writes one region's complete .js file: JSON data, module results,
        and the rendered HTML fragment.
    """
    js_path = os.path.join(regions_dir, f"{anchor}.js")
    with open(js_path, "w", encoding="utf-8") as fh:
        fh.write(f"all_regions[{json.dumps(anchor)}] = {json.dumps(region_json)};\n")
        if anchor in module_results:
            fh.write(f"resultsData[{json.dumps(anchor)}] = {json.dumps(module_results[anchor])};\n")
        content = fragment_template.render(
            record=record_layer, region=region_layer, options=options_layer,
            results_by_record_id=results_by_record_id,
            sections=html_sections,
            svg_tooltip=svg_tooltip,
            tta_name=tta.__name__, tfbs_name=tfbs.__name__,
        )
        fh.write(f"regionHTML[{json.dumps(anchor)}] = {json.dumps(content)};\n")


def _render_all_region_html(record_layers: List[RecordLayer],
                            per_region_json: Dict[str, Dict[str, Any]],
                            results_by_record_id: Dict[str, Dict[str, ModuleResults]],
                            html_sections: Dict[str, Dict[int, List[HTMLSections]]],
                            options_layer: OptionsLayer,
                            svg_tooltip: str) -> Dict[str, str]:
    """ Renders HTML fragments for all regions.

        Returns:
            a dict mapping anchor to rendered HTML string
    """
    fragment_template = FileTemplate(os.path.join(TEMPLATE_PATH, "region_fragment.html"))
    region_html: Dict[str, str] = {}

    for record_layer in record_layers:
        for region_layer in record_layer.regions:
            anchor = region_layer.anchor_id
            if anchor not in per_region_json:
                continue
            content = fragment_template.render(
                record=record_layer, region=region_layer, options=options_layer,
                results_by_record_id=results_by_record_id,
                sections=html_sections,
                svg_tooltip=svg_tooltip,
                tta_name=tta.__name__, tfbs_name=tfbs.__name__,
            )
            region_html[anchor] = content

    return region_html


def _write_all_region_data(record_layers: List[RecordLayer],
                            per_region_json: Dict[str, Dict[str, Any]],
                            module_results: dict[str, dict[str, json.JSONCompatible]],
                            results_by_record_id: Dict[str, Dict[str, ModuleResults]],
                            html_sections: Dict[str, Dict[int, List[HTMLSections]]],
                            options: ConfigType,
                            options_layer: OptionsLayer,
                            svg_tooltip: str) -> None:
    """ Renders all region HTML and writes a single compressed regions_data.js. """
    region_html = _render_all_region_html(
        record_layers, per_region_json, results_by_record_id,
        html_sections, options_layer, svg_tooltip,
    )
    _write_compressed_data_file(options.output_dir, per_region_json, module_results, region_html)


def write_regions_js(records: List[Dict[str, Any]], output_dir: str,
                     module_results: dict[str, dict[str, json.JSONCompatible]]) -> None:
    """ Writes a lightweight regions.js index plus a compressed regions_data.js.

        NOTE: This is the legacy entry point used by the classic (non-streaming)
        path when called without HTML fragment data. Per-region data written here
        will NOT include regionHTML; use generate_webpage() for the full flow.
    """
    lightweight_records, regions_index, per_region_json = _build_regions_index(records)
    _write_regions_index(lightweight_records, regions_index, output_dir)

    _write_compressed_data_file(output_dir, per_region_json, module_results, {})


def generate_html_sections(records: List[RecordLayer], results: Dict[str, Dict[str, ModuleResults]],
                           options: ConfigType) -> Dict[str, Dict[int, List[HTMLSections]]]:
    """ Generates a mapping of record->region->HTMLSections for each record, region and module

        Arguments:
            records: a list of RecordLayers to pass through to the modules
            results: a dictionary mapping record name to
                        a dictionary mapping each module name to its results object
            options: the current antiSMASH config

        Returns:
            a dictionary mapping record id to
                a dictionary mapping region number to
                    a list of HTMLSections, one for each module
    """
    details = {}
    for record in records:
        record_details = {}
        record_result = results[record.id]
        for region in record.regions:
            sections = []
            for handler in region.handlers:
                if handler.will_handle(region.products, region.product_categories):
                    handler_results = record_result.get(handler.__name__)
                    if handler_results is None:
                        continue
                    sections.append(handler.generate_html(region, handler_results, record, options))
            for aggregator in VISUALISERS:
                if not hasattr(aggregator, "generate_html"):
                    continue
                if aggregator.has_enough_results(record.seq_record, region.region_feature, record_result):
                    section = aggregator.generate_html(region, record_result, record, options)
                    # as a special case, the first section of a region should always be the gene table
                    if aggregator is gene_table:
                        sections.insert(0, section)
                    else:
                        sections.append(section)
            record_details[region.get_region_number()] = sections
        details[record.id] = record_details
    return details


def find_local_antismash_js_path(options: ConfigType) -> Optional[str]:
    """ Finds the a path to a local copy of antismash.js, if possible,
        otherwise returns None.
    """
    # is a copy in the js directory?
    js_path = path.locate_file(path.get_full_path(__file__, "js", "antismash.js"), silent=True)
    if js_path:
        return js_path

    # is it in the databases?
    version = get_antismash_js_version()
    js_path = path.locate_file(os.path.join(options.database_dir, "as-js", version, "antismash.js"), silent=True)
    if js_path:
        return js_path

    # then it doesn't exist
    return None


def build_antismash_js_url(options: ConfigType) -> str:
    """ Build the URL to the javascript that will be embedded in the HTML.
        If a local version is available, it will be copied into the output directory,
        otherwise a full remote URL will be used.

        Arguments:
            options: the antiSMASH config

        Returns:
            a string of the URL, whether relative or absolute
    """
    if find_local_antismash_js_path(options):
        return "js/antismash.js"  # generic local path after copy
    return get_antismash_js_url()


def _build_record_layers(records: List[Record], results: List[Dict[str, ModuleResults]],
                         options_layer: OptionsLayer,
                         ) -> Tuple[List[RecordLayer], List[RecordLayer], Dict[str, Dict[str, ModuleResults]]]:
    """ Constructs RecordLayer objects, split by whether they have regions.

        Returns:
            a tuple of
                record layers with regions,
                record layers without regions,
                results keyed by record id
    """
    record_layers_with_regions = []
    record_layers_without_regions = []
    results_by_record_id: Dict[str, Dict[str, ModuleResults]] = {}
    for record, record_results in zip(records, results):
        if record.get_regions():
            record_layers_with_regions.append(RecordLayer(record, None, options_layer))
        else:
            record_layers_without_regions.append(RecordLayer(record, None, options_layer))
        results_by_record_id[record.id] = record_results
    return record_layers_with_regions, record_layers_without_regions, results_by_record_id


_CACHED_FRAGMENT_TEMPLATE: Optional[FileTemplate] = None
_CACHED_SVG_TOOLTIP: Optional[str] = None
_REGIONS_DIR_CREATED: Optional[str] = None


def _get_fragment_template() -> FileTemplate:
    """Return a cached FileTemplate for region_fragment.html."""
    global _CACHED_FRAGMENT_TEMPLATE
    if _CACHED_FRAGMENT_TEMPLATE is None:
        _CACHED_FRAGMENT_TEMPLATE = FileTemplate(os.path.join(TEMPLATE_PATH, "region_fragment.html"))
    return _CACHED_FRAGMENT_TEMPLATE


def _get_svg_tooltip() -> str:
    """Return a cached SVG tooltip string."""
    global _CACHED_SVG_TOOLTIP
    if _CACHED_SVG_TOOLTIP is None:
        _CACHED_SVG_TOOLTIP = _build_svg_tooltip()
    return _CACHED_SVG_TOOLTIP


def _ensure_regions_dir(output_dir: str) -> str:
    """Create the regions/ subdirectory once, return its path."""
    global _REGIONS_DIR_CREATED
    regions_dir = os.path.join(output_dir, "regions")
    if _REGIONS_DIR_CREATED != regions_dir:
        os.makedirs(regions_dir, exist_ok=True)
        _REGIONS_DIR_CREATED = regions_dir
    return regions_dir


def _build_svg_tooltip() -> str:
    """ Builds the standard SVG tooltip text used across region rendering. """
    svg_tooltip = ("Shows the layout of the region, marking coding sequences and areas of interest. "
                   "Clicking a gene will select it and show any relevant details. "
                   "Clicking an area feature (e.g. a candidate cluster) will select all coding "
                   "sequences within that area. Double clicking an area feature will zoom to that area. "
                   "Multiple genes and area features can be selected by clicking them while holding the Ctrl key."
                   )
    doc_target = "understanding_output/#the-antismash-5-region-concept"
    svg_tooltip += f"<br>More detailed help is available {docs_link('here', doc_target)}."
    return svg_tooltip


def _build_json_data_for_record(record: Record, record_results: Dict[str, ModuleResults],
                                 options: ConfigType, all_modules: List[AntismashModule],
                                 ) -> Tuple[Dict[str, Any], dict[str, dict[str, json.JSONCompatible]]]:
    """ Single-record version of build_json_data.

        Returns:
            a tuple of
                the JSON-friendly record dict,
                module JS results keyed by region anchor
    """
    json_record = js.convert_record(record, options, record_results)
    json_record['seq_id'] = "".join(
        char for char in json_record['seq_id'] if char in string.printable
    )

    js_results: dict[str, dict[str, json.JSONCompatible]] = {}
    for region, json_region in zip(record.get_regions(), json_record['regions']):
        handlers = find_plugins_for_cluster(all_modules, json_region)
        region_results = {}
        for handler in handlers:
            if handler.__name__ not in record_results:
                continue
            if hasattr(handler, "generate_javascript_data"):
                data = handler.generate_javascript_data(record, region, record_results[handler.__name__])
                region_results[handler.__name__] = data

        for aggregator in VISUALISERS:
            if not hasattr(aggregator, "generate_javascript_data"):
                continue
            if aggregator.has_enough_results(record, region, record_results):
                data = aggregator.generate_javascript_data(record, region, record_results)
                region_results[aggregator.__name__] = data

        if region_results:
            js_results[RegionLayer.build_anchor_id(region)] = region_results

    return json_record, js_results


def _generate_html_sections_for_record(record_layer: RecordLayer,
                                        record_results: Dict[str, ModuleResults],
                                        options: ConfigType,
                                        ) -> Dict[int, List[HTMLSections]]:
    """ Single-record version of generate_html_sections.

        Returns:
            a dict mapping region number to a list of HTMLSections
    """
    record_details: Dict[int, List[HTMLSections]] = {}
    for region in record_layer.regions:
        sections = []
        for handler in region.handlers:
            if handler.will_handle(region.products, region.product_categories):
                handler_results = record_results.get(handler.__name__)
                if handler_results is None:
                    continue
                sections.append(handler.generate_html(region, handler_results, record_layer, options))
        for aggregator in VISUALISERS:
            if not hasattr(aggregator, "generate_html"):
                continue
            if aggregator.has_enough_results(record_layer.seq_record, region.region_feature, record_results):
                section = aggregator.generate_html(region, record_results, record_layer, options)
                if aggregator is gene_table:
                    sections.insert(0, section)
                else:
                    sections.append(section)
        record_details[region.get_region_number()] = sections
    return record_details


def generate_region_files_for_record(record: Record, record_results: Dict[str, ModuleResults],
                                      options: ConfigType, all_modules: List[AntismashModule],
                                      options_layer: OptionsLayer,
                                      data_writer: Optional[StreamingRegionDataWriter] = None,
                                      ) -> Tuple[Dict[str, Any], Optional[RecordLayer]]:
    """ Streaming entry point: process a single record and write its per-region
        data immediately. Returns lightweight record metadata for the index
        and the RecordLayer (or None if no regions).

        When data_writer is provided, region data is written to the single
        compressed regions_data.js file. Otherwise falls back to writing
        individual per-region .js files.

        This function is called incrementally as each record completes analysis,
        keeping memory bounded by writing and discarding per-region JSON data
        before the next record is processed.

        Returns:
            a tuple of
                lightweight JSON record dict (for the regions.js index),
                RecordLayer if the record has regions (else None)
    """
    # Build JSON data for this record
    json_record, js_results = _build_json_data_for_record(record, record_results, options, all_modules)

    # Build lightweight record for the index
    light_record = {k: v for k, v in json_record.items() if k != "regions"}
    light_record["regions"] = [_make_lightweight_region(r) for r in json_record["regions"]]

    if not record.get_regions():
        return light_record, None

    # Build record layer and HTML sections
    record_layer = RecordLayer(record, None, options_layer)
    results_by_record_id = {record.id: record_results}
    html_sections = {record.id: _generate_html_sections_for_record(record_layer, record_results, options)}

    svg_tooltip = _get_svg_tooltip()
    fragment_template = _get_fragment_template()

    # Write region data for this record
    for region_layer, json_region in zip(record_layer.regions, json_record["regions"]):
        anchor = region_layer.anchor_id
        content = fragment_template.render(
            record=record_layer, region=region_layer, options=options_layer,
            results_by_record_id=results_by_record_id,
            sections=html_sections,
            svg_tooltip=svg_tooltip,
            tta_name=tta.__name__, tfbs_name=tfbs.__name__,
        )
        if data_writer is not None:
            data_writer.write_region(
                anchor, json_region,
                js_results.get(anchor, {}),
                content,
            )
        else:
            # Fallback: write individual files
            regions_dir = _ensure_regions_dir(options.output_dir)
            _write_single_region_file(
                regions_dir, anchor, json_region, js_results,
                fragment_template, record_layer, region_layer,
                results_by_record_id, html_sections,
                options_layer, svg_tooltip,
            )

    return light_record, record_layer


def finalize_html_output(lightweight_records: List[Dict[str, Any]],
                          record_layers_with_regions: List[RecordLayer],
                          record_layers_without_regions: List[RecordLayer],
                          results_by_record_id: Dict[str, Dict[str, ModuleResults]],
                          options: ConfigType, all_modules: List[AntismashModule],
                          taxonomy_mapping: Optional[Dict[str, str]] = None,
                          skipped_record_count: int = 0,
                          ) -> None:
    """ Streaming finalization: write the regions.js index, regions.html overview,
        and index.html dashboard after all records have been processed.

        Must be called after all generate_region_files_for_record() calls.
    """
    options_layer = OptionsLayer(options, all_modules)

    # Build the regions index from lightweight records
    regions_index: Dict[str, Any] = {"order": []}
    for light_record in lightweight_records:
        for light_region in light_record["regions"]:
            anchor = light_region["anchor"]
            regions_index[anchor] = light_region
            regions_index["order"].append(anchor)

    _write_regions_index(lightweight_records, regions_index, options.output_dir)

    # In streaming mode, _generate_html_sections_for_record() already ran
    # for every record, populating most_related_area as a side effect.
    # The html_sections dict itself is unused by overview.html with lazy_loading=True.
    # We only need a valid (possibly empty) dict for the template context.
    html_sections: Dict[str, Dict[int, List[HTMLSections]]] = {}

    regions_written = sum(len(rl.regions) for rl in record_layers_with_regions)
    job_id = os.path.basename(options.output_dir)
    page_title = options.output_basename
    if options.html_title:
        page_title = options.html_title

    svg_tooltip = _build_svg_tooltip()
    as_js_url = build_antismash_js_url(options)

    template = FileTemplate(os.path.join(TEMPLATE_PATH, "overview.html"))
    regions_content = template.render(
        records=record_layers_with_regions, options=options_layer,
        version=options.version,
        regions_written=regions_written, sections=html_sections,
        results_by_record_id=results_by_record_id,
        config=options, job_id=job_id, page_title=page_title,
        records_without_regions=record_layers_without_regions,
        skipped_record_count=skipped_record_count,
        svg_tooltip=svg_tooltip, get_region_css=js.get_region_css,
        as_js_url=as_js_url, tta_name=tta.__name__, tfbs_name=tfbs.__name__,
        lazy_loading=True,
    )
    regions_content = _strip_leading_whitespace(regions_content)
    with open(os.path.join(options.output_dir, "regions.html"), "w", encoding="utf-8") as fh:
        fh.write(regions_content)

    # Generate the dashboard
    dashboard_content = generate_dashboard(
        record_layers_with_regions, record_layers_without_regions,
        options_layer, page_title, taxonomy_mapping=taxonomy_mapping,
        skipped_record_count=skipped_record_count,
    )
    dashboard_content = _strip_leading_whitespace(dashboard_content)
    with open(os.path.join(options.output_dir, "index.html"), "w", encoding="utf-8") as fh:
        fh.write(dashboard_content)


def generate_webpage(records: List[Record], results: List[Dict[str, ModuleResults]],
                     options: ConfigType, all_modules: List[AntismashModule],
                     skipped_record_count: int = 0) -> str:
    """ Generates the HTML itself """

    # Phase 1: Build all data structures (sequential)
    json_records, js_results = build_json_data(records, results, options, all_modules)

    lightweight_records, regions_index, per_region_json = _build_regions_index(json_records)
    _write_regions_index(lightweight_records, regions_index, options.output_dir)

    template = FileTemplate(os.path.join(TEMPLATE_PATH, "overview.html"))

    options_layer = OptionsLayer(options, all_modules)
    record_layers_with_regions, record_layers_without_regions, results_by_record_id = \
        _build_record_layers(records, results, options_layer)

    regions_written = sum(len(record.get_regions()) for record in records)
    job_id = os.path.basename(options.output_dir)
    page_title = options.output_basename
    if options.html_title:
        page_title = options.html_title

    html_sections = generate_html_sections(record_layers_with_regions, results_by_record_id, options)

    svg_tooltip = _build_svg_tooltip()

    as_js_url = build_antismash_js_url(options)

    content = template.render(records=record_layers_with_regions, options=options_layer,
                              version=options.version,
                              regions_written=regions_written, sections=html_sections,
                              results_by_record_id=results_by_record_id,
                              config=options, job_id=job_id, page_title=page_title,
                              records_without_regions=record_layers_without_regions,
                              skipped_record_count=skipped_record_count,
                              svg_tooltip=svg_tooltip, get_region_css=js.get_region_css,
                              as_js_url=as_js_url, tta_name=tta.__name__, tfbs_name=tfbs.__name__,
                              lazy_loading=True,
                              )

    # Phase 2: Render HTML and write compressed data file
    _write_all_region_data(record_layers_with_regions, per_region_json, js_results,
                           results_by_record_id, html_sections, options, options_layer, svg_tooltip)

    return content, record_layers_with_regions, record_layers_without_regions, options_layer, page_title


def generate_dashboard(record_layers_with_regions: List[RecordLayer],
                       record_layers_without_regions: List[RecordLayer],
                       options_layer: OptionsLayer, page_title: str,
                       taxonomy_mapping: Optional[Dict[str, str]] = None,
                       skipped_record_count: int = 0) -> str:
    """ Generates the dashboard HTML page.

        Must be called AFTER generate_webpage / generate_html_sections,
        because those populate most_related_area on each RegionLayer.
    """
    from markupsafe import Markup
    from antismash.outputs.html.dashboard import build_dashboard_data

    template = FileTemplate(os.path.join(TEMPLATE_PATH, "dashboard.html"))
    data = build_dashboard_data(record_layers_with_regions, record_layers_without_regions,
                                taxonomy_mapping=taxonomy_mapping or {},
                                skipped_record_count=skipped_record_count)
    has_taxonomy = bool(taxonomy_mapping)

    content = template.render(
        options=options_layer,
        page_title=page_title,
        data=data,
        dashboard_json=Markup(json.dumps(data)),
        has_taxonomy=has_taxonomy,
    )
    return content


def find_plugins_for_cluster(plugins: List[AntismashModule],
                             cluster: Dict[str, Any]) -> List[AntismashModule]:
    "Find a specific plugin responsible for a given gene cluster type"
    products = cluster['products']
    categories = set(cluster['product_categories'])
    handlers = []
    for plugin in plugins:
        if not hasattr(plugin, 'will_handle'):
            continue
        if plugin.will_handle(products, categories):
            handlers.append(plugin)
    return handlers
