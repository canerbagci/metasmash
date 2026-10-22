# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import base64
import gzip
import json
import os
import unittest
from unittest.mock import Mock, patch

from helperlibs.wrappers.io import TemporaryDirectory

from antismash.main import get_all_modules
from antismash.common import html_renderer, path
from antismash.common.layers import OptionsLayer
from antismash.common.secmet.test.helpers import DummyCDS, DummyProtocluster, DummySubRegion
from antismash.common.test.helpers import DummyRecord
from antismash.config import build_config, destroy_config, update_config
from antismash.outputs.html import generator
from antismash.outputs.html.dashboard import build_dashboard_data_from_summaries
from antismash.outputs.html.generator import (
    generate_webpage,
    build_antismash_js_url,
    finalize_streaming_html_output,
)
from antismash.outputs.html.streaming_summary import StreamingRecordSummary, StreamingRegionSummary, StreamingRelatedAreaSummary


class TestOutput(unittest.TestCase):
    def setUp(self):
        self.options = build_config([], modules=get_all_modules())

    def tearDown(self):
        destroy_config()

    def test_subregions_only(self):
        # construct a region with only subregions
        record = DummyRecord(seq="A" * 10)
        record.add_subregion(DummySubRegion())
        record.create_regions()

        # add information normally added by record processing
        record.record_index = 1

        # make sure HTML is generated without errors
        with TemporaryDirectory() as temp_dir:
            update_config({"output_dir": temp_dir})
            assert generate_webpage([record], [{}], self.options, get_all_modules())

    def test_generate_webpage_writes_region_data_bundles(self):
        record = DummyRecord(seq="A" * 200)
        record.record_index = 1
        record.add_cds_feature(DummyCDS(start=20, end=80, locus_tag="cds1", translation="A" * 20))
        record.add_protocluster(
            DummyProtocluster(core_start=30, core_end=50,
                              product="test-product", product_category="TEST-CATEGORY")
        )
        record.create_candidate_clusters()
        record.create_regions()

        with TemporaryDirectory() as temp_dir:
            update_config({"output_dir": temp_dir})
            assert generate_webpage([record], [{}], self.options, get_all_modules())
            assert _read_bundle_index(temp_dir) == {"r1c1": 1}
            with open(os.path.join(temp_dir, "regions_data", "0001.js"), encoding="utf-8") as handle:
                assert '_regionsData["r1c1"] = "' in handle.read()
            assert not os.path.exists(os.path.join(temp_dir, "regions"))
            assert not os.path.exists(os.path.join(temp_dir, "regions_data.js"))
            with open(os.path.join(temp_dir, "region.html"), encoding="utf-8") as handle:
                page = handle.read()
            assert "js/region_page.js" in page
            assert "regions_data/index.js" in page
            assert "regions.js" not in page


class TestJavascriptPresence(unittest.TestCase):
    def setUp(self):
        self.version = html_renderer.get_antismash_js_version()

    def build(self, return_values):
        options = Mock()
        options.database_dir = "dummy_data"
        with patch.object(path, "locate_file", side_effect=return_values) as patched:
            url = build_antismash_js_url(options)
            call_args = patched.call_args_list

        # ensure it wasn't called more times than expected
        assert len(call_args) == len(return_values)

        expected_args = [  # enough to cover all cases, only those used are checked
            path.get_full_path(generator.__file__, "js", "antismash.js"),
            os.path.join(options.database_dir, "as-js", self.version, "antismash.js"),
        ]
        for called, expected in zip(call_args, expected_args):
            assert os.path.realpath(called[0][0]) == os.path.realpath(expected)
        return url

    def test_fallback_to_web(self):
        url = self.build([None, None])
        assert url == html_renderer.get_antismash_js_url()

    def test_fallback_to_data(self):
        url = self.build([None, "data_dir"])
        assert url == "js/antismash.js"

    def test_fallback_to_in_tree(self):
        url = self.build(["local"])
        assert url == "js/antismash.js"


class TestStreamingHtmlSummaries(unittest.TestCase):
    class FakeArea:
        identifier = "BGC000001"
        description = "Known cluster"
        product = "test-product"
        similarity_percentage = 87
        url = "https://example.invalid/cluster"

    def setUp(self):
        self.options = build_config([], modules=get_all_modules())

    def tearDown(self):
        destroy_config()

    @staticmethod
    def _make_record():
        record = DummyRecord(seq="A" * 200, record_id="record-1")
        record.record_index = 7
        record.original_id = "original-record-1"
        record.annotations["source"] = "test-source"
        # CDS end aligns with the protocluster neighbourhood (core 30-50 + range 10 => region 20-60),
        # so secmet populates region.cds_children — required for gene_count to be non-zero.
        record.add_cds_feature(DummyCDS(start=20, end=60, locus_tag="cds1", translation="A" * 20))
        record.add_protocluster(
            DummyProtocluster(core_start=30, core_end=50,
                              product="test-product", product_category="TEST-CATEGORY")
        )
        record.create_candidate_clusters()
        record.create_regions()
        return record

    def test_generate_region_files_returns_summary(self):
        record = self._make_record()
        options_layer = OptionsLayer(self.options, get_all_modules())
        data_writer = Mock()

        def fake_sections(record_layer, _record_results, _options):
            record_layer.regions[0].most_related_area = self.FakeArea()
            return {1: []}

        with TemporaryDirectory() as temp_dir:
            update_config({"output_dir": temp_dir})
            with patch.object(generator, "_generate_html_sections_for_record", side_effect=fake_sections), \
                 patch.object(generator, "_get_fragment_template",
                              return_value=Mock(render=Mock(return_value="<div>fragment</div>"))):
                light_record, summary = generator.generate_region_files_for_record(
                    record, {}, self.options, get_all_modules(),
                    options_layer, data_writer=data_writer,
                )

        assert summary is not None
        assert summary.id == "record-1"
        assert summary.record_index == 7
        assert summary.original_id == "original-record-1"
        assert summary.source == "test-source"
        assert summary.get_name() == "record-1"
        assert len(summary.regions) == 1
        region = summary.regions[0]
        assert region.anchor_id == light_record["regions"][0]["anchor"]
        assert region.css_class == "TEST-CATEGORY test-product"
        assert region.gene_count == 1
        assert region.product_to_category["test-product"] == "TEST-CATEGORY"
        assert region.most_related_area.description == "Known cluster"
        data_writer.write_region_encoded.assert_called_once()

    def test_generate_region_files_writes_to_data_bundles(self):
        record = self._make_record()
        options_layer = OptionsLayer(self.options, get_all_modules())

        with TemporaryDirectory() as temp_dir:
            update_config({"output_dir": temp_dir})
            data_writer = generator.StreamingRegionDataWriter(temp_dir)
            with patch.object(generator, "_generate_html_sections_for_record", return_value={1: []}), \
                 patch.object(generator, "_get_fragment_template",
                              return_value=Mock(render=Mock(return_value="<div>fragment</div>"))):
                light_record, summary = generator.generate_region_files_for_record(
                    record, {}, self.options, get_all_modules(),
                    options_layer, data_writer=data_writer,
                )
            data_writer.finalize()

            assert light_record["regions"][0]["anchor"] == "r7c1"
            assert summary is not None
            # no per-region files
            assert not os.path.exists(os.path.join(temp_dir, "regions"))
            payload = _read_bundled_region(temp_dir, "r7c1")
            assert payload["region"]["anchor"] == "r7c1"
            assert payload["html"] == "<div>fragment</div>"


    def test_build_dashboard_data_from_summaries(self):
        summary = StreamingRecordSummary(
            id="record-1",
            record_index=3,
            original_id="",
            source="",
            combined_source_count=1,
            regions=[
                StreamingRegionSummary(
                    anchor_id="r3c1",
                    region_number=1,
                    start=10,
                    end=120,
                    length=110,
                    products=["test-product"],
                    product_categories=["TEST-CATEGORY"],
                    css_class="TEST-CATEGORY test-product",
                    contig_edge=False,
                    gene_count=4,
                    product_to_category={"test-product": "TEST-CATEGORY"},
                    most_related_area=StreamingRelatedAreaSummary(
                        identifier="BGC000001",
                        description="Known cluster",
                        product="test-product",
                        similarity_percentage=87,
                        url="https://example.invalid/cluster",
                    ),
                ),
            ],
        )
        data = build_dashboard_data_from_summaries(
            [summary], taxonomy_mapping={"record-1": "Root; Test"}, skipped_record_count=2,
        )
        assert data["summary"]["total_records"] == 3
        assert data["summary"]["records_with_bgcs"] == 1
        assert data["summary"]["total_bgcs"] == 1
        assert data["bgc_table"][0]["gene_count"] == 4
        assert data["bgc_table"][0]["most_similar_known"] == "Known cluster"
        assert data["bgc_table"][0]["taxonomy"] == "Root; Test"

    def test_finalize_streaming_html_output_defers_hidden_overview(self):
        # 51 records: two overview pages of 50, plus the compact (multi-record) table
        summaries, lightweight_records = [], []
        for i in range(1, 52):
            anchor = f"r{i}c1"
            summaries.append(StreamingRecordSummary(
                id=f"record-{i}", record_index=i, original_id="", source="", combined_source_count=1,
                regions=[StreamingRegionSummary(
                    anchor_id=anchor, region_number=1, start=5, end=45, length=40,
                    products=["test-product"], product_categories=["TEST-CATEGORY"],
                    css_class="TEST-CATEGORY test-product", contig_edge=False, gene_count=2,
                    product_to_category={"test-product": "TEST-CATEGORY"},
                )],
            ))
            lightweight_records.append({"seq_id": f"record-{i}", "regions": [{
                "anchor": anchor, "start": 5, "end": 45, "idx": 1, "type": "test-product",
                "products": ["test-product"], "product_categories": ["TEST-CATEGORY"],
                "cssClass": "TEST-CATEGORY test-product"}]})

        with TemporaryDirectory() as temp_dir:
            update_config({"output_dir": temp_dir})
            finalize_streaming_html_output(lightweight_records, summaries, self.options, get_all_modules())
            with open(os.path.join(temp_dir, "regions.html"), encoding="utf-8") as handle:
                content = handle.read()

        deferred_start = '<script type="text/x-template" class="overview-deferred">'
        blocks = [block.split("</script>", 1)[0] for block in content.split(deferred_start)[1:]]
        # the second overview page and the compact table are deferred, the first page is not
        assert len(blocks) == 2
        assert 'data-anchor="#r51c1"' in blocks[0] and 'data-anchor="#r50c1"' not in blocks[0]
        assert all(f'data-anchor="#r{i}c1"' in blocks[1] for i in range(1, 52))
        first_page = content.split('data-page="0"', 1)[1].split('data-page="1"', 1)[0]
        assert 'data-anchor="#r1c1"' in first_page and deferred_start not in first_page

    def test_finalize_streaming_html_output_sorts_record_summaries(self):
        record_b = StreamingRecordSummary(
            id="record-b",
            record_index=2,
            original_id="",
            source="",
            combined_source_count=1,
            regions=[
                StreamingRegionSummary(
                    anchor_id="r2c1",
                    region_number=1,
                    start=10,
                    end=50,
                    length=40,
                    products=["test-product"],
                    product_categories=["TEST-CATEGORY"],
                    css_class="TEST-CATEGORY test-product",
                    contig_edge=False,
                    gene_count=2,
                    product_to_category={"test-product": "TEST-CATEGORY"},
                ),
            ],
        )
        record_a = StreamingRecordSummary(
            id="record-a",
            record_index=1,
            original_id="",
            source="",
            combined_source_count=1,
            regions=[
                StreamingRegionSummary(
                    anchor_id="r1c1",
                    region_number=1,
                    start=5,
                    end=45,
                    length=40,
                    products=["test-product"],
                    product_categories=["TEST-CATEGORY"],
                    css_class="TEST-CATEGORY test-product",
                    contig_edge=False,
                    gene_count=2,
                    product_to_category={"test-product": "TEST-CATEGORY"},
                ),
            ],
        )
        lightweight_records = [
            {"seq_id": "record-b", "regions": [{"anchor": "r2c1", "start": 10, "end": 50,
                                                 "idx": 1, "type": "test-product",
                                                 "products": ["test-product"],
                                                 "product_categories": ["TEST-CATEGORY"],
                                                 "cssClass": "TEST-CATEGORY test-product"}]},
            {"seq_id": "record-a", "regions": [{"anchor": "r1c1", "start": 5, "end": 45,
                                                 "idx": 1, "type": "test-product",
                                                 "products": ["test-product"],
                                                 "product_categories": ["TEST-CATEGORY"],
                                                 "cssClass": "TEST-CATEGORY test-product"}]},
        ]

        with TemporaryDirectory() as temp_dir:
            update_config({"output_dir": temp_dir})
            finalize_streaming_html_output(
                lightweight_records, [record_b, record_a], self.options, get_all_modules(),
            )
            with open(os.path.join(temp_dir, "regions.html"), encoding="utf-8") as handle:
                content = handle.read()
            with open(os.path.join(temp_dir, "regions.js"), encoding="utf-8") as handle:
                regions_js = handle.read()

        assert content.index("record-a") < content.index("record-b")
        # the region order (used by the prev/next arrows) follows the same record order
        order = json.loads(regions_js.split("var all_regions = ", 1)[1].split(";\n", 1)[0])["order"]
        assert order == ["r1c1", "r2c1"]
        records = json.loads(regions_js.split("var recordData = ", 1)[1].split(";\n", 1)[0])
        assert [record["seq_id"] for record in records] == ["record-a", "record-b"]


class TestRegionDataBundles(unittest.TestCase):
    def test_bundles_roll_over_at_max_bytes(self):
        with TemporaryDirectory() as temp_dir:
            writer = generator.StreamingRegionDataWriter(temp_dir, max_bytes=100)
            for i in range(5):
                writer.write_region(f"r{i}c1", {"anchor": f"r{i}c1"}, {}, "x" * 200)
            writer.finalize()
            bundle_dir = os.path.join(temp_dir, generator.REGIONS_DATA_DIR)
            # each region is larger than max_bytes, so each starts a new bundle
            assert sorted(os.listdir(bundle_dir)) == ["0001.js", "0002.js", "0003.js",
                                                      "0004.js", "0005.js", "index.js"]
            for i in range(5):
                assert _read_bundled_region(temp_dir, f"r{i}c1")["region"] == {"anchor": f"r{i}c1"}

    def test_bundles_share_files_below_max_bytes(self):
        with TemporaryDirectory() as temp_dir:
            writer = generator.StreamingRegionDataWriter(temp_dir)
            for i in range(5):
                writer.write_region(f"r{i}c1", {"anchor": f"r{i}c1"}, {"mod": i}, "html")
            writer.finalize()
            index = _read_bundle_index(temp_dir)
            assert index == {f"r{i}c1": 1 for i in range(5)}
            assert _read_bundled_region(temp_dir, "r3c1")["results"] == {"mod": 3}

    def test_empty_writer_still_writes_index(self):
        with TemporaryDirectory() as temp_dir:
            generator.StreamingRegionDataWriter(temp_dir).finalize()
            assert _read_bundle_index(temp_dir) == {}
            assert os.listdir(os.path.join(temp_dir, generator.REGIONS_DATA_DIR)) == ["index.js"]


def _read_bundle_index(output_dir):
    with open(os.path.join(output_dir, generator.REGIONS_DATA_DIR, "index.js"), encoding="utf-8") as handle:
        content = handle.read()
    prefix = "var _regionsDataIndex = "
    assert content.startswith(prefix)
    return json.loads(content[len(prefix):].rstrip().rstrip(";"))


def _read_bundled_region(output_dir, anchor):
    """ Decodes one region's payload from its bundle, as lazy_loader.js does """
    bundle = _read_bundle_index(output_dir)[anchor]
    path = os.path.join(output_dir, generator.REGIONS_DATA_DIR, f"{bundle:04d}.js")
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(f'_regionsData["{anchor}"] = "'):
                encoded = line.split('"')[3]
                return json.loads(gzip.decompress(base64.b64decode(encoded)))
    raise AssertionError(f"{anchor} not in bundle {bundle}")
