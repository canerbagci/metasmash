# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

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
        record.add_cds_feature(DummyCDS(start=20, end=80, locus_tag="cds1", translation="A" * 20))
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
        data_writer.write_region.assert_called_once()

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

        assert content.index("record-a") < content.index("record-b")
