# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

"""Tests for the streaming pipeline functions:
    - should_use_streaming (main.py)
    - _safe_process_record_full (main.py)
"""

import os
import tempfile
import unittest
from unittest import mock

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from antismash import config


class TestShouldUseStreaming(unittest.TestCase):
    """Tests for main.should_use_streaming()"""

    def setUp(self):
        config.build_config([], isolated=True)

    def tearDown(self):
        config.destroy_config()

    def _make_options(self, **overrides):
        """Build a simple namespace with the attributes should_use_streaming reads."""
        defaults = {
            "streaming": "auto",
            "reuse_results": None,
            "genefinding_gff3": None,
            "taxon": "bacteria",
            "minlength": -1,
            "start": -1,
            "end": -1,
            "abort_on_invalid_records": True,
        }
        defaults.update(overrides)
        return type("Options", (), defaults)()

    def test_streaming_off(self):
        from antismash.main import should_use_streaming
        opts = self._make_options(streaming="off")
        result, records = should_use_streaming(opts, "some_file.gbk")
        assert result is False
        assert records is None

    def test_streaming_on(self):
        from antismash.main import should_use_streaming
        opts = self._make_options(streaming="on")
        result, records = should_use_streaming(opts, "some_file.gbk")
        assert result is True
        assert records is None

    def test_streaming_on_with_reuse_results(self):
        from antismash.main import should_use_streaming
        opts = self._make_options(streaming="on", reuse_results="/some/path")
        result, records = should_use_streaming(opts, "some_file.gbk")
        assert result is False
        assert records is None

    def test_streaming_on_with_gff3(self):
        from antismash.main import should_use_streaming
        opts = self._make_options(streaming="on", genefinding_gff3="/some/gff3")
        result, records = should_use_streaming(opts, "some_file.gbk")
        assert result is False
        assert records is None

    def test_auto_no_sequence_file(self):
        from antismash.main import should_use_streaming
        opts = self._make_options(streaming="auto")
        result, records = should_use_streaming(opts, None)
        assert result is False
        assert records is None
        result, records = should_use_streaming(opts, "")
        assert result is False
        assert records is None

    def test_auto_with_reuse_results(self):
        from antismash.main import should_use_streaming
        opts = self._make_options(streaming="auto", reuse_results="/some/path")
        result, records = should_use_streaming(opts, "some_file.gbk")
        assert result is False
        assert records is None

    def test_auto_with_gff3(self):
        from antismash.main import should_use_streaming
        opts = self._make_options(streaming="auto", genefinding_gff3="/some/gff3")
        result, records = should_use_streaming(opts, "some_file.gbk")
        assert result is False
        assert records is None

    def test_auto_many_records(self):
        from antismash.main import should_use_streaming, _STREAMING_RECORD_THRESHOLD
        opts = self._make_options(streaming="auto")
        count = _STREAMING_RECORD_THRESHOLD + 1
        fake_metadata = [(f"rec{i}", 400) for i in range(count)]
        with mock.patch("antismash.common.record_processing.scan_record_metadata",
                        return_value=fake_metadata):
            result, metadata = should_use_streaming(opts, "some_file.gbk")
        assert result is True
        assert isinstance(metadata, list)
        assert len(metadata) == count

    def test_auto_few_records(self):
        from antismash.main import should_use_streaming
        opts = self._make_options(streaming="auto")
        fake_metadata = [(f"rec{i}", 400) for i in range(3)]
        with mock.patch("antismash.common.record_processing.scan_record_metadata",
                        return_value=fake_metadata):
            result, metadata = should_use_streaming(opts, "some_file.gbk")
        assert result is False
        assert metadata is None

    def test_auto_scan_failure_falls_back(self):
        from antismash.main import should_use_streaming
        opts = self._make_options(streaming="auto")
        with mock.patch("antismash.common.record_processing.scan_record_metadata",
                        side_effect=Exception("scan error")):
            result, metadata = should_use_streaming(opts, "some_file.gbk")
        assert result is False
        assert metadata is None


class TestSafeProcessRecordFull(unittest.TestCase):
    """Tests for main._safe_process_record_full()"""

    def test_normal_return_passes_through(self):
        from antismash.main import _safe_process_record_full
        expected = ("record", {"module": "result"}, {"timing": 1.0})
        record_tuple = (SeqRecord(Seq("ACGT"), id="test"), 1)
        opts = mock.MagicMock()
        opts.abort_on_invalid_records = True
        with mock.patch("antismash.main.process_record_full", return_value=expected):
            result = _safe_process_record_full(record_tuple, opts)
        assert result == expected

    def test_exception_propagates_when_abort(self):
        from antismash.main import _safe_process_record_full
        record_tuple = (SeqRecord(Seq("ACGT"), id="test"), 1)
        opts = mock.MagicMock()
        opts.abort_on_invalid_records = True
        with mock.patch("antismash.main.process_record_full", side_effect=RuntimeError("boom")):
            with self.assertRaises(RuntimeError):
                _safe_process_record_full(record_tuple, opts)

    def test_exception_returns_sentinel_when_no_abort(self):
        from antismash.main import _safe_process_record_full, _RECORD_FAILED
        record_tuple = (SeqRecord(Seq("ACGT"), id="test"), 1)
        opts = mock.MagicMock()
        opts.abort_on_invalid_records = False
        with mock.patch("antismash.main.process_record_full", side_effect=RuntimeError("boom")):
            result = _safe_process_record_full(record_tuple, opts)
        assert result[0] == _RECORD_FAILED
        assert result[1] == "test"
        assert "Traceback" in result[2]
        assert "boom" in result[2]


class TestStreamingWindowHelpers(unittest.TestCase):
    def test_compute_phase1_batch_size(self):
        from antismash.main import _STREAMING_PHASE2_WINDOW_MIN, _compute_streaming_phase1_batch_size
        opts = type("Options", (), {
            "workers": 1,
            "streaming_phase1_batch_size": 0,
            "streaming_phase2_window_size": 0,
        })()
        assert _compute_streaming_phase1_batch_size(opts) == _STREAMING_PHASE2_WINDOW_MIN
        opts.workers = 8
        assert _compute_streaming_phase1_batch_size(opts) == _STREAMING_PHASE2_WINDOW_MIN
        opts.workers = max(300, (_STREAMING_PHASE2_WINDOW_MIN // 4) + 1)
        assert _compute_streaming_phase1_batch_size(opts) == max(
            _STREAMING_PHASE2_WINDOW_MIN, opts.workers * 4,
        )
        opts.streaming_phase1_batch_size = 96
        assert _compute_streaming_phase1_batch_size(opts) == 96

    def test_compute_phase1_batch_size_rejects_negative_override(self):
        from antismash.main import _compute_streaming_phase1_batch_size
        opts = type("Options", (), {
            "workers": 4,
            "streaming_phase1_batch_size": -1,
            "streaming_phase2_window_size": 0,
        })()
        with self.assertRaisesRegex(ValueError, "streaming_phase1_batch_size must be >= 0"):
            _compute_streaming_phase1_batch_size(opts)

    def test_compute_phase2_window_size(self):
        from antismash.main import _STREAMING_PHASE2_WINDOW_MIN, _compute_streaming_phase2_window_size
        opts = type("Options", (), {
            "workers": 1,
            "streaming_phase1_batch_size": 0,
            "streaming_phase2_window_size": 0,
        })()
        assert _compute_streaming_phase2_window_size(opts) == _STREAMING_PHASE2_WINDOW_MIN
        opts.workers = 8
        assert _compute_streaming_phase2_window_size(opts) == _STREAMING_PHASE2_WINDOW_MIN
        opts.workers = max(300, (_STREAMING_PHASE2_WINDOW_MIN // 4) + 1)
        assert _compute_streaming_phase2_window_size(opts) == max(
            _STREAMING_PHASE2_WINDOW_MIN, opts.workers * 4,
        )
        opts.streaming_phase2_window_size = 96
        assert _compute_streaming_phase2_window_size(opts) == 96

    def test_compute_phase2_window_size_rejects_negative_override(self):
        from antismash.main import _compute_streaming_phase2_window_size
        opts = type("Options", (), {
            "workers": 4,
            "streaming_phase1_batch_size": 0,
            "streaming_phase2_window_size": -1,
        })()
        with self.assertRaisesRegex(ValueError, "streaming_phase2_window_size must be >= 0"):
            _compute_streaming_phase2_window_size(opts)

    def test_take_record_batch_full_and_partial(self):
        from antismash.main import _take_record_batch
        records = iter([
            (SeqRecord(Seq("ACGT"), id="a"), 1),
            (SeqRecord(Seq("ACGT"), id="b"), 2),
            (SeqRecord(Seq("ACGT"), id="c"), 3),
        ])
        batch = _take_record_batch(records, 2)
        assert [record.id for record, _ in batch] == ["a", "b"]
        batch = _take_record_batch(records, 2)
        assert [record.id for record, _ in batch] == ["c"]
        assert _take_record_batch(records, 2) == []


class TestWindowedStreamingPipeline(unittest.TestCase):
    def setUp(self):
        config.build_config([], isolated=True)

    def tearDown(self):
        config.destroy_config()

    @staticmethod
    def _make_options(output_dir: str):
        return type("Options", (), {
            "minlength": -1,
            "abort_on_invalid_records": True,
            "output_dir": output_dir,
            "taxon": "bacteria",
            "cpus": 2,
            "workers": 1,
            "streaming_phase1_batch_size": 0,
            "streaming_phase2_window_size": 0,
            "output_skip_records_without_regions": False,
            "summary_gbk": False,
            "region_gbks": False,
            "zip_output": False,
            "html_taxonomy": "",
            "html_title": "",
            "output_basename": "streaming-test",
            "html_description": "",
            "debug": False,
            "version": "8.dev",
        })()

    def test_streaming_drains_multiple_phase2_windows(self):
        from antismash.main import _run_antismash_streaming

        class FakeRecord:
            def __init__(self, record_id, record_index, has_regions):
                self.id = record_id
                self.record_index = record_index
                self.original_id = ""
                self._has_regions = has_regions

            def get_regions(self):
                return [object()] if self._has_regions else []

        class FakeWriter:
            instances = []

            def __init__(self, *_args, **_kwargs):
                self.records = []
                self.finalized = False
                FakeWriter.instances.append(self)

            def write_record(self, record, _results):
                self.records.append(record.id)

            def finalize(self, _timings=None):
                self.finalized = True

        def fake_parallel(function, args, cpus=None):  # pylint: disable=unused-argument
            for argset in args:
                yield function(*argset)

        def fake_detection(record_tuple, _options):
            bio_record, record_index = record_tuple
            has_regions = not bio_record.id.endswith("_nr")
            return FakeRecord(bio_record.id, record_index, has_regions), {}, {"detection": 1.0}

        phase2_window_sizes = []

        def fake_phase2_window(phase2_inputs, _options, _picklable_options_p2, _user_workers,
                               json_writer, _all_modules, _options_layer, _data_writer, _gbk_handle,
                               _timings_by_record, _lightweight_records, _record_summaries,
                               _detection_module_names, phase2_seen, regions_count, failed_count,
                               trace_snapshot, _window_index, _window_size):
            phase2_window_sizes.append(len(phase2_inputs))
            for record_id in list(phase2_inputs):
                record, mod_results = phase2_inputs.pop(record_id)
                json_writer.write_record(record, mod_results)
                regions_count += 1
                phase2_seen += 1
            return phase2_seen, regions_count, failed_count, trace_snapshot

        records = [
            (SeqRecord(Seq("ACGT"), id="rec0"), 1),
            (SeqRecord(Seq("ACGT"), id="rec1"), 2),
            (SeqRecord(Seq("ACGT"), id="rec2"), 3),
            (SeqRecord(Seq("ACGT"), id="rec3"), 4),
            (SeqRecord(Seq("ACGT"), id="rec4"), 5),
            (SeqRecord(Seq("ACGT"), id="rec5_nr"), 6),
        ]
        accepted_ids = {record.id: (record.id, index) for record, index in records}

        with tempfile.TemporaryDirectory() as temp_dir:
            options = self._make_options(temp_dir)
            with mock.patch("antismash.main.prepare_output_directory"), \
                 mock.patch("antismash.main._preload_pfam_caches"), \
                 mock.patch("antismash.main._preload_analysis_caches"), \
                 mock.patch("antismash.main.get_all_modules", return_value=[]), \
                 mock.patch("antismash.main.html.is_enabled", return_value=False), \
                 mock.patch("antismash.main.record_processing.resolve_record_ids",
                            return_value=(accepted_ids, len(records))), \
                 mock.patch("antismash.main.record_processing.iter_accepted_records",
                            return_value=iter(records)), \
                 mock.patch("antismash.common.subprocessing.parallel_function_lazy",
                            side_effect=fake_parallel), \
                 mock.patch("antismash.main._safe_process_record_detection_streaming",
                            side_effect=fake_detection), \
                 mock.patch("antismash.main._compute_streaming_phase1_batch_size",
                            return_value=2), \
                 mock.patch("antismash.main._run_phase2_window",
                            side_effect=fake_phase2_window), \
                 mock.patch("antismash.main._compute_streaming_phase2_window_size",
                            return_value=2), \
                 mock.patch("antismash.main.serialiser.StreamingJsonWriter", FakeWriter):
                result = _run_antismash_streaming("input.gbk", options,
                                                  prefetched_metadata=[("rec0", 10)])

        assert result == 0
        assert phase2_window_sizes == [2, 2, 1]
        writer = FakeWriter.instances[0]
        assert writer.finalized is True
        assert writer.records.count("rec5_nr") == 1

    def test_streaming_writes_empty_html_assets_without_phase2(self):
        from antismash.main import _run_antismash_streaming

        class FakeRecord:
            def __init__(self, record_id, record_index):
                self.id = record_id
                self.record_index = record_index
                self.original_id = ""

            def get_regions(self):
                return []

        class FakeWriter:
            instances = []

            def __init__(self, *_args, **_kwargs):
                self.records = []
                self.finalized = False
                FakeWriter.instances.append(self)

            def write_record(self, record, _results):
                self.records.append(record.id)

            def finalize(self, _timings=None):
                self.finalized = True

        def fake_parallel(function, args, cpus=None):  # pylint: disable=unused-argument
            for argset in args:
                yield function(*argset)

        def fake_detection(record_tuple, _options):
            bio_record, record_index = record_tuple
            return FakeRecord(bio_record.id, record_index), {}, {"detection": 1.0}

        records = [(SeqRecord(Seq("ACGT"), id="rec0_nr"), 1)]
        accepted_ids = {record.id: (record.id, index) for record, index in records}

        with tempfile.TemporaryDirectory() as temp_dir:
            options = self._make_options(temp_dir)
            finalize_html = mock.Mock()
            with mock.patch("antismash.main.prepare_output_directory"), \
                 mock.patch("antismash.main._preload_pfam_caches"), \
                 mock.patch("antismash.main._preload_analysis_caches"), \
                 mock.patch("antismash.main.get_all_modules", return_value=[]), \
                 mock.patch("antismash.main.html.is_enabled", return_value=True), \
                 mock.patch("antismash.main.html.copy_template_dir"), \
                 mock.patch("antismash.main.html.find_local_antismash_js_path", return_value=None), \
                 mock.patch("antismash.main.record_processing.resolve_record_ids",
                            return_value=(accepted_ids, len(records))), \
                 mock.patch("antismash.main.record_processing.iter_accepted_records",
                            return_value=iter(records)), \
                 mock.patch("antismash.common.subprocessing.parallel_function_lazy",
                            side_effect=fake_parallel), \
                 mock.patch("antismash.main._safe_process_record_detection_streaming",
                            side_effect=fake_detection), \
                 mock.patch("antismash.outputs.html.generator.finalize_streaming_html_output",
                            finalize_html), \
                 mock.patch("antismash.main.serialiser.StreamingJsonWriter", FakeWriter):
                result = _run_antismash_streaming("input.gbk", options,
                                                  prefetched_metadata=[("rec0_nr", 10)])

            assert result == 0
            assert os.path.exists(os.path.join(temp_dir, "regions_data.js"))

        finalize_html.assert_called_once()
        assert finalize_html.call_args.args[0] == []
        assert finalize_html.call_args.args[1] == []

    def test_phase2_window_writes_detection_only_json_on_analysis_failure(self):
        from antismash.main import _RECORD_FAILED, _run_phase2_window

        class FakeRecord:
            def __init__(self, record_id):
                self.id = record_id

        def fake_parallel(_function, _args, cpus=None):  # pylint: disable=unused-argument
            yield (_RECORD_FAILED, "rec0", "boom")

        options = self._make_options("/tmp")
        json_writer = mock.Mock()
        phase2_inputs = {"rec0": (FakeRecord("rec0"), {"detection": "results"})}

        with mock.patch("antismash.common.subprocessing.parallel_function_lazy",
                        side_effect=fake_parallel):
            phase2_seen, regions_count, failed_count, trace_snapshot = _run_phase2_window(
                phase2_inputs, options, mock.Mock(), 1,
                json_writer, [], None, None, None,
                {}, [], [], set(),
                0, 0, 0, None, 1, 2,
            )

        assert phase2_seen == 1
        assert regions_count == 0
        assert failed_count == 1
        assert trace_snapshot is None
        assert phase2_inputs == {}
        json_writer.write_record.assert_called_once()
