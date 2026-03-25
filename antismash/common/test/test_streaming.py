# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

"""Tests for the streaming pipeline functions:
    - should_use_streaming (main.py)
    - _safe_process_record_full (main.py)
"""

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
