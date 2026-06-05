# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest
from unittest import mock
from unittest.mock import patch

from antismash.common import subprocessing
from antismash.common.subprocessing import hmmscan as hmmscan_module
from antismash.common.subprocessing import hmmscan

from .helpers import DummyConfig, DummyResult


@patch.object(hmmscan, "get_config", return_value=DummyConfig(paths={"hmmscan": "/some/non-path"}))
class TestFailureMessaging(unittest.TestCase):
    def check(self, expected, stdout="", stderr=""):
        result = DummyResult(stdout=stdout, stderr=stderr, return_code=1)
        with patch.object(hmmscan, "execute", return_value=result) as mocked:
            with patch.object(hmmscan, "run_hmmscan_help", return_value="# hmmscan"):
                with self.assertRaisesRegex(RuntimeError, f"hmmscan returned .*: '{expected}'"):
                    subprocessing.run_hmmscan("reference.hmm", ">query\nAAAA")
            mocked.assert_called_once()

    def test_stderr_before_stdout(self, dummy_config):
        assert dummy_config.executables.hmmscan
        stdout = "some output"
        stderr = "some error text"
        # no stdout should use stderr
        self.check(expected=stderr, stderr=stderr)
        # no stderr should use stdout
        self.check(expected=stdout, stdout=stdout)
        # both should use stderr
        self.check(expected=stderr, stdout=stdout, stderr=stderr)
        # and with neither the default should be used
        self.check(expected="unknown error")

    def test_empty_lines_skipped(self, dummy_config):
        assert dummy_config.executables.hmmscan
        lines = ["first", "second"]
        self.check(expected="first", stderr="\n".join(lines))
        # leading full whitespace lines should be skipped
        self.check(expected="first", stderr="\n".join(["", " "] + lines))

    def test_specific_lines(self, dummy_config):
        assert dummy_config.executables.hmmscan
        # any line starting with "Error: " is important and should be used
        error = "Error: some error message"
        lines = [" ", "first", error, "next", "last"]
        self.check(expected=f"{error} next", stderr="\n".join(lines))


def test_run_hmmscan_dispatches_to_pyhmmer_when_configured():
    sentinel = ["pyhmmer-result"]
    with mock.patch.object(hmmscan_module, "get_config") as cfg, \
         mock.patch("antismash.common.subprocessing._pyhmmer_backend.run_hmmscan_pyhmmer",
                    return_value=sentinel) as backend:
        cfg.return_value.hmmer_engine = "pyhmmer"
        result = hmmscan_module.run_hmmscan("/data/Pfam-A.hmm", ">q\nMAGIC", opts=["--cut_tc"])
    backend.assert_called_once_with("/data/Pfam-A.hmm", ">q\nMAGIC", opts=["--cut_tc"])
    assert result is sentinel


def test_run_hmmscan_uses_subprocess_by_default():
    fake = mock.Mock(); fake.successful.return_value = True; fake.stdout = ""
    with mock.patch.object(hmmscan_module, "run_hmmscan_help", return_value=" --cpu "), \
         mock.patch.object(hmmscan_module, "get_config") as cfg, \
         mock.patch.object(hmmscan_module, "get_effective_cpus", return_value=1), \
         mock.patch.object(hmmscan_module, "execute", return_value=fake) as ex:
        cfg.return_value.executables.hmmscan = "hmmscan"
        cfg.return_value.hmmer_engine = "subprocess"
        hmmscan_module.run_hmmscan("/data/Pfam-A.hmm", ">q\nMAGIC")
    ex.assert_called_once()
