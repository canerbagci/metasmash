# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest
from unittest.mock import mock_open, patch

from antismash.common import cpu_diagnostics


def _stats(parent_cpu: float, descendants: dict) -> dict:
    """Build a process_tree_cpu_seconds-shaped dict for tests."""
    descendant_cpu = sum(descendants.values())
    return {
        "parent_cpu_seconds": parent_cpu,
        "descendant_cpu_seconds": descendant_cpu,
        "total_cpu_seconds": parent_cpu + descendant_cpu,
        "descendant_count": float(len(descendants)),
        "descendant_known_count": float(len(descendants)),
    }


class CpuSamplerTest(unittest.TestCase):
    def setUp(self):
        self.sampler = cpu_diagnostics.CpuSampler()
        self.wall_times = iter([100.0, 110.0, 120.0, 130.0, 140.0])
        self._wall_patcher = patch.object(
            cpu_diagnostics.time, "monotonic", side_effect=lambda: next(self.wall_times)
        )
        self._wall_patcher.start()
        self.addCleanup(self._wall_patcher.stop)

    def _sample_with(self, stats):
        with patch.object(cpu_diagnostics, "process_tree_cpu_seconds",
                          return_value=stats):
            return self.sampler.sample()

    def test_first_sample_returns_baseline_with_none_deltas(self):
        result = self._sample_with(_stats(parent_cpu=2.0,
                                          descendants={10: 3.0, 11: 4.0}))
        assert result["wall_delta"] is None
        assert result["cpu_delta"] is None
        assert result["parent_cpu_delta"] is None
        assert result["descendant_cpu_delta"] is None
        assert result["effective_cpus"] is None
        assert result["descendants"] == 2

    def test_growing_descendant_set_yields_positive_delta(self):
        # Baseline at t=100
        self._sample_with(_stats(parent_cpu=2.0,
                                 descendants={10: 3.0, 11: 4.0}))
        # New worker spawned (PID 12), existing workers progressed
        result = self._sample_with(_stats(parent_cpu=4.0,
                                          descendants={10: 8.0, 11: 9.0, 12: 6.0}))
        assert result["wall_delta"] == 10.0
        assert result["parent_cpu_delta"] == 2.0  # 4 - 2
        # descendants snapshot grew from 7 to 23 → +16
        assert result["descendant_cpu_delta"] == 16.0
        assert result["cpu_delta"] == 18.0
        assert result["effective_cpus"] == 1.8
        assert result["descendants"] == 3

    def test_pool_reap_with_parent_absorbing_worker_cpu(self):
        # Regression test for the original negative-delta bug AND for
        # the follow-up 1 double-count bug at once.
        #
        # Sample 1: 3 workers each at inclusive 50/60/70, parent at 2.0.
        # Snapshot total = 182.
        self._sample_with(_stats(parent_cpu=2.0,
                                 descendants={10: 50.0, 11: 60.0, 12: 70.0}))
        # Sample 2: pool joined. All workers reaped, so the kernel
        # folded their full lifetime CPU (180 cs) into the parent's
        # cutime+cstime. Plus the parent did 1.5 cs of own work in the
        # gap, so the parent's inclusive value becomes 2.0 + 180 + 1.5
        # = 183.5. Workers gone from the walk.
        result = self._sample_with(_stats(parent_cpu=183.5, descendants={}))
        assert result["wall_delta"] == 10.0
        # Parent's inclusive jumped by 181.5 (180 absorbed + 1.5 own).
        assert result["parent_cpu_delta"] == 181.5
        # Descendant side dropped by 180 (they vanished from the walk).
        assert result["descendant_cpu_delta"] == -180.0
        # The conserved quantity is the sum: 1.5 cs of NEW work,
        # NOT 181.5 (which would be double-counting the workers' CPU)
        # and NOT a negative number (which the pre-cutime sampler
        # produced when workers exited).
        assert result["cpu_delta"] == 1.5
        assert result["effective_cpus"] == 0.15
        assert result["descendants"] == 0

    def test_phase1_batch_does_not_double_count_at_reap(self):
        # Explicit regression test for the follow-up 1 / follow-up 2
        # interaction bug observed in slurm-3222873.err, where every
        # batch-complete sample re-attributed the workers' full lifetime
        # CPU on top of what we'd already captured per-PID during the
        # batch's progress samples — yielding effective_cpus=33800.95
        # at one phase1-batch-complete line.
        #
        # Simulate four samples spanning a single batch:
        #   t=100: baseline, 3 workers each at 0 cs
        #   t=110: workers progressed to 5 cs each (15 total new)
        #   t=120: workers progressed to 10 cs each (15 more)
        #   t=130: pool joined, parent absorbed 30 cs of worker CPU,
        #          plus did 0.1 cs of own work
        self._sample_with(_stats(parent_cpu=1.0,
                                 descendants={10: 0.0, 11: 0.0, 12: 0.0}))
        mid1 = self._sample_with(_stats(parent_cpu=1.05,
                                        descendants={10: 5.0, 11: 5.0, 12: 5.0}))
        mid2 = self._sample_with(_stats(parent_cpu=1.10,
                                        descendants={10: 10.0, 11: 10.0, 12: 10.0}))
        complete = self._sample_with(_stats(parent_cpu=1.10 + 30.0 + 0.1,
                                            descendants={}))

        # Each progress sample captured 15 cs of new work (3 workers
        # × 5 cs each) plus 0.05 cs of parent own work.
        assert abs(mid1["cpu_delta"] - 15.05) < 1e-9
        assert abs(mid2["cpu_delta"] - 15.05) < 1e-9

        # The reap-complete sample sees the parent jump by 30.0 (the
        # absorbed worker total) + 0.1 (parent own work), and the
        # descendant side drop by 30.0 (workers gone). Sum = 0.1.
        assert abs(complete["parent_cpu_delta"] - 30.1) < 1e-9
        assert complete["descendant_cpu_delta"] == -30.0
        assert abs(complete["cpu_delta"] - 0.1) < 1e-9

        # The summed cpu_delta across the three intervals must equal
        # the actual incremental work, not 2× it.
        # Real new work since baseline: 30 cs of worker work + 0.2 cs
        # of parent own work = 30.2.
        total = mid1["cpu_delta"] + mid2["cpu_delta"] + complete["cpu_delta"]
        assert abs(total - 30.2) < 1e-6
        # Pre-fix, the same scenario would have summed to ~60 cs because
        # the reap event re-counted all 30 cs of worker CPU.

    def test_two_pool_recycles_keep_summed_delta_honest(self):
        # Two batches in sequence. After both reaps, the summed
        # cpu_delta across all reported intervals should equal the
        # true new work, not 2x the reaped CPU.
        # Baseline: parent only.
        self._sample_with(_stats(parent_cpu=1.0, descendants={}))
        # Pool 1 alive: 64 workers each at 50 cs.
        a = self._sample_with(_stats(parent_cpu=1.1,
                                     descendants={i: 50.0 for i in range(100, 164)}))
        # Pool 1 reaped: parent absorbs 64 * 50 = 3200 cs.
        b = self._sample_with(_stats(parent_cpu=1.2 + 3200,
                                     descendants={}))
        # Pool 2 alive: 64 workers each at 25 cs.
        c = self._sample_with(_stats(parent_cpu=1.3 + 3200,
                                     descendants={i: 25.0 for i in range(200, 264)}))
        # Pool 2 reaped: parent absorbs another 64 * 25 = 1600 cs.
        d = self._sample_with(_stats(parent_cpu=1.4 + 3200 + 1600,
                                     descendants={}))

        # All deltas non-negative.
        for r in (a, b, c, d):
            assert r["cpu_delta"] >= 0

        # Real new work since baseline = 0.4 cs of parent own work
        # + 3200 cs (pool 1) + 1600 cs (pool 2) = 4800.4
        total = sum(r["cpu_delta"] for r in (a, b, c, d))
        assert abs(total - 4800.4) < 1e-6


class ProcessTreeCpuSecondsShapeTest(unittest.TestCase):
    def test_returned_dict_includes_summary_fields(self):
        # The sampler reads parent_cpu_seconds, descendant_cpu_seconds,
        # total_cpu_seconds, and descendant_count from the returned
        # dict. The shape test guards against accidentally renaming
        # any of those.
        result = cpu_diagnostics.process_tree_cpu_seconds()
        for key in ("parent_cpu_seconds", "descendant_cpu_seconds",
                    "total_cpu_seconds", "descendant_count",
                    "descendant_known_count"):
            assert key in result
        # total = parent + descendant (sanity check)
        assert (result["total_cpu_seconds"]
                == result["parent_cpu_seconds"] + result["descendant_cpu_seconds"])


class ReadProcStatCpuSecondsTest(unittest.TestCase):
    """Verify the /proc/<pid>/stat parser captures cutime + cstime."""

    def test_inclusive_cpu_includes_cutime_and_cstime(self):
        # Field layout per proc(5), 1-indexed:
        #   1   pid
        #   2   (comm)
        #   3   state
        #   4   ppid
        #  ... (state through cmajflt: 12 fields after comm)
        #  14   utime    -> 11th element after stripping "(comm) "
        #  15   stime    -> 12th
        #  16   cutime   -> 13th
        #  17   cstime   -> 14th
        # Build a stat line where the parenthesised comm contains a
        # space and a stray ')' to exercise the rfind(')') path.
        comm = "(hmm scan)thing)"
        # Pre-comm: pid then comm. Post-comm: 12 placeholders then
        # utime=100, stime=20, cutime=300, cstime=40, then a few extras.
        post_comm_fields = [
            "S",   # state
            "1",   # ppid
            "1",   # pgrp
            "1",   # session
            "0",   # tty_nr
            "-1",  # tpgid
            "0",   # flags
            "0", "0", "0", "0",   # minflt..cmajflt
            "100", "20", "300", "40",  # utime, stime, cutime, cstime
            "0", "0",  # priority, nice (extras to exceed the guard)
        ]
        line = f"42 {comm} " + " ".join(post_comm_fields) + "\n"
        with patch("builtins.open", mock_open(read_data=line)):
            result = cpu_diagnostics._read_proc_stat_cpu_seconds(42)
        expected = (100 + 20 + 300 + 40) / cpu_diagnostics._CLK_TCK
        assert result == expected

    def test_too_few_fields_returns_none(self):
        # A stat line missing cutime/cstime should be rejected,
        # because the inclusive read depends on having all four fields.
        line = "1 (proc) S 0 1 1 1 0 -1 0 0 0 0 0 100 20\n"
        with patch("builtins.open", mock_open(read_data=line)):
            assert cpu_diagnostics._read_proc_stat_cpu_seconds(1) is None


class CpuSamplerSubprocessReapTest(unittest.TestCase):
    """Regression test for the cutime/cstime undercount fix.

    When a worker reaps a transient subprocess between two of our
    walks, the subprocess's full lifetime CPU flows into the worker's
    cutime, and the snapshot delta on that interval should reflect the
    jump. This test exercises the case end-to-end through CpuSampler,
    so it catches both the parser change (per-PID inclusive read) and
    the snapshot-difference behaviour against an inclusive value.
    """

    def setUp(self):
        self.sampler = cpu_diagnostics.CpuSampler()
        self.wall_times = iter([100.0, 110.0, 120.0])
        self._wall_patcher = patch.object(
            cpu_diagnostics.time, "monotonic",
            side_effect=lambda: next(self.wall_times),
        )
        self._wall_patcher.start()
        self.addCleanup(self._wall_patcher.stop)

    def _sample_with(self, stats):
        with patch.object(cpu_diagnostics, "process_tree_cpu_seconds",
                          return_value=stats):
            return self.sampler.sample()

    def test_worker_cputime_jump_via_reaped_subprocess_is_captured(self):
        # Sample 1: parent + 1 worker, worker has done a tiny bit of
        # work itself (reflected in its inclusive value of 2.0).
        self._sample_with(_stats(parent_cpu=0.5, descendants={10: 2.0}))
        # Sample 2: between samples, the worker called subprocess.run()
        # for an hmmsearch that ran for 8 wall seconds and used ~64
        # core-seconds (e.g. --cpu 8 burned 8 cores for 8 seconds).
        # The worker reaped it, so the kernel folded those 64 cs into
        # the worker's cutime/cstime — the worker's inclusive value
        # jumps from 2.0 to 66.5 even though the worker spent the
        # entire interval blocked in wait().
        result = self._sample_with(_stats(parent_cpu=0.6,
                                          descendants={10: 66.5}))
        assert result["wall_delta"] == 10.0
        assert abs(result["parent_cpu_delta"] - 0.1) < 1e-9   # 0.6 - 0.5
        assert result["descendant_cpu_delta"] == 64.5  # 66.5 - 2.0
        assert abs(result["cpu_delta"] - 64.6) < 1e-9
        # 64.6 cpu-seconds in 10 wall seconds → ~6.46 effective cores,
        # which is the right kind of order of magnitude — pre-fix the
        # same scenario would have reported ~0.06 effective cores
        # because the subprocess CPU was completely invisible.
        assert abs(result["effective_cpus"] - 6.46) < 1e-6
