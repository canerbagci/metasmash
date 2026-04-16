# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Utilities for opt-in CPU-load diagnostics during large streaming runs.

Mirrors :mod:`antismash.common.memory` (the ``memdiag`` facility): reads
``/proc/<pid>/stat`` for the parent process and its descendants, computes
the effective number of CPUs used between successive samples, and exposes
helpers that emit ``cpudiag parent <label> ...`` log lines parallel to the
existing ``memdiag parent <label> ...`` lines.

The point of this module is to make stretches of single-core time during
streaming runs visible in the log: when ``effective_cpus`` between two
samples drops well below ``options.cpus``, that window was bottlenecked on
a single thread (typically Biopython parsing on the parent process).

Accounting model
----------------
The sampler uses *snapshot-difference* over the live process tree.
Every call to :meth:`CpuSampler.sample` walks the parent process and
its descendants, sums each PID's *inclusive* CPU time (own
``utime + stime`` plus ``cutime + cstime`` — i.e. transitively
including any descendants the PID has already waited on), and reports
the delta against the previous snapshot.

Inclusive reads make the snapshot total monotonic across worker-pool
reaps, which is the trick that lets a simple snapshot-difference
sampler work without ever going negative or double-counting:

* When a worker exits between samples, it vanishes from our walk and
  its contribution to the snapshot drops by ``worker.inclusive``. But
  the kernel folds the worker's full lifetime CPU into the parent's
  ``cutime + cstime`` at reap time, so the parent's contribution to
  the next snapshot rises by exactly the same amount. Net change to
  the snapshot total from this single transition is zero — no
  negative delta, no lost CPU.
* When a worker reaps a transient ``hmmsearch``/``hmmscan``/``prodigal``
  grandchild between two of our walks, the grandchild's lifetime CPU
  flows into the worker's ``cutime + cstime`` and the worker's
  inclusive value rises by exactly that amount. We never had to
  observe the grandchild directly: its CPU is captured via the
  worker's next snapshot reading.

Two known limitations:

* CPU consumed by a descendant that **both starts and exits between
  two consecutive walks of an ancestor that we never observe**
  is lost. In practice this is bounded by sample-interval and is
  small relative to the captured total; the run-total
  ``/usr/bin/time -v`` output remains the gold standard for absolute
  numbers.
* Boundary samples (``phase1-batch-start``, ``phase1-batch-complete``,
  ``phase2-window-start``, ``phase2-window-complete``,
  ``streaming-complete``) are taken when no worker pool is alive, so
  their ``descendants=0`` and ``effective_cpus`` close to ``1.0`` is
  *correct*: only the parent process ran in those gaps. Use
  ``parent_cpu_delta`` on those lines to read off the single-threaded
  parent CPU burned in the gap. Note that across a worker-pool reap
  the per-side ``parent_cpu_delta`` and ``descendant_cpu_delta``
  components can have opposite signs (CPU "moves" from the descendant
  side of the ledger to the parent side as workers are reaped). Their
  sum, ``cpu_delta``, is the conserved non-negative quantity.
"""

from __future__ import annotations

import logging
import os
import resource
import time
from typing import Any, Dict, Optional

from antismash.common import memory as memory_diagnostics


_CLK_TCK: float = float(os.sysconf("SC_CLK_TCK")) if hasattr(os, "sysconf") else 100.0


def _self_cpu_seconds_via_rusage() -> float:
    """Return total CPU seconds (user + system) for the calling process.

    This is the cross-platform fallback used when ``/proc/<pid>/stat`` is
    not available (e.g. on macOS dev hosts). On Linux production hosts the
    /proc reader is preferred because it also works for descendant PIDs.
    """
    usage = resource.getrusage(resource.RUSAGE_SELF)
    return usage.ru_utime + usage.ru_stime


def diagnostics_enabled(options: Any) -> bool:
    """Return whether streaming CPU diagnostics are enabled."""
    return bool(getattr(options, "cpu_diagnostics", False))


def diagnostics_interval(options: Any) -> int:
    """Return the configured CPU diagnostics interval, clamped to >= 1."""
    return max(1, int(getattr(options, "cpu_diagnostics_interval", 100)))


def _read_proc_stat_cpu_seconds(pid: int) -> Optional[float]:
    """Return inclusive CPU seconds for a pid via /proc/<pid>/stat.

    "Inclusive" means: this PID's own user+system CPU plus the
    user+system CPU of every descendant the PID has already waited on
    (recursively, since a child's own ``cutime/cstime`` rolls up into
    its parent's at exit time). Mirrors the accounting that
    ``/usr/bin/time -v`` performs at the very top of the process tree,
    applied per PID so that subprocess CPU is not lost when a transient
    grandchild starts and exits between two of our sample intervals.

    Returns ``None`` if the stat file cannot be read or parsed (e.g. the
    process exited between sampling steps, or we're not on Linux).
    """
    path = f"/proc/{pid}/stat"
    try:
        with open(path, "r", encoding="utf-8") as handle:
            data = handle.read()
    except (FileNotFoundError, PermissionError, ProcessLookupError):
        return None
    # The 'comm' field is wrapped in parentheses and may itself contain
    # spaces or parentheses; split on the *last* ')' to be safe.
    rparen = data.rfind(")")
    if rparen < 0:
        return None
    rest = data[rparen + 2:].split()
    # Per proc(5), after '(comm) ' the fields are: state, ppid, pgrp,
    # session, tty_nr, tpgid, flags, minflt, cminflt, majflt, cmajflt,
    # utime(11), stime(12), cutime(13), cstime(14), ...
    # cutime/cstime are the inclusive total CPU of all reaped
    # descendants (recursively) — including them is what makes us
    # capture transient subprocesses we never directly observed.
    if len(rest) < 15:
        return None
    try:
        utime_jiffies = int(rest[11])
        stime_jiffies = int(rest[12])
        cutime_jiffies = int(rest[13])
        cstime_jiffies = int(rest[14])
    except ValueError:
        return None
    return (utime_jiffies + stime_jiffies
            + cutime_jiffies + cstime_jiffies) / _CLK_TCK


def process_tree_cpu_seconds(pid: Optional[int] = None) -> Dict[str, float]:
    """Return summed CPU seconds for the parent process and its descendants.

    Uses :func:`antismash.common.memory.descendant_pids` to discover the
    process tree, so the descendant-walking logic stays consistent with
    the existing memdiag accounting. Each PID's CPU time is read
    inclusively (own ``utime + stime`` plus ``cutime + cstime``) so that
    transient subprocesses reaped between samples are not lost. On
    hosts without ``/proc`` (e.g. macOS) the parent's CPU time falls
    back to ``resource.getrusage``; descendants are unavailable in
    that case.
    """
    if pid is None:
        pid = os.getpid()
    parent_cpu_proc = _read_proc_stat_cpu_seconds(pid)
    if parent_cpu_proc is None and pid == os.getpid():
        parent_cpu = _self_cpu_seconds_via_rusage()
    else:
        parent_cpu = parent_cpu_proc or 0.0
    descendants = memory_diagnostics.descendant_pids(pid)
    descendant_cpu = 0.0
    known = 0
    for child in descendants:
        cpu = _read_proc_stat_cpu_seconds(child)
        if cpu is None:
            continue
        descendant_cpu += cpu
        known += 1
    return {
        "parent_cpu_seconds": parent_cpu,
        "descendant_cpu_seconds": descendant_cpu,
        "total_cpu_seconds": parent_cpu + descendant_cpu,
        "descendant_count": float(len(descendants)),
        "descendant_known_count": float(known),
    }


class CpuSampler:
    """Snapshot-difference CPU sampler over the process tree.

    Walks the live process tree on each :meth:`sample` call, sums each
    PID's *inclusive* CPU time (own ``utime + stime`` plus
    ``cutime + cstime``, i.e. transitively including any reaped
    descendants), and reports the delta against the previous snapshot.
    Inclusive reads make the sum monotonic across worker-pool reaps:
    a worker's lifetime CPU stays accounted for as it transitions from
    "live in our walk" to "rolled into the parent's ``cutime``", so
    the snapshot total never drops and the reported delta is never
    negative.

    The first call records a baseline and returns ``effective_cpus=None``.
    Subsequent calls return the delta in CPU seconds (parent +
    descendants) divided by the wall-clock delta, which gives the
    effective number of CPUs used during the interval between samples.
    Compare against ``options.cpus`` to spot single-core stalls.

    The ``parent_cpu_delta`` and ``descendant_cpu_delta`` fields are
    direct snapshot diffs of the parent's and descendants' inclusive
    sums respectively. Across a reap, the descendant side can go
    negative on its own with the parent side correspondingly positive
    — that's CPU "moving" from one side of the ledger to the other
    rather than new work — but their sum (``cpu_delta``) is the
    conserved, non-negative quantity.
    """

    def __init__(self) -> None:
        self._previous_total: Optional[float] = None
        self._previous_parent: Optional[float] = None
        self._previous_descendant: Optional[float] = None
        self._previous_wall: Optional[float] = None

    def reset(self) -> None:
        """Forget the previous sample so the next call sets a new baseline."""
        self.__init__()

    def sample(self) -> Dict[str, Any]:
        """Take a sample and return a dict describing it relative to the previous one."""
        stats = process_tree_cpu_seconds()
        wall_now = time.monotonic()
        total_now: float = stats["total_cpu_seconds"]
        parent_now: float = stats["parent_cpu_seconds"]
        descendant_now: float = stats["descendant_cpu_seconds"]
        descendant_count = int(stats["descendant_count"])

        first_sample = self._previous_wall is None
        wall_delta: Optional[float] = None
        cpu_delta: Optional[float] = None
        parent_cpu_delta: Optional[float] = None
        descendant_cpu_delta: Optional[float] = None
        effective: Optional[float] = None
        if not first_sample:
            assert self._previous_wall is not None  # for type-checker
            assert self._previous_total is not None
            assert self._previous_parent is not None
            assert self._previous_descendant is not None
            wall_delta = wall_now - self._previous_wall
            cpu_delta = total_now - self._previous_total
            parent_cpu_delta = parent_now - self._previous_parent
            descendant_cpu_delta = descendant_now - self._previous_descendant
            if wall_delta > 0:
                effective = cpu_delta / wall_delta

        self._previous_wall = wall_now
        self._previous_total = total_now
        self._previous_parent = parent_now
        self._previous_descendant = descendant_now

        return {
            "wall_delta": wall_delta,
            "cpu_delta": cpu_delta,
            "parent_cpu_delta": parent_cpu_delta,
            "descendant_cpu_delta": descendant_cpu_delta,
            "effective_cpus": effective,
            "descendants": descendant_count,
        }


# Module-level default sampler used by the helpers in main.py so callers
# don't have to thread an instance through the streaming orchestration.
_DEFAULT_SAMPLER = CpuSampler()


def default_sampler() -> CpuSampler:
    """Return the module-level default sampler."""
    return _DEFAULT_SAMPLER


def reset_default_sampler() -> None:
    """Forget any prior baseline on the default sampler."""
    _DEFAULT_SAMPLER.reset()


def format_seconds(value: Optional[float]) -> str:
    """Format a CPU/wall delta in seconds for logging."""
    if value is None:
        return "n/a"
    return f"{value:.2f}s"


def format_effective_cpus(value: Optional[float]) -> str:
    """Format an effective-CPUs value for logging."""
    if value is None:
        return "n/a"
    return f"{value:.2f}"


def log_parent_sample(label: str, options: Any,
                      extra: Optional[Dict[str, Any]] = None) -> None:
    """Take a sample on the default sampler and emit a ``cpudiag parent ...`` log line.

    No-op if ``options.cpu_diagnostics`` is false. The set of fields and
    the formatting convention deliberately mirror
    :func:`antismash.common.memory` style so the two facilities can be
    grepped together.
    """
    if not diagnostics_enabled(options):
        return
    sample = _DEFAULT_SAMPLER.sample()
    configured_cpus = getattr(options, "cpus", None)
    details = [
        f"effective_cpus={format_effective_cpus(sample['effective_cpus'])}",
        f"cpu_delta={format_seconds(sample['cpu_delta'])}",
        f"parent_cpu_delta={format_seconds(sample['parent_cpu_delta'])}",
        f"descendant_cpu_delta={format_seconds(sample['descendant_cpu_delta'])}",
        f"wall_delta={format_seconds(sample['wall_delta'])}",
        f"descendants={sample['descendants']}",
        f"configured_cpus={configured_cpus if configured_cpus is not None else '?'}",
    ]
    if extra:
        details.extend(f"{key}={value}" for key, value in extra.items())
    logging.info("cpudiag parent %s %s", label, " ".join(details))
