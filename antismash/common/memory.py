# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Utilities for opt-in memory diagnostics during large streaming runs."""

from __future__ import annotations

import gc
import logging
import os
import resource
import sys
import tracemalloc
from typing import Any, Dict, Iterable, List, Optional, Set


_MB = 1024 * 1024


def diagnostics_enabled(options: Any) -> bool:
    """Return whether streaming memory diagnostics are enabled."""
    return bool(getattr(options, "memory_diagnostics", False))


def diagnostics_interval(options: Any) -> int:
    """Return the configured diagnostics interval, clamped to >= 1."""
    return max(1, int(getattr(options, "memory_diagnostics_interval", 100)))


def tracemalloc_enabled(options: Any) -> bool:
    """Return whether tracemalloc snapshots should be captured."""
    return bool(getattr(options, "memory_diagnostics_tracemalloc", False))


def format_mb(value: Optional[float]) -> str:
    """Format an RSS value in MB for logging."""
    if value is None:
        return "n/a"
    return f"{value:.1f} MB"


def _read_proc_status_value(pid: int, key: str) -> Optional[int]:
    """Read a kB value from /proc/<pid>/status."""
    path = f"/proc/{pid}/status"
    try:
        with open(path, "r", encoding="utf-8") as handle:
            for line in handle:
                if not line.startswith(key):
                    continue
                parts = line.split()
                if len(parts) < 2:
                    return None
                return int(parts[1]) * 1024
    except (FileNotFoundError, PermissionError, ProcessLookupError, ValueError):
        return None
    return None


def current_rss_mb(pid: Optional[int] = None) -> Optional[float]:
    """Return current RSS in MB for a PID, best-effort."""
    if pid is None:
        pid = os.getpid()
    rss_bytes = _read_proc_status_value(pid, "VmRSS:")
    if rss_bytes is not None:
        return rss_bytes / _MB
    if pid != os.getpid():
        return None
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform == "darwin":
        return peak / _MB
    return peak / 1024


def peak_rss_mb() -> float:
    """Return peak RSS for the current process in MB."""
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform == "darwin":
        return peak / _MB
    return peak / 1024


def _read_child_pids(pid: int) -> List[int]:
    """Return direct child PIDs from /proc, best-effort."""
    path = f"/proc/{pid}/task/{pid}/children"
    try:
        with open(path, "r", encoding="utf-8") as handle:
            raw = handle.read().strip()
    except (FileNotFoundError, PermissionError, ProcessLookupError):
        return []
    if not raw:
        return []
    try:
        return [int(part) for part in raw.split()]
    except ValueError:
        return []


def descendant_pids(pid: Optional[int] = None) -> List[int]:
    """Return all descendant PIDs for a process, best-effort."""
    if pid is None:
        pid = os.getpid()
    seen: Set[int] = set()
    pending = list(_read_child_pids(pid))
    descendants: List[int] = []
    while pending:
        child = pending.pop()
        if child in seen:
            continue
        seen.add(child)
        descendants.append(child)
        pending.extend(_read_child_pids(child))
    return descendants


def process_tree_stats(pid: Optional[int] = None) -> Dict[str, Optional[float]]:
    """Return current RSS for the parent process and its descendants."""
    if pid is None:
        pid = os.getpid()
    children = descendant_pids(pid)
    children_rss = 0.0
    known_children = 0
    for child in children:
        rss = current_rss_mb(child)
        if rss is None:
            continue
        children_rss += rss
        known_children += 1
    return {
        "parent_rss_mb": current_rss_mb(pid),
        "descendant_count": float(len(children)),
        "descendant_rss_mb": children_rss,
        "descendant_rss_known_count": float(known_children),
    }


def tracked_object_counts() -> Dict[str, int]:
    """Count selected GC-tracked antiSMASH objects."""
    from antismash.common.layers import RecordLayer, RegionLayer
    from antismash.common.module_results import ModuleResults
    from antismash.common.secmet import Record

    counts = {
        "Record": 0,
        "RecordLayer": 0,
        "RegionLayer": 0,
        "ModuleResults": 0,
    }
    for obj in gc.get_objects():
        try:
            if isinstance(obj, Record):
                counts["Record"] += 1
            if isinstance(obj, RecordLayer):
                counts["RecordLayer"] += 1
            if isinstance(obj, RegionLayer):
                counts["RegionLayer"] += 1
            if isinstance(obj, ModuleResults):
                counts["ModuleResults"] += 1
        except ReferenceError:
            continue
    return counts


def maybe_start_tracemalloc(options: Any) -> bool:
    """Start tracemalloc if requested."""
    if not tracemalloc_enabled(options):
        return False
    if not tracemalloc.is_tracing():
        tracemalloc.start(25)
    return True


def maybe_log_tracemalloc(label: str, previous: Any = None, limit: int = 5) -> Any:
    """Log a tracemalloc snapshot diff and return the new snapshot."""
    if not tracemalloc.is_tracing():
        return previous
    snapshot = tracemalloc.take_snapshot()
    if previous is None:
        stats: Iterable[Any] = snapshot.statistics("lineno")[:limit]
        for stat in stats:
            frame = stat.traceback[0]
            logging.info("memdiag tracemalloc %s %s:%d size=%s count=%d",
                         label, frame.filename, frame.lineno,
                         format_mb(stat.size / _MB), stat.count)
        return snapshot

    for stat in snapshot.compare_to(previous, "lineno")[:limit]:
        frame = stat.traceback[0]
        logging.info("memdiag tracemalloc %s %s:%d size_diff=%s count_diff=%+d",
                     label, frame.filename, frame.lineno,
                     format_mb(stat.size_diff / _MB), stat.count_diff)
    return snapshot
