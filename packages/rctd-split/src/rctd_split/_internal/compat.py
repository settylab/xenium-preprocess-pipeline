"""Small compatibility shims used across stages.

Kept minimal on purpose; if something grows past ~30 lines it probably
wants its own module.
"""
from __future__ import annotations

from pathlib import Path


def sentinel_exists(path: Path, force_rerun: bool) -> bool:
    """Return True if `path` already exists and `force_rerun` is False.

    Every stage guards its work behind a call to this helper so that
    re-running the pipeline picks up where the last run left off.
    """
    return path.exists() and not force_rerun
