"""Drive tests/test_preserve_meta.R via Rscript.

Skipped when Rscript is unavailable OR Seurat is missing — this keeps
the pytest suite runnable in a plain Python env. The R script itself is
self-contained (no SPLIT/spacexr dependency) so it runs anywhere Seurat
loads.

Exercises _preserve_meta.R (the helper split_purify.R uses to carry
proseg + Step 1 enrichment metadata onto the purified variant — fix for
TracyY123-nexus#14).
"""
from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest


_TESTS_DIR = Path(__file__).resolve().parent
_R_SCRIPT = _TESTS_DIR / "test_preserve_meta.R"


def _rscript_bin() -> str | None:
    return shutil.which("Rscript")


def _seurat_available(rscript: str) -> bool:
    """Try `Rscript -e 'library(Seurat)'` — skip the test if this fails,
    which lets the suite run in a bare R install without the pipeline's
    R stack. Prepend the user-local lib the pipeline uses so we find
    Seurat wherever the pipeline finds it."""
    probe = (
        'user_lib <- file.path(Sys.getenv("HOME"), ".claude", "r_libs", "4.4.1"); '
        'if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths())); '
        'suppressPackageStartupMessages(library(Seurat)); '
        'cat("OK")'
    )
    try:
        r = subprocess.run(
            [rscript, "--vanilla", "-e", probe],
            capture_output=True, text=True, timeout=60,
        )
    except (OSError, subprocess.TimeoutExpired):
        return False
    return r.returncode == 0 and "OK" in r.stdout


def test_preserve_meta_r_helper():
    """Run test_preserve_meta.R end-to-end; assert stdout carries the
    OK markers and the script exits 0."""
    rscript = _rscript_bin()
    if rscript is None:
        pytest.skip("Rscript not on PATH — module load fhR/4.4.1-foss-2023b first")
    if not _seurat_available(rscript):
        pytest.skip("Seurat not loadable in this R install")

    env = os.environ.copy()
    r = subprocess.run(
        [rscript, "--vanilla", str(_R_SCRIPT)],
        capture_output=True, text=True, timeout=180, env=env,
    )
    combined = r.stdout + "\n---stderr---\n" + r.stderr
    assert r.returncode == 0, (
        f"test_preserve_meta.R exited {r.returncode}:\n{combined}"
    )
    for marker in (
        "OK: all",  # assertion 1: all expected cols preserved
        "OK: per-cell values match",  # assertion 2
        "OK: RCTD-added columns not duplicated",  # assertion 3
        "OK: Seurat-managed columns intact",  # assertion 4
        "OK: purified cell set unchanged",  # assertion 5
        "OK: cells not present in unpurified get NA",  # assertion 6
        "ALL ASSERTIONS PASSED",
    ):
        assert marker in r.stdout, (
            f"missing marker {marker!r} in R stdout:\n{combined}"
        )
