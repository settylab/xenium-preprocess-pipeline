"""Merge-config invariant: each pipeline writes ONLY its own top-level
key, preserves all siblings.

The invariant is load-bearing for the three-pipeline (step1/step3/step4)
sbatch chain (TracyY123-nexus#26 comment 5251080220 §2). A regression
here silently clobbers a sibling pipeline's config section.
"""
from __future__ import annotations

from pathlib import Path

import pytest
import yaml


def test_merge_config_writes_only_own_key(tmp_path: Path):
    from xenium_preprocess._internal.merge_config import merge_config

    resolved = tmp_path / "config.yaml"

    # Step 1 writes first.
    merge_config(resolved, "step1", {"sample_id": "MH10", "run_id": "42"})
    got = yaml.safe_load(resolved.read_text())
    assert set(got) == {"step1"}
    assert got["step1"]["sample_id"] == "MH10"

    # Simulate step 3 writing (would happen in the sister repo).
    merge_config(resolved, "step3", {"flex_h5ad_path": "/some/path.h5ad"})
    got = yaml.safe_load(resolved.read_text())
    assert set(got) == {"step1", "step3"}
    assert got["step1"]["sample_id"] == "MH10"     # step1 survived
    assert got["step3"]["flex_h5ad_path"] == "/some/path.h5ad"

    # Now step 1 re-runs — its own section is REPLACED, step 3's is preserved.
    merge_config(resolved, "step1", {"sample_id": "MH10", "run_id": "43"})
    got = yaml.safe_load(resolved.read_text())
    assert set(got) == {"step1", "step3"}
    assert got["step1"]["run_id"] == "43"          # step1 replaced
    assert got["step3"]["flex_h5ad_path"] == "/some/path.h5ad"  # step3 survived


def test_merge_config_rejects_unknown_top_level_key(tmp_path: Path):
    from xenium_preprocess._internal.merge_config import merge_config

    resolved = tmp_path / "config.yaml"
    with pytest.raises(SystemExit):
        merge_config(resolved, "step2", {"foo": 1})
    with pytest.raises(SystemExit):
        merge_config(resolved, "Step1", {"foo": 1})  # case-sensitive


def test_merge_config_atomic_write_leaves_no_tmp(tmp_path: Path):
    from xenium_preprocess._internal.merge_config import merge_config

    resolved = tmp_path / "config.yaml"
    merge_config(resolved, "step1", {"k": "v"})
    # No stray .tmp sidecar.
    siblings = list(resolved.parent.iterdir())
    assert [p.name for p in siblings] == [resolved.name]


def test_merge_config_creates_parent_dirs(tmp_path: Path):
    from xenium_preprocess._internal.merge_config import merge_config

    deep = tmp_path / "a" / "b" / "config.yaml"
    merge_config(deep, "step1", {"k": "v"})
    assert deep.exists()
