"""Merge-config invariant: each pipeline writes ONLY its own top-level
key, preserves all siblings.

The invariant is load-bearing for the three-pipeline
(xenium-preprocess / ref-build / rctd-split) sbatch chain
(TracyY123-nexus#26 comment 5251080220 §2). A regression here silently
clobbers a sibling pipeline's config section.
"""
from __future__ import annotations

from pathlib import Path

import pytest
import yaml


def test_merge_config_writes_only_own_key(tmp_path: Path):
    from xenium_preprocess._internal.merge_config import merge_config

    resolved = tmp_path / "config.yaml"

    # xenium-preprocess writes first.
    merge_config(resolved, "xenium_preprocess", {"sample_id": "MH10", "run_id": "42"})
    got = yaml.safe_load(resolved.read_text())
    assert set(got) == {"xenium_preprocess"}
    assert got["xenium_preprocess"]["sample_id"] == "MH10"

    # Simulate ref-build writing (would happen in the sister repo).
    merge_config(resolved, "ref_build", {"flex_h5ad_path": "/some/path.h5ad"})
    got = yaml.safe_load(resolved.read_text())
    assert set(got) == {"xenium_preprocess", "ref_build"}
    assert got["xenium_preprocess"]["sample_id"] == "MH10"     # xenium_preprocess survived
    assert got["ref_build"]["flex_h5ad_path"] == "/some/path.h5ad"

    # Now xenium-preprocess re-runs — its own section is REPLACED,
    # ref_build's is preserved.
    merge_config(resolved, "xenium_preprocess", {"sample_id": "MH10", "run_id": "43"})
    got = yaml.safe_load(resolved.read_text())
    assert set(got) == {"xenium_preprocess", "ref_build"}
    assert got["xenium_preprocess"]["run_id"] == "43"          # xenium_preprocess replaced
    assert got["ref_build"]["flex_h5ad_path"] == "/some/path.h5ad"  # ref_build survived


def test_merge_config_rejects_unknown_top_level_key(tmp_path: Path):
    from xenium_preprocess._internal.merge_config import merge_config

    resolved = tmp_path / "config.yaml"
    with pytest.raises(SystemExit):
        merge_config(resolved, "unknown_stage", {"foo": 1})
    with pytest.raises(SystemExit):
        merge_config(resolved, "Xenium_Preprocess", {"foo": 1})  # case-sensitive


def test_merge_config_atomic_write_leaves_no_tmp(tmp_path: Path):
    from xenium_preprocess._internal.merge_config import merge_config

    resolved = tmp_path / "config.yaml"
    merge_config(resolved, "xenium_preprocess", {"k": "v"})
    # No stray .tmp sidecar.
    siblings = list(resolved.parent.iterdir())
    assert [p.name for p in siblings] == [resolved.name]


def test_merge_config_creates_parent_dirs(tmp_path: Path):
    from xenium_preprocess._internal.merge_config import merge_config

    deep = tmp_path / "a" / "b" / "config.yaml"
    merge_config(deep, "xenium_preprocess", {"k": "v"})
    assert deep.exists()


def test_merge_config_drops_legacy_numeric_keys(tmp_path: Path, capsys):
    """A resolved config left over from a pre-semantic-migration run
    folder can contain numeric `step1:` / `step3:` / `step4:` top-level
    keys. Re-running the current pipeline against that folder MUST
    drop the legacy keys (so the file self-heals) while preserving the
    semantic siblings."""
    from xenium_preprocess._internal.merge_config import merge_config

    resolved = tmp_path / "config.yaml"
    resolved.write_text(yaml.safe_dump({
        "step1": {"legacy_junk": True},
        "step4": {"legacy_junk": True},
        "ref_build": {"flex_h5ad_path": "/kept.h5ad"},
    }))
    merged = merge_config(resolved, "xenium_preprocess", {"sample_id": "MH10"})
    assert set(merged) == {"xenium_preprocess", "ref_build"}
    assert merged["ref_build"]["flex_h5ad_path"] == "/kept.h5ad"
    on_disk = yaml.safe_load(resolved.read_text())
    assert set(on_disk) == {"xenium_preprocess", "ref_build"}
    err = capsys.readouterr().err
    assert "step1" in err and "step4" in err
