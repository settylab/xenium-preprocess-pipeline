"""Merge-config invariant tests.

Enforces: each pipeline reads current YAML, replaces ONLY its own
top-level key (`step4:`), preserves every sibling. Cross-repo consistency
gate — the sister repos (`xenium-preprocess`, `ref-build`) MUST ship a
parallel test so a regression here is caught locally.
"""
from __future__ import annotations

from pathlib import Path


def test_merge_config_creates_new_file(tmp_path: Path):
    from rctd_split._internal.merge_config import merge_config

    p = tmp_path / "resolved_config.yaml"
    merged = merge_config(p, "step4", {"foo": 1})
    assert p.exists()
    assert merged == {"step4": {"foo": 1}}


def test_merge_config_preserves_siblings(tmp_path: Path):
    """Step 4 writing under `step4:` must not touch existing step1/step3."""
    import yaml
    from rctd_split._internal.merge_config import merge_config

    p = tmp_path / "resolved_config.yaml"
    # Simulate step 1 writing first.
    p.write_text(yaml.safe_dump({
        "step1": {"sample_id": "MH10", "output_root": "/tmp/o"},
        "step3": {"primary_h5ad": "/tmp/primary.h5ad"},
    }))

    merged = merge_config(p, "step4", {"sample_id": "MH10", "run_id": "42"})
    assert merged["step1"] == {"sample_id": "MH10", "output_root": "/tmp/o"}
    assert merged["step3"] == {"primary_h5ad": "/tmp/primary.h5ad"}
    assert merged["step4"] == {"sample_id": "MH10", "run_id": "42"}

    reread = yaml.safe_load(p.read_text())
    assert reread == merged


def test_merge_config_overwrites_own_key(tmp_path: Path):
    """Rerunning step 4 replaces the step4: block wholesale."""
    from rctd_split._internal.merge_config import merge_config

    p = tmp_path / "resolved_config.yaml"
    merge_config(p, "step4", {"a": 1})
    merged = merge_config(p, "step4", {"b": 2})
    assert merged["step4"] == {"b": 2}


def test_merge_config_rejects_unknown_step_key(tmp_path: Path):
    import pytest
    from rctd_split._internal.merge_config import merge_config

    p = tmp_path / "resolved_config.yaml"
    with pytest.raises(SystemExit) as exc:
        merge_config(p, "step2", {"foo": 1})
    assert "step_key" in str(exc.value)


def test_merge_config_atomic_write(tmp_path: Path):
    """No .tmp file left after a successful write."""
    from rctd_split._internal.merge_config import merge_config

    p = tmp_path / "resolved_config.yaml"
    merge_config(p, "step4", {"a": 1})
    tmp = p.with_suffix(p.suffix + ".tmp")
    assert not tmp.exists()
