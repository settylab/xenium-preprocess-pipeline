"""Merge-config invariant tests.

Enforces: each pipeline reads current YAML, replaces ONLY its own
top-level key (`rctd_split:`), preserves every sibling. Cross-repo
consistency gate — the sister repos (`xenium-preprocess`, `ref-build`)
MUST ship a parallel test so a regression here is caught locally.
"""
from __future__ import annotations

from pathlib import Path


def test_merge_config_creates_new_file(tmp_path: Path):
    from rctd_split._internal.merge_config import merge_config

    p = tmp_path / "config.yaml"
    merged = merge_config(p, "rctd_split", {"foo": 1})
    assert p.exists()
    assert merged == {"rctd_split": {"foo": 1}}


def test_merge_config_preserves_siblings(tmp_path: Path):
    """rctd-split writing under `rctd_split:` must not touch existing
    xenium_preprocess / ref_build sections."""
    import yaml
    from rctd_split._internal.merge_config import merge_config

    p = tmp_path / "config.yaml"
    # Simulate xenium-preprocess writing first.
    p.write_text(yaml.safe_dump({
        "xenium_preprocess": {"sample_id": "MH10", "output_root": "/tmp/o"},
        "ref_build": {"primary_h5ad": "/tmp/primary.h5ad"},
    }))

    merged = merge_config(p, "rctd_split", {"sample_id": "MH10", "run_id": "42"})
    assert merged["xenium_preprocess"] == {"sample_id": "MH10", "output_root": "/tmp/o"}
    assert merged["ref_build"] == {"primary_h5ad": "/tmp/primary.h5ad"}
    assert merged["rctd_split"] == {"sample_id": "MH10", "run_id": "42"}

    reread = yaml.safe_load(p.read_text())
    assert reread == merged


def test_merge_config_overwrites_own_key(tmp_path: Path):
    """Rerunning rctd-split replaces the rctd_split: block wholesale."""
    from rctd_split._internal.merge_config import merge_config

    p = tmp_path / "config.yaml"
    merge_config(p, "rctd_split", {"a": 1})
    merged = merge_config(p, "rctd_split", {"b": 2})
    assert merged["rctd_split"] == {"b": 2}


def test_merge_config_rejects_unknown_step_key(tmp_path: Path):
    import pytest
    from rctd_split._internal.merge_config import merge_config

    p = tmp_path / "config.yaml"
    with pytest.raises(SystemExit) as exc:
        merge_config(p, "unknown_stage", {"foo": 1})
    assert "step_key" in str(exc.value)


def test_merge_config_atomic_write(tmp_path: Path):
    """No .tmp file left after a successful write."""
    from rctd_split._internal.merge_config import merge_config

    p = tmp_path / "config.yaml"
    merge_config(p, "rctd_split", {"a": 1})
    tmp = p.with_suffix(p.suffix + ".tmp")
    assert not tmp.exists()


def test_merge_config_drops_legacy_numeric_keys(tmp_path: Path, capsys):
    """A resolved config left over from a pre-semantic-migration run
    folder can contain numeric `step1:` / `step3:` / `step4:` top-level
    keys. Re-running rctd-split against that folder MUST drop the
    legacy keys (so the file self-heals) while preserving the semantic
    siblings."""
    import yaml
    from rctd_split._internal.merge_config import merge_config

    p = tmp_path / "config.yaml"
    p.write_text(yaml.safe_dump({
        "step1": {"legacy_junk": True},
        "step3": {"legacy_junk": True},
        "step4": {"writeback_to_step1_raw": {"stale": True}},
        "xenium_preprocess": {"sample_id": "MH10"},
    }))
    merged = merge_config(p, "rctd_split", {"a": 1})
    assert set(merged) == {"xenium_preprocess", "rctd_split"}
    assert merged["xenium_preprocess"] == {"sample_id": "MH10"}
    on_disk = yaml.safe_load(p.read_text())
    assert set(on_disk) == {"xenium_preprocess", "rctd_split"}
    err = capsys.readouterr().err
    for legacy in ("step1", "step3", "step4"):
        assert legacy in err
