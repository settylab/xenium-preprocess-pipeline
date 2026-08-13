"""Smoke test: default YAML loads + carries the required top-level keys."""
from __future__ import annotations


def test_default_yaml_loads():
    from xenium_preprocess.config import load_default

    cfg = load_default()
    assert isinstance(cfg, dict)
    assert cfg  # non-empty


def test_default_yaml_has_expected_top_level_sections():
    from xenium_preprocess.config import load_default

    cfg = load_default()
    for section in (
        "proseg_to_anndata", "enrich_xenium_id", "qc_filter",
        "xenium_ranger_to_anndata", "preprocess", "split_prep", "rctd_prep",
    ):
        assert section in cfg, f"missing section: {section}"


def test_default_stages_excludes_preprocess():
    """`preprocess` is retained in VALID_STAGES for reversibility but is
    NOT part of DEFAULT_STAGES (roadmap §2.5)."""
    from xenium_preprocess.config import DEFAULT_STAGES, VALID_STAGES

    assert "preprocess" in VALID_STAGES
    assert "preprocess" not in DEFAULT_STAGES
    for new_stage in ("qc_filter", "xenium_ranger_to_anndata"):
        assert new_stage in DEFAULT_STAGES
        assert new_stage in VALID_STAGES


def test_default_stages_split_out_rctd_in():
    """Post-2026-08-12 rectification (TracyY123-nexus#26 comment
    5273369399): `rctd_prep` is BACK in DEFAULT_STAGES — Tracy still
    needs the `test_object.rds`. `split_prep` (the persisted 10x-style
    mtx bundle) stays OUT of defaults; the pipeline auto-produces it
    as a transient prerequisite for `rctd_prep` and deletes it after
    the RDS is written. Both remain in VALID_STAGES for opt-in."""
    from xenium_preprocess.config import DEFAULT_STAGES, VALID_STAGES

    assert "split_prep" in VALID_STAGES
    assert "rctd_prep" in VALID_STAGES
    assert "split_prep" not in DEFAULT_STAGES, (
        "split_prep should NOT be in DEFAULT_STAGES "
        "(mtx bundle is a transient prerequisite, not a persisted output)"
    )
    assert "rctd_prep" in DEFAULT_STAGES, (
        "rctd_prep should be in DEFAULT_STAGES "
        "(Tracy still needs test_object.rds — comment 5273369399)"
    )


def test_default_stages_rctd_runs_after_prereqs():
    """`rctd_prep` reads the raw AnnData (via a transient split_prep
    bundle), so `proseg_to_anndata` must dispatch before it."""
    from xenium_preprocess.config import DEFAULT_STAGES

    stages = list(DEFAULT_STAGES)
    assert stages.index("proseg_to_anndata") < stages.index("rctd_prep"), (
        f"unexpected order: {stages}"
    )


def test_default_stages_enrich_runs_after_xenium_ranger():
    """`enrich_xenium_id` reads from the xenium-ranger h5ad, so it
    must dispatch AFTER `xenium_ranger_to_anndata`."""
    from xenium_preprocess.config import DEFAULT_STAGES

    stages = list(DEFAULT_STAGES)
    assert "enrich_xenium_id" in stages
    assert "xenium_ranger_to_anndata" in stages
    assert stages.index("xenium_ranger_to_anndata") < stages.index(
        "enrich_xenium_id"
    ), f"unexpected order: {stages}"


def test_default_yaml_has_required_stubs():
    """The three required-at-runtime keys are present in default.yaml but
    null — user must set them via CLI/user YAML."""
    from xenium_preprocess.config import REQUIRED_KEYS, load_default

    cfg = load_default()
    for key in REQUIRED_KEYS:
        assert key in cfg, f"missing required-key stub: {key}"


def test_validate_rejects_missing_required():
    """Validation raises when a required key is null."""
    import pytest
    from xenium_preprocess.config import validate

    with pytest.raises(SystemExit):
        validate({})


def test_validate_accepts_all_required():
    from xenium_preprocess.config import validate

    validate({
        "sample_id": "MH10",
        "proseg_dir": "/tmp/proseg",
        "output_root": "/tmp/out",
    })


def test_deep_update_merges_recursively():
    from xenium_preprocess.config import deep_update

    base = {"a": 1, "nested": {"x": 10, "y": 20}}
    override = {"nested": {"y": 99, "z": 100}, "b": 2}
    got = deep_update(base, override)
    assert got == {"a": 1, "b": 2, "nested": {"x": 10, "y": 99, "z": 100}}
