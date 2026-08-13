"""Smoke test: default YAML loads + carries the required top-level keys."""
from __future__ import annotations


def test_default_yaml_loads():
    from ref_build.config import load_default

    cfg = load_default()
    assert isinstance(cfg, dict)
    assert cfg  # non-empty


def test_default_yaml_has_expected_top_level_sections():
    from ref_build.config import load_default

    cfg = load_default()
    for section in (
        "load_primary_and_donors",
        "census",
        "assemble",
        "export_mtx",
        "rctd_reference_build",
    ):
        assert section in cfg, f"missing section: {section}"


def test_default_yaml_has_required_stubs():
    """The four required-at-runtime keys are present in default.yaml but
    null — user must set them via CLI/user YAML."""
    from ref_build.config import REQUIRED_KEYS, load_default

    cfg = load_default()
    for key in REQUIRED_KEYS:
        assert key in cfg, f"missing required-key stub: {key}"


def test_validate_rejects_missing_required():
    """Validation raises when required keys are null."""
    import pytest
    from ref_build.config import validate

    with pytest.raises(SystemExit):
        validate({})


def test_validate_accepts_all_required():
    from ref_build.config import validate

    validate({
        "sample_id": "MH7",
        "primary_h5ad": "/tmp/primary.h5ad",
        "output_root": "/tmp/out",
        "celltype_marker_json": "/tmp/markers.json",
    })


def test_deep_update_merges_recursively():
    from ref_build.config import deep_update

    base = {"a": 1, "nested": {"x": 10, "y": 20}}
    override = {"nested": {"y": 99, "z": 100}, "b": 2}
    got = deep_update(base, override)
    assert got == {"a": 1, "b": 2, "nested": {"x": 10, "y": 99, "z": 100}}


def test_default_carries_hybrid_rule_knobs_and_primary_only():
    """Hybrid rule set knobs (TracyY123-nexus#21 comment 5137958844)."""
    from ref_build.config import load_default

    cfg = load_default()
    cen = cfg["census"]
    assert cen["donor_borrow_cap"] == 100, (
        "single cap serving Rule 1 gate, Rule 2 LOW total, and Rule 5 top-up"
    )
    assert cen["cell_min_instance"] == 20, (
        "Rule 2 branch selector + Rule 5 trigger"
    )
    assert "tumor" in cen["primary_only_celltypes"], "summary caveat §1 guard"
    assert "liver" in cen["primary_only_celltypes"], "summary caveat §1 guard"
    # `variant` and `rule1_threshold` are retired.
    assert "variant" not in cen, "census.variant was retired in #21"


def test_default_carries_rctd_reference_knobs():
    """Stage C knobs from ref-build-summary v3 lines 264-272."""
    from ref_build.config import load_default

    cfg = load_default()
    rr = cfg["rctd_reference_build"]
    assert rr["min_UMI"] == 10
    assert rr["require_int"] is True
