"""Smoke test: default YAML loads + carries the required top-level keys."""
from __future__ import annotations


def test_default_yaml_loads():
    from rctd_split.config import load_default

    cfg = load_default()
    assert isinstance(cfg, dict)
    assert cfg  # non-empty


def test_default_yaml_has_expected_top_level_sections():
    from rctd_split.config import load_default

    cfg = load_default()
    for section in ("rctd_run", "split_purify", "export_mtx", "mtx_to_h5ad"):
        assert section in cfg, f"missing section: {section}"


def test_default_yaml_has_required_stubs():
    """The required-at-runtime keys are present in default.yaml but null."""
    from rctd_split.config import REQUIRED_KEYS, load_default

    cfg = load_default()
    for key in REQUIRED_KEYS:
        assert key in cfg, f"missing required-key stub: {key}"


def test_validate_rejects_missing_required():
    """Validation raises when a required key is null."""
    import pytest
    from rctd_split.config import validate

    with pytest.raises(SystemExit):
        validate({})


def test_validate_accepts_all_required():
    from rctd_split.config import validate

    validate({
        "sample_id": "MH10",
        "test_object": "/tmp/test_object.rds",
        "reference_rds": "/tmp/ref.rds",
        "output_root": "/tmp/out",
    })


def test_validate_no_longer_requires_test_object_or_reference_rds():
    """These are auto-derived from the run-folder layout by
    `rctd_split.pipeline.run`, so they are OPTIONAL at the config layer
    (settylab/TracyY123-nexus#26 comment 5260778217).
    """
    from rctd_split.config import REQUIRED_KEYS, validate

    assert "test_object" not in REQUIRED_KEYS
    assert "reference_rds" not in REQUIRED_KEYS

    # validate() accepts a cfg with sample_id + output_root only.
    validate({"sample_id": "MH10", "output_root": "/tmp/out"})


def test_validate_still_rejects_missing_sample_or_output_root():
    """The two remaining REQUIRED_KEYS are still enforced."""
    import pytest
    from rctd_split.config import validate

    with pytest.raises(SystemExit):
        validate({"sample_id": "MH10"})            # missing output_root
    with pytest.raises(SystemExit):
        validate({"output_root": "/tmp/out"})      # missing sample_id


def test_deep_update_merges_recursively():
    from rctd_split.config import deep_update

    base = {"a": 1, "nested": {"x": 10, "y": 20}}
    override = {"nested": {"y": 99, "z": 100}, "b": 2}
    got = deep_update(base, override)
    assert got == {"a": 1, "b": 2, "nested": {"x": 10, "y": 99, "z": 100}}


def test_rctd_run_defaults_match_rmd():
    """Verify the create.RCTD / run.RCTD defaults are the Rmd values,
    except for CELL_MIN_INSTANCE which the Setty Lab lowered from 25 → 20
    on 2026-07-28 (TracyY123-nexus#21) to accommodate detailed-annotation
    references with smaller per-type counts.
    """
    from rctd_split.config import load_default

    cfg = load_default()
    rr = cfg["rctd_run"]
    # Rmd lines 344-353, 362.
    assert rr["UMI_min"] == 10
    assert rr["counts_MIN"] == 10
    assert rr["UMI_min_sigma"] == 100
    assert rr["max_cores"] == 4
    # Rmd used 25 (line 351); Setty Lab default lowered to 20 per #21.
    assert rr["CELL_MIN_INSTANCE"] == 20
    assert rr["doublet_mode"] == "doublet"

    # Rmd line 512.
    assert cfg["split_purify"]["DO_purify_singlets"] is True
