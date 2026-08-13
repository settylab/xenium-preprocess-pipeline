"""Auto-derivation of `test_object` / `reference_rds` from the run-folder
layout when the caller doesn't pass them explicitly
(settylab/TracyY123-nexus#26 comment 5260778217).
"""
from __future__ import annotations


def _touch(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"")


def _base_cfg(tmp_path):
    return {
        "sample_id": "MH8_2",
        "run_id": "demo_v1",
        "output_root": str(tmp_path),
        "r_lib_paths": [],
    }


def _expected_paths(tmp_path, sample_id="MH8_2", run_id="demo_v1"):
    rctd_dir = tmp_path / sample_id / f"{sample_id}_{run_id}" / "rctd"
    return (
        rctd_dir / f"{sample_id}_test_object.rds",
        rctd_dir / f"{sample_id}_reference.rds",
    )


def test_pipeline_run_auto_derives_from_layout_when_unset(tmp_path, capsys):
    """When cfg has no test_object / reference_rds, `pipeline.run` MUST
    resolve them to the run-folder-layout paths — matching the paths
    step 1 and step 3 write to in the same run folder.
    """
    from rctd_split.pipeline import run

    exp_test, exp_ref = _expected_paths(tmp_path)
    _touch(exp_test)
    _touch(exp_ref)

    cfg = _base_cfg(tmp_path)
    # No test_object / reference_rds → auto-derive.
    # Empty stages list to avoid actually running any R/Python stage.
    exit_code = run(cfg, stages=[], argv=["rctd-split"])
    assert exit_code == 0

    stdout = capsys.readouterr().out
    assert f"test_object:   {exp_test}" in stdout
    assert f"reference_rds: {exp_ref}" in stdout
    assert "test_object   source: layout" in stdout
    assert "reference_rds source: layout" in stdout


def test_pipeline_run_honors_explicit_test_object_and_reference(tmp_path, capsys):
    """Explicit paths in cfg take precedence over auto-derivation
    (backward compat + foreign-input experiments)."""
    from rctd_split.pipeline import run

    foreign_test = tmp_path / "foreign" / "custom_test.rds"
    foreign_ref = tmp_path / "foreign" / "custom_ref.rds"
    _touch(foreign_test)
    _touch(foreign_ref)

    cfg = _base_cfg(tmp_path)
    cfg["test_object"] = str(foreign_test)
    cfg["reference_rds"] = str(foreign_ref)

    exit_code = run(cfg, stages=[], argv=["rctd-split"])
    assert exit_code == 0

    stdout = capsys.readouterr().out
    assert f"test_object:   {foreign_test}" in stdout
    assert f"reference_rds: {foreign_ref}" in stdout
    assert "test_object   source: config" in stdout
    assert "reference_rds source: config" in stdout


def test_pipeline_run_fails_loud_when_layout_derived_test_object_missing(tmp_path):
    """Auto-derive + the derived file doesn't exist → SystemExit with a
    message that names both the derived path AND the override flag."""
    import pytest
    from rctd_split.pipeline import run

    exp_test, exp_ref = _expected_paths(tmp_path)
    # reference exists but test_object does NOT.
    _touch(exp_ref)

    cfg = _base_cfg(tmp_path)
    with pytest.raises(SystemExit) as exc:
        run(cfg, stages=[], argv=["rctd-split"])
    msg = str(exc.value)
    assert str(exp_test) in msg
    assert "layout-derived" in msg
    assert "--test-object" in msg


def test_pipeline_run_fails_loud_when_layout_derived_reference_missing(tmp_path):
    """Symmetric to the test_object case for the reference RDS."""
    import pytest
    from rctd_split.pipeline import run

    exp_test, exp_ref = _expected_paths(tmp_path)
    _touch(exp_test)  # test_object exists; reference does NOT.

    cfg = _base_cfg(tmp_path)
    with pytest.raises(SystemExit) as exc:
        run(cfg, stages=[], argv=["rctd-split"])
    msg = str(exc.value)
    assert str(exp_ref) in msg
    assert "layout-derived" in msg
    assert "--reference-rds" in msg


def test_pipeline_run_fails_loud_when_explicit_test_object_missing(tmp_path):
    """Caller passed --test-object explicitly but the file doesn't
    exist — the message MUST reflect the explicit source (not
    misattribute to the layout convention)."""
    import pytest
    from rctd_split.pipeline import run

    _, exp_ref = _expected_paths(tmp_path)
    _touch(exp_ref)  # layout-derived reference exists.

    cfg = _base_cfg(tmp_path)
    cfg["test_object"] = str(tmp_path / "does-not-exist.rds")

    with pytest.raises(SystemExit) as exc:
        run(cfg, stages=[], argv=["rctd-split"])
    msg = str(exc.value)
    assert "does-not-exist.rds" in msg
    assert "source: explicit" in msg
