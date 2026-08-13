"""`--run-id` precedence: CLI flag > $SLURM_JOB_ID > YYYYMMDD_HHMMSS."""
from __future__ import annotations


def test_run_id_cli_flag_wins(monkeypatch):
    from xenium_preprocess.cli import _resolve_run_id

    monkeypatch.setenv("SLURM_JOB_ID", "999999")
    assert _resolve_run_id("my_explicit_id") == "my_explicit_id"


def test_run_id_falls_back_to_slurm_job_id(monkeypatch):
    from xenium_preprocess.cli import _resolve_run_id

    monkeypatch.setenv("SLURM_JOB_ID", "12345")
    assert _resolve_run_id(None) == "12345"


def test_run_id_falls_back_to_timestamp(monkeypatch):
    from xenium_preprocess.cli import _resolve_run_id

    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    result = _resolve_run_id(None)
    # `YYYYMMDD_HHMMSS` — 8 digits, underscore, 6 digits.
    assert len(result) == 15
    assert result[8] == "_"
    assert result[:8].isdigit()
    assert result[9:].isdigit()


def test_run_id_empty_string_falls_through(monkeypatch):
    """An empty --run-id string should fall through to $SLURM_JOB_ID
    (not silently used as a run id — an empty folder name would break
    the layout)."""
    from xenium_preprocess.cli import _resolve_run_id

    monkeypatch.setenv("SLURM_JOB_ID", "12345")
    assert _resolve_run_id("") == "12345"


def test_run_id_landed_in_cfg_by_resolve_config(monkeypatch):
    """`_resolve_config` should stamp the resolved run_id into cfg."""
    from xenium_preprocess.cli import _resolve_config, build_parser

    monkeypatch.setenv("SLURM_JOB_ID", "42")
    parser = build_parser()
    args = parser.parse_args([
        "run",
        "--sample-id", "MH10",
        "--proseg-dir", "/tmp/proseg",
        "--output-root", "/tmp/out",
    ])
    cfg = _resolve_config(args)
    assert cfg["run_id"] == "42"


def test_run_id_cli_wins_in_cfg(monkeypatch):
    from xenium_preprocess.cli import _resolve_config, build_parser

    monkeypatch.setenv("SLURM_JOB_ID", "999")
    parser = build_parser()
    args = parser.parse_args([
        "run",
        "--sample-id", "MH10",
        "--proseg-dir", "/tmp/proseg",
        "--output-root", "/tmp/out",
        "--run-id", "my_experiment_v2",
    ])
    cfg = _resolve_config(args)
    assert cfg["run_id"] == "my_experiment_v2"
