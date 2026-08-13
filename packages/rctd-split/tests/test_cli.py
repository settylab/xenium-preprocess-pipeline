"""Smoke test: the argparse tree assembles and `--help` runs."""
from __future__ import annotations


def test_top_level_help_runs():
    import pytest
    from rctd_split.cli import build_parser

    parser = build_parser()
    with pytest.raises(SystemExit) as exc:
        parser.parse_args(["--help"])
    assert exc.value.code == 0


def test_run_subcommand_help_runs():
    import pytest
    from rctd_split.cli import build_parser

    parser = build_parser()
    with pytest.raises(SystemExit) as exc:
        parser.parse_args(["run", "--help"])
    assert exc.value.code == 0


def test_run_subcommand_parses_minimal_args():
    from rctd_split.cli import build_parser

    parser = build_parser()
    args = parser.parse_args([
        "run",
        "--sample-id", "MH10",
        "--test-object", "/tmp/test.rds",
        "--reference-rds", "/tmp/ref.rds",
        "--output-root", "/tmp/out",
    ])
    assert args.cmd == "run"
    assert args.sample_id == "MH10"


def test_version_flag_reports_package_version():
    import pytest
    from rctd_split import __version__
    from rctd_split.cli import build_parser

    parser = build_parser()
    with pytest.raises(SystemExit) as exc:
        parser.parse_args(["--version"])
    assert exc.value.code == 0
    assert __version__


def test_resolve_run_id_precedence(monkeypatch):
    """--run-id > $SLURM_JOB_ID > YYYYMMDD_HHMMSS."""
    from rctd_split.cli import _resolve_run_id

    # 1. --run-id wins.
    monkeypatch.setenv("SLURM_JOB_ID", "999")
    assert _resolve_run_id("mycustom") == "mycustom"

    # 2. $SLURM_JOB_ID wins when --run-id is None.
    monkeypatch.setenv("SLURM_JOB_ID", "12345")
    assert _resolve_run_id(None) == "12345"

    # 3. Timestamp fallback when neither is set.
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    got = _resolve_run_id(None)
    # A YYYYMMDD_HHMMSS string is 8+1+6=15 chars, all digits+underscore.
    assert len(got) == 15
    assert got[8] == "_"


def test_stages_flag_accepts_writeback_only_subset():
    """Tracy on `settylab/TracyY123-nexus#26` comment 5260916505: she
    wants to rerun ONLY the writeback stages after a completed SPLIT.
    Prove the CLI parses `--stages writeback_to_step1_raw
    celltype_writeback` cleanly."""
    from rctd_split.cli import build_parser

    parser = build_parser()
    args = parser.parse_args([
        "run",
        "--sample-id", "MH10",
        "--test-object", "/tmp/test.rds",
        "--reference-rds", "/tmp/ref.rds",
        "--output-root", "/tmp/out",
        "--stages", "writeback_to_step1_raw", "celltype_writeback",
    ])
    assert args.stages == ["writeback_to_step1_raw", "celltype_writeback"]


def test_stages_flag_rejects_unknown_stage():
    """The stage vocabulary is closed — a typo must fail argparse
    with a clear message, not silently run a smaller pipeline."""
    import pytest
    from rctd_split.cli import build_parser

    parser = build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args([
            "run",
            "--sample-id", "MH10",
            "--test-object", "/tmp/test.rds",
            "--reference-rds", "/tmp/ref.rds",
            "--output-root", "/tmp/out",
            "--stages", "typo_stage",
        ])


def test_run_id_cli_flag_parses():
    from rctd_split.cli import build_parser

    parser = build_parser()
    args = parser.parse_args([
        "run",
        "--sample-id", "MH10",
        "--test-object", "/tmp/test.rds",
        "--reference-rds", "/tmp/ref.rds",
        "--output-root", "/tmp/out",
        "--run-id", "my-run-42",
    ])
    assert args.run_id == "my-run-42"
