"""Smoke test: the argparse tree assembles and `--help` runs."""
from __future__ import annotations


def test_top_level_help_runs():
    """`ref-build --help` — invoked in-process — exits 0 without exploding."""
    import pytest
    from ref_build.cli import build_parser

    parser = build_parser()
    with pytest.raises(SystemExit) as exc:
        parser.parse_args(["--help"])
    assert exc.value.code == 0


def test_run_subcommand_help_runs():
    """`ref-build run --help` also exits 0 (dodges an argparse subparser
    regression that xenium-preprocess had before)."""
    import pytest
    from ref_build.cli import build_parser

    parser = build_parser()
    with pytest.raises(SystemExit) as exc:
        parser.parse_args(["run", "--help"])
    assert exc.value.code == 0


def test_run_subcommand_parses_minimal_args():
    """The full CLI surface accepts a minimal invocation without raising
    during arg-parsing (validation happens later, in _resolve_config)."""
    from ref_build.cli import build_parser

    parser = build_parser()
    args = parser.parse_args([
        "run",
        "--sample-id", "MH7",
        "--primary-h5ad", "/tmp/primary.h5ad",
        "--donor-h5ad", "/tmp/donor1.h5ad",
        "--donor-h5ad", "/tmp/donor2.h5ad",
        "--celltype-marker-json", "/tmp/markers.json",
        "--tumor-type", "Bladder",
        "--output-root", "/tmp/out",
    ])
    assert args.cmd == "run"
    assert args.sample_id == "MH7"
    assert args.tumor_type == "Bladder"
    # `action="append"` gives us both donor paths.
    assert len(args.donor_h5ads) == 2


def test_census_variant_flag_no_longer_accepted():
    """The `--census-variant` flag was retired (locked spec on
    TracyY123-nexus#21 comment 5137958844). Passing it must fail loudly
    so a legacy invocation doesn't silently drop the flag."""
    import pytest
    from ref_build.cli import build_parser

    parser = build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args([
            "run",
            "--sample-id", "MH7",
            "--primary-h5ad", "/tmp/primary.h5ad",
            "--celltype-marker-json", "/tmp/markers.json",
            "--output-root", "/tmp/out",
            "--census-variant", "balanced",
        ])


def test_cell_min_instance_flag_parses():
    """New --cell-min-instance flag parses and lands on args."""
    from ref_build.cli import build_parser

    parser = build_parser()
    args = parser.parse_args([
        "run",
        "--sample-id", "MH7",
        "--primary-h5ad", "/tmp/primary.h5ad",
        "--celltype-marker-json", "/tmp/markers.json",
        "--output-root", "/tmp/out",
        "--cell-min-instance", "50",
        "--donor-borrow-cap", "200",
    ])
    assert args.cell_min_instance == 50
    assert args.donor_borrow_cap == 200


def test_rule1_threshold_alias_still_parses():
    """`--rule1-threshold` still parses (deprecated alias)."""
    from ref_build.cli import build_parser

    parser = build_parser()
    args = parser.parse_args([
        "run",
        "--sample-id", "MH7",
        "--primary-h5ad", "/tmp/primary.h5ad",
        "--celltype-marker-json", "/tmp/markers.json",
        "--output-root", "/tmp/out",
        "--rule1-threshold", "77",
    ])
    assert args.rule1_threshold == 77


def test_version_flag_reports_package_version():
    import pytest
    from ref_build import __version__
    from ref_build.cli import build_parser

    parser = build_parser()
    with pytest.raises(SystemExit) as exc:
        parser.parse_args(["--version"])
    assert exc.value.code == 0
    assert __version__  # non-empty
