"""Command-line entrypoint for ref-build.

Exposes `ref-build run …` via the `[project.scripts]` entry point in
pyproject.toml. Design mirrors xenium-preprocess.cli — argparse builds
the surface, the config machinery in ref_build.config layers default
YAML + user YAML + CLI overrides, and the pipeline stage-loop in
ref_build.pipeline runs the requested subset.

`--run-id` precedence (bound HERE so the pipeline module sees an
already-resolved id):
    `--run-id` > `$SLURM_JOB_ID` > `YYYYMMDD_HHMMSS` timestamp fallback.
"""
from __future__ import annotations

import argparse
import os
import time
from pathlib import Path

from ref_build import __version__
from ref_build.config import (
    VALID_STAGES,
    deep_update,
    load_default,
    load_yaml,
    validate,
)


def _resolve_run_id(cli_run_id: str | None) -> str:
    """`--run-id` > `$SLURM_JOB_ID` > `YYYYMMDD_HHMMSS` timestamp.

    Public so tests can pin the precedence directly.
    """
    if cli_run_id:
        return str(cli_run_id)
    env = os.environ.get("SLURM_JOB_ID")
    if env:
        return str(env)
    return time.strftime("%Y%m%d_%H%M%S")


def _str2bool(v):
    if isinstance(v, bool):
        return v
    s = str(v).strip().lower()
    if s in {"true", "1", "yes", "y", "on"}:
        return True
    if s in {"false", "0", "no", "n", "off"}:
        return False
    raise argparse.ArgumentTypeError(f"expected boolean, got {v!r}")


def _add_run_args(p: argparse.ArgumentParser) -> None:
    """Attach every `run` subcommand flag to `p`."""
    p.add_argument("--config", type=Path, default=None,
                   help="User config YAML (overrides config/default.yaml entries).")
    p.add_argument("--stages", nargs="+", choices=VALID_STAGES,
                   default=list(VALID_STAGES),
                   help="Which stages to run.")
    p.add_argument("--run-id", default=None,
                   help="Explicit run identifier. Precedence: --run-id > "
                        "$SLURM_JOB_ID > YYYYMMDD_HHMMSS. Names the run "
                        "folder <output_root>/<sample>/<sample>_<run_id>/.")
    p.add_argument("--flex-h5ad", type=Path, default=None,
                   help="Path to the flex-preprocessed scRNA h5ad (reference-"
                        "in-place). The pipeline records this path verbatim "
                        "under `step3.flex_h5ad_path` in the merged "
                        "config.yaml — no copy, no symlink. When "
                        "--primary-h5ad is not passed the flex path is used "
                        "as the primary input.")

    # I/O
    p.add_argument("--sample-id",
                   help="Primary sample identifier (e.g. SAMPLE1). Names the output subdirectory + .rds.")
    p.add_argument("--primary-h5ad", type=Path,
                   help="Path to the primary's preprocessed scRNA h5ad. Historically "
                        "produced by step 2 (flex-preprocess), which was deprecated on "
                        "2026-08-10 (internal issue review); typically manual scRNA "
                        "prep now. Input contract: .obs[--celltype-col] populated; counts "
                        "in layers['counts'] / layers['raw_count'] / .X (first present "
                        "wins) and MUST be integer.")
    p.add_argument("--donor-h5ad", type=Path, action="append", dest="donor_h5ads",
                   default=None,
                   help="Path to a donor's preprocessed scRNA h5ad (same input contract "
                        "as --primary-h5ad). Pass repeatedly for multiple donors (must "
                        "share the primary's tumor type).")
    p.add_argument("--fallback-donor-h5ad", type=Path, action="append",
                   dest="fallback_donor_h5ads", default=None,
                   help="Path to a FALLBACK donor h5ad (Rule 5; locked spec on "
                        "(internal issue review). Pass repeatedly "
                        "for multiple fallbacks. Consulted when the Rule 2 LOW "
                        "total (primary + donor supplement) is still below "
                        "--cell-min-instance for a target-list celltype; the "
                        "assemble stage tops up from fallback to reach "
                        "--donor-borrow-cap. Intended source for immune cells "
                        "(B/Plasma, T/NK, rbc) that flex data commonly misses "
                        "across an entire sample set. Column resolution on the "
                        "fallback h5ad is auto-detected.")
    p.add_argument("--celltype-marker-json", type=Path,
                   help="Path to the marker-gene JSON declaring the expected celltype set. "
                        "Same shape as step 1's --global-non-tumor-json.")
    p.add_argument("--celltype-target-list", type=Path, default=None,
                   help="Path to a JSON that declares the CANONICAL celltype "
                        "list to enforce completeness over (Rule 5). Keys are "
                        "celltype names (with optional `_marker`/`_markers` "
                        "suffix). If omitted, --celltype-marker-json's keys "
                        "are used. Useful when the marker JSON's keys don't "
                        "match the labeling JSON one-to-one.")
    p.add_argument("--celltype-target-key", default=None,
                   help="Top-level key to read from a NESTED target-list JSON "
                        "of shape `{tissue: {celltype: markers}}`. Omit for a "
                        "flat `{celltype: markers}` JSON (the common case).")
    p.add_argument("--tumor-type", default=None,
                   help="Cosmetic label for logging + provenance (e.g. Bladder, Prostate, "
                        "Breast, LiverBenign).")
    p.add_argument("--output-root", type=Path,
                   help="Root output directory. Per-sample results land at "
                        "<output_root>/<sample_id>/.")
    p.add_argument("--force-rerun", action="store_true",
                   help="Re-run all stages even if sentinel outputs exist.")

    # load_primary_and_donors
    p.add_argument("--celltype-col", default=None,
                   help="Per-cell celltype-annotation column in .obs "
                        "(default: Final_level1_celltype_annotation).")

    # census
    p.add_argument("--donor-borrow-cap", type=int, default=None,
                   help="Cap for the user rule 1 (primary_count above this → "
                        "primary_only), Rule 2 LOW (total donor supplement), "
                        "and Rule 5 (fallback top-up target). Default 100 "
                        "(internal issue review).")
    p.add_argument("--cell-min-instance", type=int, default=None,
                   help="Cell minimum instance for Rule 2 branch selection "
                        "and Rule 5 trigger. primary_count >= this → Rule 2 "
                        "HIGH (per-donor cap = primary_count); "
                        "primary_count <  this → Rule 2 LOW (total cap = "
                        "donor_borrow_cap). Rule 5 fires when total after "
                        "Rule 2 LOW is still below this. Default 20 "
                        "(borrowed from downstream RCTD's CELL_MIN_INSTANCE).")
    p.add_argument("--rule1-threshold", type=int, default=None,
                   help="DEPRECATED alias for --donor-borrow-cap; kept for "
                        "backwards compat with pre-#21 configs.")
    p.add_argument("--primary-only-celltype", action="append", default=None,
                   dest="primary_only_celltypes",
                   help="Celltype always sourced exclusively from the primary. "
                        "Pass repeatedly. Default: tumor, liver.")
    p.add_argument("--random-seed", type=int, default=None,
                   help="Seed for the balanced/borrow samplers (default 42). "
                        "Also accepted as --random-state for backwards compat.")
    p.add_argument("--random-state", type=int, default=None,
                   help="DEPRECATED alias for --random-seed.")

    # export_mtx
    p.add_argument("--export-layer", default=None,
                   help="AnnData layer to export as mtx (default: 'counts', with "
                        "fallback to 'raw_count' then .X).")

    # rctd_reference_build
    p.add_argument("--rscript-bin", default=None,
                   help="Path to Rscript. Default: 'Rscript' (PATH lookup).")
    p.add_argument("--min-umi", type=int, default=None,
                   help="spacexr::Reference min_UMI (default: 10).")
    p.add_argument("--require-int", type=_str2bool, default=None,
                   help="spacexr::Reference require_int (default: true).")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="ref-build",
        description="scRNA reference-dataset builder for RCTD (spacexr) — "
                    "primary + donor pool → per-celltype adaptive migration → "
                    "10X-style mtx bundle → spacexr Reference .rds.",
    )
    parser.add_argument("--version", action="version",
                        version=f"ref-build {__version__}")
    sub = parser.add_subparsers(dest="cmd", required=True)

    run_p = sub.add_parser(
        "run",
        help="Run the five-stage reference-build pipeline.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    _add_run_args(run_p)
    return parser


def _resolve_config(args: argparse.Namespace) -> dict:
    """Load default + optional user YAML, then apply CLI overrides."""
    cfg = load_default()
    if args.config is not None:
        cfg = deep_update(cfg, load_yaml(args.config))

    overrides: dict = {}
    if args.sample_id is not None:
        overrides["sample_id"] = args.sample_id
    # run_id: resolve precedence HERE so downstream sees a concrete id.
    overrides["run_id"] = _resolve_run_id(getattr(args, "run_id", None))
    # flex-in-place: record the flex path verbatim; downstream reads
    # from it directly. When --primary-h5ad is not passed, use the
    # flex path as the primary input as well.
    if args.flex_h5ad is not None:
        overrides["flex_h5ad_path"] = str(args.flex_h5ad)
        if args.primary_h5ad is None:
            overrides["primary_h5ad"] = str(args.flex_h5ad)
    if args.primary_h5ad is not None:
        overrides["primary_h5ad"] = str(args.primary_h5ad)
    if args.donor_h5ads is not None:
        overrides["donor_h5ads"] = [str(p) for p in args.donor_h5ads]
    if args.fallback_donor_h5ads is not None:
        overrides["fallback_donor_h5ads"] = [str(p) for p in args.fallback_donor_h5ads]
    if args.celltype_marker_json is not None:
        overrides["celltype_marker_json"] = str(args.celltype_marker_json)
    if args.tumor_type is not None:
        overrides["tumor_type"] = args.tumor_type
    if args.output_root is not None:
        overrides["output_root"] = str(args.output_root)
    if args.force_rerun:
        overrides["force_rerun"] = True

    lpd_over = {}
    if args.celltype_col is not None:
        lpd_over["celltype_col"] = args.celltype_col
    if lpd_over:
        overrides["load_primary_and_donors"] = lpd_over

    cen_over = {}
    if args.celltype_col is not None:
        cen_over["celltype_col"] = args.celltype_col
    if args.donor_borrow_cap is not None:
        cen_over["donor_borrow_cap"] = args.donor_borrow_cap
    elif args.rule1_threshold is not None:
        # Backwards-compat alias — rule1_threshold used to be the Rule-1
        # gate; under the hybrid rule set it collapses into donor_borrow_cap.
        cen_over["donor_borrow_cap"] = args.rule1_threshold
    if args.cell_min_instance is not None:
        cen_over["cell_min_instance"] = args.cell_min_instance
    if args.celltype_target_list is not None:
        cen_over["celltype_target_list"] = str(args.celltype_target_list)
    if args.celltype_target_key is not None:
        cen_over["celltype_target_key"] = args.celltype_target_key
    if args.primary_only_celltypes is not None:
        cen_over["primary_only_celltypes"] = list(args.primary_only_celltypes)
    if args.random_seed is not None:
        cen_over["random_seed"] = args.random_seed
    elif args.random_state is not None:
        cen_over["random_seed"] = args.random_state
    if cen_over:
        overrides["census"] = cen_over

    asm_over = {}
    if args.celltype_col is not None:
        asm_over["celltype_col"] = args.celltype_col
    if asm_over:
        overrides["assemble"] = asm_over

    ex_over = {}
    if args.export_layer is not None:
        ex_over["layer"] = args.export_layer
    if ex_over:
        overrides["export_mtx"] = ex_over

    rr_over = {}
    if args.rscript_bin is not None:
        rr_over["rscript_bin"] = args.rscript_bin
    if args.min_umi is not None:
        rr_over["min_UMI"] = args.min_umi
    if args.require_int is not None:
        rr_over["require_int"] = args.require_int
    if rr_over:
        overrides["rctd_reference_build"] = rr_over

    cfg = deep_update(cfg, overrides)
    validate(cfg)
    return cfg


def main(argv: list[str] | None = None) -> int:
    import sys

    args = build_parser().parse_args(argv)
    if args.cmd == "run":
        from ref_build.pipeline import run
        cfg = _resolve_config(args)
        return run(cfg, stages=list(args.stages), argv=sys.argv)
    raise SystemExit(f"unknown command: {args.cmd}")


if __name__ == "__main__":
    raise SystemExit(main())
