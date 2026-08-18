"""Command-line entrypoint for rctd-split.

Exposes `rctd-split run …` via the `[project.scripts]` entry point in
pyproject.toml. Design mirrors xenium_preprocess.cli — argparse builds
the surface, the config machinery in rctd_split.config layers default
YAML + user YAML + CLI overrides, and the pipeline stage-loop in
rctd_split.pipeline runs the requested subset.

`--run-id` precedence (bound HERE so the pipeline module sees an
already-resolved id):
    `--run-id` > `$SLURM_JOB_ID` > `YYYYMMDD_HHMMSS` timestamp fallback.
"""
from __future__ import annotations

import argparse
import os
import time
from pathlib import Path

from rctd_split import __version__
from rctd_split.config import (
    DEFAULT_STAGES,
    VALID_STAGES,
    deep_update,
    load_default,
    load_yaml,
    validate,
)


def _str2bool(v):
    if isinstance(v, bool):
        return v
    s = str(v).strip().lower()
    if s in {"true", "1", "yes", "y", "on"}:
        return True
    if s in {"false", "0", "no", "n", "off"}:
        return False
    raise argparse.ArgumentTypeError(f"expected boolean, got {v!r}")


def _parse_extra_reports(raw: list[str]) -> list[dict]:
    """Parse each ``--extra-report`` CLI value into a validated dict.

    Accepts the shorthand ``PATH,NAME`` (single-comma split); anything
    after the first comma is treated as part of the name so display
    names can themselves contain commas. Both fields are required —
    fail-loud on either missing.
    """
    out: list[dict] = []
    for token in raw:
        if "," not in token:
            raise SystemExit(
                f"[cli] --extra-report expected 'PATH,NAME'; got {token!r}"
            )
        path, name = token.split(",", 1)
        path, name = path.strip(), name.strip()
        if not path or not name:
            raise SystemExit(
                f"[cli] --extra-report has empty PATH or NAME: {token!r}"
            )
        out.append({"path": path, "name": name})
    return out


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


def _add_run_args(p: argparse.ArgumentParser) -> None:
    """Attach every `run` subcommand flag to `p`."""
    p.add_argument("--config", type=Path, default=None,
                   help="User config YAML (overrides config/default.yaml entries).")
    p.add_argument("--stages", nargs="+", choices=VALID_STAGES,
                   default=list(DEFAULT_STAGES),
                   help="Which stages to run.")
    p.add_argument("--run-id", default=None,
                   help="Explicit run identifier. Precedence: --run-id > "
                        "$SLURM_JOB_ID > YYYYMMDD_HHMMSS. Names the run "
                        "folder <output_root>/<sample>/<sample>_<run_id>/.")

    # I/O
    p.add_argument("--sample-id",
                   help="Sample identifier (e.g. SAMPLE1).")
    p.add_argument("--test-object", type=Path,
                   help="Path to the spatial test object RDS "
                        "(step 1 rctd_prep output). Optional — when "
                        "omitted, resolves to "
                        "<output_root>/<sample>/<sample>_<run_id>/rctd/"
                        "<sample>_test_object.rds (the standard "
                        "workflow layout).")
    p.add_argument("--reference-rds", type=Path,
                   help="Path to the scRNA reference RDS "
                        "(step 3 rctd_reference_build output — a "
                        "spacexr::Reference object). Optional — when "
                        "omitted, resolves to "
                        "<output_root>/<sample>/<sample>_<run_id>/rctd/"
                        "<sample>_reference.rds (the standard workflow "
                        "layout).")
    p.add_argument("--output-root", type=Path,
                   help="Root output directory.")
    p.add_argument("--force-rerun", action="store_true",
                   help="Re-run all stages even if sentinel outputs exist.")
    p.add_argument("--keep-intermediate", action="store_true", default=None,
                   help="Keep <run_dir>/intermediate/ after the terminal "
                        "stage (celltype_writeback) finishes. Default: drop "
                        "it (matches ref-build's post-pipeline cleanup). "
                        "Pass this to retain the intermediate R RDS + mtx "
                        "bundles + unpurified h5ad for debugging / QC.")

    # Cross-stage R config
    p.add_argument("--r-lib-paths", nargs="+", default=None,
                   help="R library paths prepended to R_LIBS_USER for every "
                        "R stage. Typically the renv library dir(s) that hold "
                        "spacexr + SPLIT + Seurat.")

    # rctd_run
    p.add_argument("--rscript-bin", default=None,
                   help="Path to Rscript. Default: 'Rscript' (PATH lookup).")
    p.add_argument("--max-cores", type=int, default=None,
                   help="RCTD max_cores parameter (default: 4, matches Rmd).")
    p.add_argument("--umi-min", type=int, default=None,
                   help="RCTD UMI_min parameter (default: 10, matches Rmd).")
    p.add_argument("--counts-min", type=int, default=None,
                   help="RCTD counts_MIN parameter (default: 10, matches Rmd).")
    p.add_argument("--umi-min-sigma", type=int, default=None,
                   help="RCTD UMI_min_sigma parameter (default: 100, matches Rmd).")
    p.add_argument("--cell-min-instance", type=int, default=None,
                   help="RCTD CELL_MIN_INSTANCE parameter (default: 20 "
                        "— Pipeline default; Rmd used 25).")
    p.add_argument("--doublet-mode", default=None,
                   help="run.RCTD doublet_mode (default: 'doublet', matches Rmd).")
    p.add_argument("--assay-name", default=None,
                   help="Seurat assay name on the spatial test object (default: 'Proseg').")

    # split_purify
    p.add_argument("--do-purify-singlets", type=_str2bool, default=None,
                   help="SPLIT::purify DO_purify_singlets parameter "
                        "(default: true, matches Rmd).")

    # filter_status
    p.add_argument("--filter-status-write-inplace", type=_str2bool, default=None,
                   help="filter_status: also augment unpurified.h5ad's .obs "
                        "with the three filter-status columns "
                        "(default: true). Sidecar CSV is written regardless.")

    # postprocess
    p.add_argument("--postprocess-min-counts", type=int, default=None,
                   help="postprocess: sc.pp.filter_cells min_counts threshold "
                        "(default: 50, matches the user's notebook).")
    p.add_argument("--postprocess-minprop", type=float, default=None,
                   help="postprocess: MINPROP in log(clip(X_row_normalized, "
                        "MINPROP, 1)) (default: 1e-3, matches notebook).")
    p.add_argument("--postprocess-positive-x", type=_str2bool, default=None,
                   help="postprocess: shift X_clipped by -log(minprop) so all "
                        "values are ≥ 0 (default: true).")
    p.add_argument("--postprocess-n-neighbors", type=int, default=None,
                   help="postprocess: NNEIGHBORS for KNN + UMAP + Leiden "
                        "(default: 15, matches notebook).")
    p.add_argument("--postprocess-leiden-resolutions", nargs="+", type=float,
                   default=None,
                   help="postprocess: leiden resolution(s) to run "
                        "(default: [0.5, 0.7], matches notebook).")
    p.add_argument("--postprocess-process-unpurified-layer", type=_str2bool,
                   default=None,
                   help="postprocess: also compute PCA/UMAP on "
                        "layers['unpurified_counts'] (default: true).")

    # writeback_to_step1_raw
    p.add_argument("--step1-raw-h5ad", type=Path, default=None,
                   help="writeback_to_step1_raw: path to step-1's "
                        "<S>_proseg_raw.h5ad. Default resolves to "
                        "<output_root>/<sample>/<sample>_<run_id>/"
                        "spatial_adata/<sample>_proseg_raw.h5ad.")

    # celltype_writeback
    p.add_argument("--step1-xenium-ranger-h5ad", type=Path, default=None,
                   help="celltype_writeback: path to step-1's "
                        "<S>_xenium_ranger.h5ad. Default resolves to "
                        "<output_root>/<sample>/<sample>_<run_id>/"
                        "spatial_adata/<sample>_xenium_ranger.h5ad.")

    # qc_report
    p.add_argument("--qc-purification-status-column", default=None,
                   help="qc_report: obs column on proseg_purified "
                        "carrying the SPLIT purification status "
                        "categorical used to color UMAP #1 "
                        "(default: purification_status).")
    p.add_argument("--qc-raw-layer", default=None,
                   help="qc_report: layer on proseg_raw.h5ad whose "
                        "matrix drives raw-side QC metrics + histogram "
                        "(default: maxpost_counts; fail-loud if absent).")
    p.add_argument("--extra-report", action="append", default=None,
                   metavar="PATH,NAME",
                   help="qc_report: append a link to an external HTML "
                        "report at the bottom of the summary_report.html. "
                        "Format is 'PATH,NAME' (single comma splits into "
                        "the HTML path and the display name). May be "
                        "repeated for multiple reports; equivalent to "
                        "`qc_report.extra_reports: [{path, name}, ...]` "
                        "in the config YAML.")


def _add_strip_purified_obs_args(p: argparse.ArgumentParser) -> None:
    """Attach `strip-purified-obs` flags to `p`."""
    p.add_argument("--purified-h5ad", type=Path, required=True,
                   help="Path to an existing proseg_purified.h5ad. The "
                        "file is rewritten in place after a timestamped "
                        ".bak copy is made next to it.")
    p.add_argument("--no-backup", action="store_true",
                   help="Skip writing the .bak copy. Not recommended.")
    p.add_argument("--dry-run", action="store_true",
                   help="Print columns that would be dropped; do not "
                        "modify or back up the input h5ad.")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="rctd-split",
        description="RCTD + SPLIT typing (Stage D) with mtx exports + h5ad conversion.",
    )
    parser.add_argument("--version", action="version",
                        version=f"rctd-split {__version__}")
    sub = parser.add_subparsers(dest="cmd", required=True)

    run_p = sub.add_parser(
        "run",
        help="Run the RCTD + SPLIT typing pipeline.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    _add_run_args(run_p)

    strip_p = sub.add_parser(
        "strip-purified-obs",
        help="Retroactively drop non-SPLIT-native obs columns from an "
             "existing proseg_purified.h5ad (idempotent; backs up first).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    _add_strip_purified_obs_args(strip_p)

    return parser


def _resolve_config(args: argparse.Namespace) -> dict:
    """Load default + optional user YAML, then apply CLI overrides."""
    cfg = load_default()
    if args.config is not None:
        cfg = deep_update(cfg, load_yaml(args.config))

    overrides: dict = {}
    if args.sample_id is not None:
        overrides["sample_id"] = args.sample_id
    if args.test_object is not None:
        overrides["test_object"] = str(args.test_object)
    if args.reference_rds is not None:
        overrides["reference_rds"] = str(args.reference_rds)
    if args.output_root is not None:
        overrides["output_root"] = str(args.output_root)
    if args.force_rerun:
        overrides["force_rerun"] = True
    if args.keep_intermediate:
        overrides["keep_intermediate"] = True
    if args.r_lib_paths is not None:
        overrides["r_lib_paths"] = [str(p) for p in args.r_lib_paths]

    # run_id: resolve precedence HERE so downstream sees a concrete id.
    overrides["run_id"] = _resolve_run_id(getattr(args, "run_id", None))

    rr_over: dict = {}
    if args.rscript_bin is not None:
        rr_over["rscript_bin"] = args.rscript_bin
    if args.max_cores is not None:
        rr_over["max_cores"] = args.max_cores
    if args.umi_min is not None:
        rr_over["UMI_min"] = args.umi_min
    if args.counts_min is not None:
        rr_over["counts_MIN"] = args.counts_min
    if args.umi_min_sigma is not None:
        rr_over["UMI_min_sigma"] = args.umi_min_sigma
    if args.cell_min_instance is not None:
        rr_over["CELL_MIN_INSTANCE"] = args.cell_min_instance
    if args.doublet_mode is not None:
        rr_over["doublet_mode"] = args.doublet_mode
    if args.assay_name is not None:
        rr_over["assay_name"] = args.assay_name
    if rr_over:
        overrides["rctd_run"] = rr_over

    sp_over: dict = {}
    if args.rscript_bin is not None:
        sp_over["rscript_bin"] = args.rscript_bin
    if args.assay_name is not None:
        sp_over["assay_name"] = args.assay_name
    if args.do_purify_singlets is not None:
        sp_over["DO_purify_singlets"] = args.do_purify_singlets
    if sp_over:
        overrides["split_purify"] = sp_over

    em_over: dict = {}
    if args.rscript_bin is not None:
        em_over["rscript_bin"] = args.rscript_bin
    if args.assay_name is not None:
        em_over["assay_name"] = args.assay_name
    if em_over:
        overrides["export_mtx"] = em_over

    fs_over: dict = {}
    if args.filter_status_write_inplace is not None:
        fs_over["write_inplace"] = args.filter_status_write_inplace
    if fs_over:
        overrides["filter_status"] = fs_over

    pp_over: dict = {}
    pp_qc_over: dict = {}
    pp_norm_over: dict = {}
    pp_nn_over: dict = {}
    pp_ld_over: dict = {}
    if args.postprocess_min_counts is not None:
        pp_qc_over["min_counts"] = args.postprocess_min_counts
    if args.postprocess_minprop is not None:
        pp_norm_over["minprop"] = args.postprocess_minprop
    if args.postprocess_positive_x is not None:
        pp_norm_over["positive_x"] = args.postprocess_positive_x
    if args.postprocess_n_neighbors is not None:
        pp_nn_over["n_neighbors"] = args.postprocess_n_neighbors
    if args.postprocess_leiden_resolutions is not None:
        pp_ld_over["resolutions"] = list(args.postprocess_leiden_resolutions)
    if args.postprocess_process_unpurified_layer is not None:
        pp_over["process_unpurified_layer"] = args.postprocess_process_unpurified_layer
    if pp_qc_over:
        pp_over["qc"] = pp_qc_over
    if pp_norm_over:
        pp_over["normalize"] = pp_norm_over
    if pp_nn_over:
        pp_over["neighbors"] = pp_nn_over
    if pp_ld_over:
        pp_over["leiden"] = pp_ld_over
    if pp_over:
        overrides["postprocess"] = pp_over

    wb_over: dict = {}
    if args.step1_raw_h5ad is not None:
        wb_over["raw_h5ad"] = str(args.step1_raw_h5ad)
    if wb_over:
        overrides["writeback_to_step1_raw"] = wb_over

    ct_over: dict = {}
    if args.step1_xenium_ranger_h5ad is not None:
        ct_over["xenium_ranger_h5ad"] = str(args.step1_xenium_ranger_h5ad)
    if ct_over:
        overrides["celltype_writeback"] = ct_over

    qc_over: dict = {}
    if args.qc_purification_status_column is not None:
        qc_over["purification_status_column"] = args.qc_purification_status_column
    if args.qc_raw_layer is not None:
        qc_over["raw_layer"] = args.qc_raw_layer
    if getattr(args, "extra_report", None):
        qc_over["extra_reports"] = _parse_extra_reports(args.extra_report)
    if qc_over:
        overrides["qc_report"] = qc_over

    cfg = deep_update(cfg, overrides)
    validate(cfg)
    return cfg


def _run_strip_purified_obs(args: argparse.Namespace) -> int:
    """`strip-purified-obs` subcommand: drop non-SPLIT-native obs cols
    from an existing proseg_purified.h5ad, in place, with a
    timestamped .bak copy."""
    import shutil

    import anndata as ad

    from rctd_split._internal.layout import atomic_write_h5ad
    from rctd_split._internal.split_native import (
        filter_purified_obs,
        is_split_native_obs_col,
    )

    p = Path(args.purified_h5ad)
    if not p.exists():
        raise SystemExit(f"[strip-purified-obs] not found: {p}")

    print(f"[strip-purified-obs] reading {p}")
    adata = ad.read_h5ad(p)
    orig_cols = list(adata.obs.columns)

    to_drop = sorted(c for c in orig_cols if not is_split_native_obs_col(c))

    print(f"[strip-purified-obs]   obs cols before: {len(orig_cols)}")
    print(f"[strip-purified-obs]   kept:    {sorted(set(orig_cols) - set(to_drop))}")
    print(f"[strip-purified-obs]   dropped: {to_drop}")

    if args.dry_run:
        print("[strip-purified-obs] --dry-run: not modifying "
              f"{p} (would drop {len(to_drop)} col(s)).")
        return 0

    if not to_drop:
        print("[strip-purified-obs] nothing to drop — file already clean; "
              "skipping backup + rewrite.")
        return 0

    if not args.no_backup:
        stamp = time.strftime("%Y%m%d-%H%M%S")
        bak = p.with_name(p.name + f".bak-{stamp}")
        print(f"[strip-purified-obs] backup: {p} -> {bak}")
        shutil.copy2(p, bak)

    dropped = filter_purified_obs(adata)
    assert dropped == to_drop  # sanity: pre-computed list matches mutation
    atomic_write_h5ad(adata, p)
    print(f"[strip-purified-obs] wrote {p} "
          f"(obs cols now: {list(adata.obs.columns)})")
    return 0


def main(argv: list[str] | None = None) -> int:
    import sys

    args = build_parser().parse_args(argv)
    if args.cmd == "run":
        from rctd_split.pipeline import run
        cfg = _resolve_config(args)
        return run(cfg, stages=list(args.stages), argv=sys.argv)
    if args.cmd == "strip-purified-obs":
        return _run_strip_purified_obs(args)
    raise SystemExit(f"unknown command: {args.cmd}")


if __name__ == "__main__":
    raise SystemExit(main())
