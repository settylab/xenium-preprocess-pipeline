"""Command-line entrypoint for xenium-preprocess.

Exposes `xenium-preprocess run …` via the `[project.scripts]` entry
point in pyproject.toml. Design mirrors hexenium.cli — argparse builds
the surface, the config machinery in xenium_preprocess.config layers
default YAML + user YAML + CLI overrides, and the pipeline stage-loop
in xenium_preprocess.pipeline runs the requested subset.

`--run-id` precedence (bound HERE so the pipeline module sees a
already-resolved id):
    `--run-id` > `$SLURM_JOB_ID` > `YYYYMMDD_HHMMSS` timestamp fallback.
"""
from __future__ import annotations

import argparse
import os
import time
from pathlib import Path

from xenium_preprocess import __version__
from xenium_preprocess.config import (
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
                   help="Which stages to run. Default omits `preprocess` "
                        "(the PCA/UMAP/Leiden/celltype pass) and "
                        "`split_prep` (the persisted 10x-style mtx "
                        "bundle). `rctd_prep` is a default; its R script "
                        "consumes a split_prep bundle produced as a "
                        "transient prerequisite and deleted after the "
                        "RDS is written. Add `split_prep` explicitly to "
                        "keep the bundle on disk.")
    p.add_argument("--run-id", default=None,
                   help="Explicit run identifier. Precedence: --run-id > "
                        "$SLURM_JOB_ID > YYYYMMDD_HHMMSS. Names the run "
                        "folder <output_root>/<sample>/<sample>_<run_id>/.")

    # I/O
    p.add_argument("--sample-id",
                   help="Sample identifier (e.g. SAMPLE1).")
    p.add_argument("--proseg-dir", type=Path,
                   help="Directory holding the proseg output (count matrix + cell metadata).")
    p.add_argument("--output-root", type=Path,
                   help="Root output directory.")
    p.add_argument("--force-rerun", action="store_true",
                   help="Re-run all stages even if sentinel outputs exist.")

    # proseg_to_anndata
    p.add_argument("--count-matrix", type=Path, default=None,
                   help="Explicit path to the proseg expected-counts matrix "
                        "(.parquet or .csv). Overrides the glob in "
                        "config.proseg_to_anndata.count_matrix_glob.")
    p.add_argument("--maxpost-matrix", type=Path, default=None,
                   help="Explicit path to the proseg maxpost-counts matrix.")
    p.add_argument("--cell-metadata", type=Path, default=None,
                   help="Explicit path to the proseg cell metadata "
                        "(.parquet or .csv).")
    p.add_argument("--proseg-run-script", type=Path, default=None,
                   help="Optional audit trail: path to the proseg-run "
                        "shell script that produced the counts. When set, "
                        "copied verbatim to spatial_adata/provenance/"
                        "<script-name> alongside the h5ads. Omitted when "
                        "unset; no empty folder is created.")
    p.add_argument("--x-source",
                   choices=("maxpost_counts", "expected_counts"),
                   default=None,
                   help="Which proseg matrix `.X` mirrors in the raw h5ad. "
                        "Default (from config): 'maxpost_counts'.")

    # xenium_ranger_to_anndata
    p.add_argument("--xenium-ranger-dir", type=Path, default=None,
                   help="Directory holding the xenium-ranger bundle "
                        "(`cell_feature_matrix.h5` + `cells.csv.gz`). "
                        "Required by the `xenium_ranger_to_anndata` stage.")
    p.add_argument("--gex-only", type=_str2bool, default=None,
                   help="`gex_only` flag threaded into sc.read_10x_h5 for the "
                        "xenium-ranger matrix (default: true — 5001 features on the reference fixture).")

    # qc_filter
    p.add_argument("--qc-min-counts-cell", type=int, default=None,
                   help="`qc_filter.min_counts_cell` — cells with .X row sum "
                        ">= this get .obs['qc_filtered']=True. Default: 10.")

    # preprocess (kept for reversibility; NOT in DEFAULT_STAGES)
    p.add_argument("--dual-matrix-mode", type=_str2bool, default=None,
                   help="preprocess: run two independent passes (expected + maxpost).")
    p.add_argument("--leiden-resolution", type=float, default=None,
                   help="preprocess: Leiden clustering resolution (default: 0.4).")
    p.add_argument("--min-counts-cell", type=int, default=None,
                   help="preprocess: minimum total counts per cell in QC filter.")
    p.add_argument("--n-neighbors", type=int, default=None,
                   help="preprocess: neighbours for KNN + UMAP.")
    p.add_argument("--min-prop", type=float, default=None,
                   help="preprocess: clip floor for clipped-log normalization.")
    p.add_argument("--positive-x", type=_str2bool, default=None,
                   help="preprocess: shift X_clipped to positive values.")
    p.add_argument("--random-state", type=int, default=None,
                   help="preprocess: random seed for KNN + UMAP + Leiden.")
    p.add_argument("--global-non-tumor-json", type=Path, default=None,
                   help="preprocess: global-non-tumor marker JSON.")
    p.add_argument("--global-tumor-json", type=Path, default=None,
                   help="preprocess: global-tumor marker JSON.")
    p.add_argument("--tumor-type", default=None,
                   help="preprocess: tumor type key into global_tumor_json.")

    # split_prep
    p.add_argument("--split-layer", default=None,
                   help="AnnData layer to export as mtx.")
    p.add_argument("--split-name-suffix", default=None,
                   help="File-name suffix for the mtx/features/barcodes triple.")

    # rctd_prep
    p.add_argument("--rscript-bin", default=None,
                   help="Path to Rscript. Default: 'Rscript' (PATH lookup).")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="xenium-preprocess",
        description="Proseg → adata → SPLIT mtx/features/barcodes → RCTD test object.",
    )
    parser.add_argument("--version", action="version",
                        version=f"xenium-preprocess {__version__}")
    sub = parser.add_subparsers(dest="cmd", required=True)

    run_p = sub.add_parser(
        "run",
        help="Run the spatial-data preprocessing pipeline.",
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
    if args.proseg_dir is not None:
        overrides["proseg_dir"] = str(args.proseg_dir)
    if args.output_root is not None:
        overrides["output_root"] = str(args.output_root)
    if args.force_rerun:
        overrides["force_rerun"] = True

    # run_id: resolve precedence HERE so downstream sees a concrete id.
    overrides["run_id"] = _resolve_run_id(getattr(args, "run_id", None))

    pta_over = {}
    if args.count_matrix is not None:
        pta_over["count_matrix_path"] = str(args.count_matrix)
    if args.maxpost_matrix is not None:
        pta_over["maxpost_matrix_path"] = str(args.maxpost_matrix)
    if args.cell_metadata is not None:
        pta_over["cell_metadata_path"] = str(args.cell_metadata)
    if args.proseg_run_script is not None:
        pta_over["proseg_run_script"] = str(args.proseg_run_script)
    if args.x_source is not None:
        pta_over["x_source"] = args.x_source
    if pta_over:
        overrides["proseg_to_anndata"] = pta_over

    xr_over = {}
    if args.xenium_ranger_dir is not None:
        xr_over["xenium_ranger_dir"] = str(args.xenium_ranger_dir)
    if args.gex_only is not None:
        xr_over["gex_only"] = args.gex_only
    if xr_over:
        overrides["xenium_ranger_to_anndata"] = xr_over

    qc_over = {}
    if args.qc_min_counts_cell is not None:
        qc_over["min_counts_cell"] = args.qc_min_counts_cell
    if qc_over:
        overrides["qc_filter"] = qc_over

    pp_over = {}
    if args.dual_matrix_mode is not None:
        pp_over["dual_matrix_mode"] = args.dual_matrix_mode
    if args.leiden_resolution is not None:
        pp_over["leiden_resolution"] = args.leiden_resolution
    if args.min_counts_cell is not None:
        pp_over["min_counts_cell"] = args.min_counts_cell
    if args.n_neighbors is not None:
        pp_over["n_neighbors"] = args.n_neighbors
    if args.min_prop is not None:
        pp_over["min_prop"] = args.min_prop
    if args.positive_x is not None:
        pp_over["positive_x"] = args.positive_x
    if args.random_state is not None:
        pp_over["random_state"] = args.random_state
    if args.global_non_tumor_json is not None:
        pp_over["global_non_tumor_json"] = str(args.global_non_tumor_json)
    if args.global_tumor_json is not None:
        pp_over["global_tumor_json"] = str(args.global_tumor_json)
    if args.tumor_type is not None:
        pp_over["tumor_type"] = args.tumor_type
    if pp_over:
        overrides["preprocess"] = pp_over

    sp_over = {}
    if args.split_layer is not None:
        sp_over["layer"] = args.split_layer
    if args.split_name_suffix is not None:
        sp_over["name_suffix"] = args.split_name_suffix
    if sp_over:
        overrides["split_prep"] = sp_over

    rp_over = {}
    if args.rscript_bin is not None:
        rp_over["rscript_bin"] = args.rscript_bin
    if rp_over:
        overrides["rctd_prep"] = rp_over

    cfg = deep_update(cfg, overrides)
    validate(cfg)
    return cfg


def main(argv: list[str] | None = None) -> int:
    import sys

    args = build_parser().parse_args(argv)
    if args.cmd == "run":
        from xenium_preprocess.pipeline import run
        cfg = _resolve_config(args)
        return run(cfg, stages=list(args.stages), argv=sys.argv)
    raise SystemExit(f"unknown command: {args.cmd}")


if __name__ == "__main__":
    raise SystemExit(main())
