#!/usr/bin/env python
"""Hotfix — fold RCTD results + per-step filter provenance onto
xenium_ranger.h5ad on an EXISTING finished v0.3-era run, without
re-running the pipeline.

Standalone sibling of the v0.2.9 celltype_writeback augmentation
(settylab/xenium-preprocess-pipeline branch v0.2.9,
settylab/TracyY123-nexus#40). Both call paths share one R-side
extractor (``packages/rctd-split/src/rctd_split/r/rctd_results_to_csv.R``)
+ one Python-side joiner (``_internal.rctd_results.extract_rctd_results_df``)
so the obs schema written here is byte-identical to what the pipeline
will write on fresh xenium_ranger-source runs.

Reads
-----
    <run-dir>/rctd/<S>_rctd_results.rds           (spacexr RCTD + SPLIT)
    <run-dir>/spatial_adata/<S>_xenium_ranger_purified.h5ad
                                                  (survived-all filter)
    <run-dir>/spatial_adata/<S>_xenium_ranger.h5ad (target: obs merged)

Writes (atomic in-place overwrite on xenium_ranger.h5ad)
--------------------------------------------------------
    obs[rctd_spot_class, rctd_first_type, rctd_second_type,
        rctd_first_class, rctd_second_class, rctd_min_score,
        rctd_singlet_score, rctd_conv_all, rctd_conv_doublet,
        rctd_purification_status, rctd_w1_larger_w2, rctd_same_class]
    obs[in_rctd_input]  (bool)
    obs[in_purified]    (bool)

Columns already present (celltype, celltype_source_distance) are NOT
touched. Columns being written are overwritten in place; re-running
the hotfix is idempotent.

Usage
-----
    module load R/4.4.0-gfbf-2023b           # or Tracy's conda-env R
    python scripts/hotfix_rctd_writeback_to_xenium_ranger.py \\
        --run-dir /path/to/<S>_<run_id>/ \\
        --sample-id <S> \\
        [--purified-h5ad-stem xenium_ranger_purified] \\
        [--rctd-results-rds ./rctd/<S>_rctd_results.rds] \\
        [--rscript-bin Rscript] \\
        [--r-lib-paths /home/<user>/R/x86_64-pc-linux-gnu-library/4.4]

Example (Tracy's MH8_2_BR_2453 v0.3 run):

    module load R/4.4.0-gfbf-2023b
    PYTHONPATH=/fh/fast/setty_m/user/ryang/nexus/nexus-code/work/\\
xenium-preprocess-pipeline-tracy-v0.2.7/packages/rctd-split/src \\
    /home/ryang/micromamba/envs/xenium_main/bin/python \\
        scripts/hotfix_rctd_writeback_to_xenium_ranger.py \\
        --run-dir /fh/fast/setty_m/user/ryang/xenium_preprocess_runs/\\
MH8_2_BR_2453_MH8_2_BR_2453_xenium_ranger_split_v0.3 \\
        --sample-id MH8_2_BR_2453 \\
        --purified-h5ad-stem xenium_ranger_purified \\
        --r-lib-paths /home/ryang/R/x86_64-pc-linux-gnu-library/4.4

Runtime: ~2-3 minutes (dominated by reading the 626 MB purified h5ad
and the 95 MB xenium_ranger h5ad + the 159 MB RCTD .rds via Rscript).
Idempotent: safe to re-run.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path


def _main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run-dir", required=True, type=Path,
                   help="Pipeline run directory (parent of spatial_adata/, "
                        "rctd/, split_prep/, summary/).")
    p.add_argument("--sample-id", required=True,
                   help="Sample ID — used to derive h5ad filenames.")
    p.add_argument("--purified-h5ad-stem", default="xenium_ranger_purified",
                   help="Stem of the purified h5ad (default: "
                        "'xenium_ranger_purified' — matches v0.2.8's "
                        "--purified-h5ad-stem override for xr-source runs).")
    p.add_argument("--rctd-results-rds", type=Path, default=None,
                   help="Explicit RCTD results .rds path. Defaults to "
                        "<run-dir>/rctd/<sample-id>_rctd_results.rds.")
    p.add_argument("--xenium-ranger-h5ad", type=Path, default=None,
                   help="Explicit target h5ad. Defaults to "
                        "<run-dir>/spatial_adata/<sample-id>_xenium_ranger.h5ad.")
    p.add_argument("--rscript-bin", default="Rscript",
                   help="Rscript binary (default 'Rscript'; load the R "
                        "module first, or pass an absolute path).")
    p.add_argument("--r-lib-paths", default=None,
                   help="Colon-separated R lib paths, prepended to "
                        "R_LIBS_USER for the Rscript subprocess.")
    p.add_argument("--rctd-obs-prefix", default="rctd_",
                   help="Prefix for the per-cell RCTD obs cols "
                        "(default 'rctd_' -> rctd_spot_class, rctd_first_type, ...).")
    args = p.parse_args(argv)

    run_dir = args.run_dir
    sample_id = args.sample_id
    rctd_rds = args.rctd_results_rds or (
        run_dir / "rctd" / f"{sample_id}_rctd_results.rds"
    )
    xr_h5ad = args.xenium_ranger_h5ad or (
        run_dir / "spatial_adata" / f"{sample_id}_xenium_ranger.h5ad"
    )
    pur_h5ad = run_dir / "spatial_adata" / (
        f"{sample_id}_{args.purified_h5ad_stem}.h5ad"
    )

    # --- preflight: all 3 inputs exist, writable target ---
    for p_in in (rctd_rds, xr_h5ad, pur_h5ad):
        if not p_in.exists():
            print(f"[hotfix] FATAL: input missing: {p_in}", file=sys.stderr)
            sys.exit(2)
    print(f"[hotfix] inputs OK")
    print(f"[hotfix]   rctd_results_rds : {rctd_rds}")
    print(f"[hotfix]   xenium_ranger.h5ad: {xr_h5ad}")
    print(f"[hotfix]   purified.h5ad    : {pur_h5ad}")

    # --- Reuse the pipeline's extract + fold helper (same code path
    # the v0.2.9 celltype_writeback stage uses). Reading the h5ad in
    # memory and writing back via atomic replace matches
    # `_internal.layout.atomic_write_h5ad`.
    import anndata as ad

    from rctd_split._internal.layout import atomic_write_h5ad
    from rctd_split.stages.celltype_writeback import _fold_rctd_results_onto_query

    r_lib_paths = args.r_lib_paths.split(":") if args.r_lib_paths else None

    print(f"[hotfix] reading target {xr_h5ad}")
    query = ad.read_h5ad(xr_h5ad)
    n_before = len(query.obs.columns)

    _fold_rctd_results_onto_query(
        query=query,
        reference_h5ad=pur_h5ad,
        rctd_results_rds=rctd_rds,
        rscript_bin=args.rscript_bin,
        r_lib_paths=r_lib_paths,
        rctd_obs_prefix=args.rctd_obs_prefix,
    )

    n_added = len(query.obs.columns) - n_before
    print(f"[hotfix] added {n_added} obs columns -> {sorted(query.obs.columns)}")

    print(f"[hotfix] atomic-writing merged h5ad -> {xr_h5ad}")
    atomic_write_h5ad(query, xr_h5ad, compression="gzip")
    print(f"[hotfix] done.")


if __name__ == "__main__":
    _main()
