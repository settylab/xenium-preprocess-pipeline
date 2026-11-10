#!/usr/bin/env python
"""Hotfix v2 — v0.2.9 Part 1 (settylab/TracyY123-nexus#40 c6755625444).

Fold RCTD results + the Part-1-corrected filter provenance onto an
existing finished xenium_ranger-source run WITHOUT re-running the
pipeline. Supersedes
`scripts/hotfix_rctd_writeback_to_xenium_ranger.py` (v1):

  * Corrects `in_rctd_input` semantics — now sourced from the
    split_prep `*_barcodes.tsv.gz` (cells SUBMITTED to RCTD), not from
    the rctd_results.rds (cells with a RESULT). Also writes a NEW
    `has_rctd_result` column = what v1 called `in_rctd_input`.
  * Writes `rctd_purification_status` sourced from the purified h5ad's
    `obs['purification_status']` (populated by SPLIT::purify,
    preserved through mtx_to_h5ad + postprocess).
  * Backs up the target h5ad BEFORE overwriting (atomic sibling copy).

Expected obs columns written (14 total, superseding v1's 13):

  * 11 `rctd_*` payload cols (same as v1) — `rctd_spot_class`,
    `rctd_first_type`, `rctd_second_type`, `rctd_first_class`,
    `rctd_second_class`, `rctd_min_score`, `rctd_singlet_score`,
    `rctd_conv_all`, `rctd_conv_doublet`, `rctd_w1_larger_w2`,
    `rctd_same_class`.
  * `rctd_purification_status` (NEW) — "" for cells not in purified,
    the SPLIT-assigned value otherwise.
  * `in_rctd_input` (semantic CORRECTED) — cells in the split_prep
    barcodes set (= cells submitted to RCTD).
  * `has_rctd_result` (NEW) — cells with a row in rctd_results_df.
  * `in_purified` (same as v1) — cells in the purified h5ad.

On MH8_2_BR_2453 v0.3:
  n_obs                   = 322,690
  in_rctd_input == True   = 299,585   (from split_prep barcodes)
  has_rctd_result == True = 255,661   (from rctd_results.rds)
  in_purified == True     = 186,080   (from purified.h5ad)

Usage::

    module load fhR/4.4.1-foss-2023b          # supplies Matrix + Seurat
    eval "$(/home/<user>/.local/bin/micromamba shell hook --shell bash)"
    micromamba activate /home/<user>/micromamba/envs/xenium_main
    unset PYTHONPATH

    PYTHONPATH=<staging>/packages/rctd-split/src \\
    python scripts/hotfix_rctd_writeback_v2.py \\
        --run-dir /path/to/<S>_<run_id>/ \\
        --sample-id <S> \\
        [--purified-h5ad-stem xenium_ranger_purified] \\
        [--rctd-input-barcodes split_prep/<S>_<suffix>_barcodes.tsv.gz] \\
        [--backup-suffix .pre-v0.2.9-part1.h5ad.bak] \\
        [--rscript-bin Rscript] \\
        [--r-lib-paths /home/<user>/R/x86_64-pc-linux-gnu-library/4.4]

Idempotent: re-running overwrites columns in place. Safe to re-run
against an h5ad already hotfixed by v1 (the stale `in_rctd_input` from
v1 is overwritten with the correct semantics; `has_rctd_result` is
added fresh).
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path


def _main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run-dir", required=True, type=Path,
                   help="Pipeline run directory.")
    p.add_argument("--sample-id", required=True)
    p.add_argument("--purified-h5ad-stem", default="xenium_ranger_purified")
    p.add_argument("--rctd-results-rds", type=Path, default=None,
                   help="Default: <run-dir>/rctd/<S>_rctd_results.rds.")
    p.add_argument("--xenium-ranger-h5ad", type=Path, default=None,
                   help="Default: <run-dir>/spatial_adata/<S>_xenium_ranger.h5ad.")
    p.add_argument("--rctd-input-barcodes", type=Path, default=None,
                   help="Default: autodiscover under "
                        "<run-dir>/split_prep/<S>*_barcodes.tsv.gz "
                        "(unique match required).")
    p.add_argument("--backup-suffix", default=".pre-v0.2.9-part1.h5ad.bak",
                   help="Suffix appended to the target h5ad for the "
                        "pre-hotfix backup copy. Set to the empty string "
                        "to skip the backup.")
    p.add_argument("--rscript-bin", default="Rscript")
    p.add_argument("--r-lib-paths", default=None)
    p.add_argument("--rctd-obs-prefix", default="rctd_")
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

    if args.rctd_input_barcodes is not None:
        rctd_input_barcodes = args.rctd_input_barcodes
    else:
        sp_dir = run_dir / "split_prep"
        matches = sorted(sp_dir.glob(f"{sample_id}*_barcodes.tsv.gz"))
        if len(matches) == 1:
            rctd_input_barcodes = matches[0]
        elif len(matches) > 1:
            print(f"[hotfix] FATAL: multiple barcodes.tsv.gz under "
                  f"{sp_dir}: {[p.name for p in matches]}. Pass "
                  f"--rctd-input-barcodes explicitly.", file=sys.stderr)
            sys.exit(4)
        else:
            print(f"[hotfix] FATAL: no barcodes.tsv.gz under {sp_dir}. "
                  f"Can't source in_rctd_input. Pass --rctd-input-barcodes "
                  f"explicitly or check that split_prep ran for this run.",
                  file=sys.stderr)
            sys.exit(4)

    for p_in in (rctd_rds, xr_h5ad, pur_h5ad, rctd_input_barcodes):
        if not p_in.exists():
            print(f"[hotfix] FATAL: input missing: {p_in}", file=sys.stderr)
            sys.exit(2)
    print("[hotfix] inputs OK")
    print(f"[hotfix]   rctd_results_rds : {rctd_rds}")
    print(f"[hotfix]   xenium_ranger.h5ad: {xr_h5ad}")
    print(f"[hotfix]   purified.h5ad    : {pur_h5ad}")
    print(f"[hotfix]   rctd_input_barcodes: {rctd_input_barcodes}")

    # Backup target BEFORE touching it. Skip if suffix is "".
    if args.backup_suffix:
        backup = xr_h5ad.with_name(xr_h5ad.name + args.backup_suffix)
        if backup.exists():
            print(f"[hotfix] backup already exists (not overwriting): {backup}")
        else:
            print(f"[hotfix] backing up target -> {backup}")
            shutil.copy2(xr_h5ad, backup)

    # Reuse the pipeline helper — single source of truth for the write
    # schema. Changes to column semantics / fill behavior land in both
    # pipeline + hotfix automatically.
    import anndata as ad

    from rctd_split._internal.layout import atomic_write_h5ad
    from rctd_split.stages.celltype_writeback import _fold_rctd_results_onto_query

    r_lib_paths = args.r_lib_paths.split(":") if args.r_lib_paths else None

    print(f"[hotfix] reading target {xr_h5ad}")
    query = ad.read_h5ad(xr_h5ad)
    cols_before = set(query.obs.columns)

    _fold_rctd_results_onto_query(
        query=query,
        reference_h5ad=pur_h5ad,
        rctd_results_rds=rctd_rds,
        rscript_bin=args.rscript_bin,
        r_lib_paths=r_lib_paths,
        rctd_obs_prefix=args.rctd_obs_prefix,
        rctd_input_barcodes_path=rctd_input_barcodes,
    )

    cols_after = set(query.obs.columns)
    added = sorted(cols_after - cols_before)
    print(f"[hotfix] added {len(added)} obs columns: {added}")

    print(f"[hotfix] atomic-writing merged h5ad -> {xr_h5ad}")
    atomic_write_h5ad(query, xr_h5ad, compression="gzip")
    print("[hotfix] done.")


if __name__ == "__main__":
    _main()
