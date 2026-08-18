"""Stage-orchestration loop for xenium-preprocess.

Runs the requested subset of the pipeline stages, each guarded by its
own sentinel-existence resume check (nuke with `force_rerun`).

Post-2026-08-11 refactor: outputs live under
`<output_root>/<sample_id>/<sample_id>_<run_id>/{spatial_adata,rctd}/`,
and `config.yaml` at that run folder is written via the shared
`_merge_config` helper so ref-build / rctd-split (spawned as separate sbatch
jobs in the driver chain) don't clobber each other's config sections.
"""
from __future__ import annotations

import os
import socket
import sys
import time
from pathlib import Path

from xenium_preprocess._internal.layout import (
    resolved_config_path,
    run_dir,
    spatial_adata_path,
)
from xenium_preprocess._internal.logging import banner, log
from xenium_preprocess._internal.merge_config import merge_config


def log_invocation_banner(argv: list[str], stages: list[str], cfg: dict) -> None:
    """Print everything anyone would want to grep out of the log if the
    job ever silently does nothing again."""
    banner("xenium-preprocess starting")
    log(f"command:      {' '.join(argv)}")
    log(f"host:         {socket.gethostname()}")
    log(f"cwd:          {os.getcwd()}")
    log(f"pid:          {os.getpid()}")
    log(f"python:       {sys.executable}")
    log(f"python ver:   {sys.version.splitlines()[0]}")
    log(f"sample_id:    {cfg.get('sample_id')}")
    log(f"run_id:       {cfg.get('run_id')}")
    log(f"proseg_dir:   {cfg.get('proseg_dir')}")
    log(f"output_root:  {cfg.get('output_root')}")
    log(f"stages:       {stages}")
    log(f"force_rerun:  {cfg.get('force_rerun', False)}")
    for name in ("numpy", "yaml", "scipy", "scanpy", "anndata"):
        try:
            mod = __import__(name)
            ver = getattr(mod, "__version__", "no __version__ attr")
            log(f"  {name:14s}{ver}")
        except Exception as e:
            log(f"  {name:14s}NOT IMPORTABLE ({type(e).__name__})")
    log("=" * 64)


def run(cfg: dict, stages: list[str], argv: list[str]) -> int:
    """Execute the pipeline. Assumes `cfg` has already been validated."""
    sample_id = cfg["sample_id"]
    run_id = cfg["run_id"]
    proseg_dir = Path(cfg["proseg_dir"]).resolve()
    output_root = Path(cfg["output_root"]).resolve()
    force_rerun = bool(cfg.get("force_rerun", False))

    log_invocation_banner(argv, stages, cfg)

    if not proseg_dir.exists():
        raise SystemExit(f"proseg_dir not found: {proseg_dir}")

    # Materialise the run folder up front. `spatial_adata/` / `rctd/`
    # are created lazily by the stages that write into them.
    the_run_dir = run_dir(output_root, sample_id, run_id)
    the_run_dir.mkdir(parents=True, exist_ok=True)
    log(f"[pipeline] run folder: {the_run_dir}")

    # Merged resolved_config: write ONLY the xenium_preprocess key,
    # preserve any ref_build/rctd_split/driver sections a sibling pipeline
    # already wrote.
    snap = resolved_config_path(output_root, sample_id, run_id)
    xp_cfg = dict(cfg)
    merge_config(snap, step_key="xenium_preprocess", step_cfg=xp_cfg)
    log(f"[pipeline] merged resolved_config -> {snap}")

    n_stages = len(stages)
    pipeline_t0 = time.time()

    # Path handles used by multiple stages.
    raw_h5ad = spatial_adata_path(output_root, sample_id, run_id, "proseg_raw")
    xr_h5ad = spatial_adata_path(output_root, sample_id, run_id, "xenium_ranger")

    # --- Stage 1: proseg → raw AnnData ------------------------------
    if "proseg_to_anndata" in stages:
        idx = stages.index("proseg_to_anndata") + 1
        banner(f"stage {idx}/{n_stages}: proseg_to_anndata — starting")
        t0 = time.time()
        from xenium_preprocess.stages.proseg_to_anndata import run_proseg_to_anndata
        pta = cfg.get("proseg_to_anndata", {}) or {}
        run_proseg_to_anndata(
            sample_id=sample_id,
            run_id=run_id,
            proseg_dir=proseg_dir,
            output_root=output_root,
            count_matrix_path=(
                Path(pta["count_matrix_path"]) if pta.get("count_matrix_path") else None
            ),
            cell_metadata_path=(
                Path(pta["cell_metadata_path"]) if pta.get("cell_metadata_path") else None
            ),
            maxpost_matrix_path=(
                Path(pta["maxpost_matrix_path"]) if pta.get("maxpost_matrix_path") else None
            ),
            count_matrix_glob=pta.get("count_matrix_glob", "expected-counts*"),
            cell_metadata_glob=pta.get("cell_metadata_glob", "cell-metadata*"),
            maxpost_matrix_glob=pta.get("maxpost_matrix_glob", "maxpost_counts*"),
            centroid_x_col=pta.get("centroid_x_col", "centroid_x"),
            centroid_y_col=pta.get("centroid_y_col", "centroid_y"),
            cell_id_col=pta.get("cell_id_col", "cell"),
            force_rerun=force_rerun,
            proseg_run_script=(
                Path(pta["proseg_run_script"]) if pta.get("proseg_run_script") else None
            ),
            x_source=pta.get("x_source", "maxpost_counts"),
        )
        banner(f"stage {idx}/{n_stages}: proseg_to_anndata — complete in {time.time()-t0:.1f}s")

    # --- Stage: qc_filter → adds .obs['qc_filtered'] ---------------
    if "qc_filter" in stages:
        idx = stages.index("qc_filter") + 1
        banner(f"stage {idx}/{n_stages}: qc_filter — starting")
        t0 = time.time()
        from xenium_preprocess.stages.qc_filter import run_qc_filter
        qc_cfg = cfg.get("qc_filter", {}) or {}
        if not raw_h5ad.exists():
            raise SystemExit(
                f"qc_filter stage requested but the raw AnnData does not "
                f"exist: {raw_h5ad}. Run proseg_to_anndata first."
            )
        run_qc_filter(
            sample_id=sample_id,
            raw_h5ad=raw_h5ad,
            min_counts_cell=int(qc_cfg.get("min_counts_cell", 10)),
            force_rerun=force_rerun,
            qc_filtered_col=qc_cfg.get("qc_filtered_col", "qc_filtered"),
        )
        banner(f"stage {idx}/{n_stages}: qc_filter — complete in {time.time()-t0:.1f}s")

    # --- Stage: xenium_ranger_to_anndata --------------------------
    if "xenium_ranger_to_anndata" in stages:
        idx = stages.index("xenium_ranger_to_anndata") + 1
        banner(f"stage {idx}/{n_stages}: xenium_ranger_to_anndata — starting")
        t0 = time.time()
        from xenium_preprocess.stages.xenium_ranger_to_anndata import (
            run_xenium_ranger_to_anndata,
        )
        xr_cfg = cfg.get("xenium_ranger_to_anndata", {}) or {}
        run_xenium_ranger_to_anndata(
            sample_id=sample_id,
            run_id=run_id,
            xenium_ranger_dir=(
                Path(xr_cfg["xenium_ranger_dir"])
                if xr_cfg.get("xenium_ranger_dir") else None
            ),
            out_h5ad=xr_h5ad,
            gex_only=bool(xr_cfg.get("gex_only", True)),
            force_rerun=force_rerun,
            cells_csv_name=xr_cfg.get("cells_csv_name", "cells.csv.gz"),
            cell_feature_matrix_name=xr_cfg.get(
                "cell_feature_matrix_name", "cell_feature_matrix.h5"
            ),
        )
        banner(f"stage {idx}/{n_stages}: xenium_ranger_to_anndata — "
               f"complete in {time.time()-t0:.1f}s")

    # --- Stage: enrich_xenium_id → NN-map proseg cell_id onto -------
    # `<S>_xenium_ranger.h5ad`. Reversed on 2026-08-11 (user request,
    # (internal issue review): previously wrote xenium
    # ids onto the proseg h5ad. Runs AFTER xenium_ranger_to_anndata so
    # both source h5ads exist.
    if "enrich_xenium_id" in stages:
        idx = stages.index("enrich_xenium_id") + 1
        exid_cfg = cfg.get("enrich_xenium_id", {}) or {}
        if not bool(exid_cfg.get("enabled", True)):
            banner(f"stage {idx}/{n_stages}: enrich_xenium_id — "
                   f"disabled in config (enrich_xenium_id.enabled=false), skipping")
        else:
            banner(f"stage {idx}/{n_stages}: enrich_xenium_id — starting")
            t0 = time.time()
            from xenium_preprocess.stages.enrich_xenium_id import run_enrich_xenium_id
            if not xr_h5ad.exists():
                raise SystemExit(
                    f"enrich_xenium_id stage requested but the xenium-ranger "
                    f"AnnData does not exist: {xr_h5ad}. Run "
                    f"xenium_ranger_to_anndata first."
                )
            if not raw_h5ad.exists():
                raise SystemExit(
                    f"enrich_xenium_id stage requested but the proseg raw "
                    f"AnnData does not exist: {raw_h5ad}. Run "
                    f"proseg_to_anndata first."
                )
            run_enrich_xenium_id(
                sample_id=sample_id,
                run_id=run_id,
                xenium_ranger_h5ad=xr_h5ad,
                proseg_h5ad=raw_h5ad,
                output_root=output_root,
                nn_k=int(exid_cfg.get("nn_k", 1)),
                distance_threshold=exid_cfg.get("distance_threshold"),
                nn_id_col=exid_cfg.get("nn_id_col", "proseg_cell_id_nn"),
                distance_col=exid_cfg.get("distance_col", "proseg_id_nn_distance"),
                note_col=exid_cfg.get("note_col", "proseg_id_nn_note"),
                write_back_h5ad=bool(exid_cfg.get("write_back_h5ad", True)),
                force_rerun=force_rerun,
            )
            banner(f"stage {idx}/{n_stages}: enrich_xenium_id — "
                   f"complete in {time.time()-t0:.1f}s")

    # --- Sub-stage: preprocess (RETAINED for reversibility, NOT default) --
    # The heavy PCA/UMAP/Leiden/celltype pass is kept in the codebase
    # but is no longer part of DEFAULT_STAGES. Users who want the old
    # 4.4-GB `<S>_preprocessed.h5ad` output can `--stages ... preprocess`.
    if "preprocess" in stages:
        idx = stages.index("preprocess") + 1
        banner(f"stage {idx}/{n_stages}: preprocess — DEPRECATED, running anyway")
        t0 = time.time()
        from xenium_preprocess.stages.preprocess import run_preprocess
        pp = cfg.get("preprocess", {}) or {}
        if not raw_h5ad.exists():
            raise SystemExit(
                f"preprocess stage requested but the raw AnnData does not "
                f"exist: {raw_h5ad}. Run proseg_to_anndata first."
            )
        # `preprocess` is DEPRECATED (not in DEFAULT_STAGES). When
        # opted-in it writes to `<raw_h5ad.parent.parent>/legacy_preprocess/`
        # — i.e. under the run folder but out of the locked
        # `spatial_adata/` / `rctd/` shape, so it can't collide with any
        # of the six default stages' outputs. The stage's signature
        # is unchanged from the pre-refactor version.
        run_preprocess(
            sample_id=sample_id,
            raw_h5ad=raw_h5ad,
            output_root=output_root,
            qc_percent_top=tuple(pp.get("qc_percent_top", [10, 20, 50, 150])),
            min_counts_cell=pp.get("min_counts_cell", 10),
            min_prop=pp.get("min_prop", 1e-3),
            positive_x=bool(pp.get("positive_x", False)),
            n_neighbors=pp.get("n_neighbors", 15),
            random_state=pp.get("random_state", 0),
            leiden_resolution=pp.get("leiden_resolution", 0.4),
            use_leiden_weights=bool(pp.get("use_leiden_weights", False)),
            global_non_tumor_json=(
                Path(pp["global_non_tumor_json"]) if pp.get("global_non_tumor_json") else None
            ),
            global_tumor_json=(
                Path(pp["global_tumor_json"]) if pp.get("global_tumor_json") else None
            ),
            tumor_type=pp.get("tumor_type"),
            threshold_global=pp.get("threshold_global", -0.005),
            low_count_label=pp.get("low_count_label", "low count"),
            global_level1_celltype_col=pp.get(
                "global_level1_celltype_col", "global_level1_celltype"
            ),
            umap_dpi=pp.get("umap_dpi", 200),
            umap_figsize=tuple(pp.get("umap_figsize", [7, 6])),
            force_rerun=force_rerun,
            legacy_symlinks=False,
            dual_matrix_mode=bool(pp.get("dual_matrix_mode", True)),
        )
        banner(f"stage {idx}/{n_stages}: preprocess — complete in {time.time()-t0:.1f}s")

    # --- Stage: SPLIT triple (mtx/features/barcodes + sidecars) ---
    if "split_prep" in stages:
        idx = stages.index("split_prep") + 1
        banner(f"stage {idx}/{n_stages}: split_prep — starting")
        t0 = time.time()
        from xenium_preprocess.stages.split_prep import run_split_prep
        sp = cfg.get("split_prep", {}) or {}
        if not raw_h5ad.exists():
            raise SystemExit(
                f"split_prep stage requested but the raw AnnData does not "
                f"exist: {raw_h5ad}. Run proseg_to_anndata first."
            )
        rp_for_split = cfg.get("rctd_prep", {}) or {}
        split_layer = rp_for_split.get("source_layer") or sp.get(
            "layer", "maxpost_counts"
        )
        run_split_prep(
            sample_id=sample_id,
            run_id=run_id,
            raw_h5ad=raw_h5ad,
            output_root=output_root,
            layer=split_layer,
            name_suffix=sp.get("name_suffix", "_lowCountThreshold"),
            gzip_outputs=bool(sp.get("gzip_outputs", True)),
            min_counts_cell=sp.get("min_counts_cell", 10),
            force_rerun=force_rerun,
        )
        banner(f"stage {idx}/{n_stages}: split_prep — complete in {time.time()-t0:.1f}s")

    # --- Stage: RCTD test object ----------------------------------
    # rctd_prep's R script consumes the mtx/features/barcodes triple
    # written by split_prep. the user's 2026-08-12 clarification
    # (internal issue review): keep producing the RDS
    # by default, but don't leave the mtx bundle on disk. So when the
    # user asked for rctd_prep WITHOUT explicit split_prep, we run
    # split_prep as a TRANSIENT prerequisite here, then delete
    # `<run>/split_prep/` after the RDS is written. When split_prep is
    # explicitly in stages, the earlier block wrote the bundle and we
    # leave it in place.
    if "rctd_prep" in stages:
        import shutil

        idx = stages.index("rctd_prep") + 1
        banner(f"stage {idx}/{n_stages}: rctd_prep — starting")
        t0 = time.time()
        from xenium_preprocess.stages.rctd_prep import run_rctd_prep
        rp = cfg.get("rctd_prep", {}) or {}
        sp = cfg.get("split_prep", {}) or {}
        split_dir = the_run_dir / "split_prep"
        split_is_transient = "split_prep" not in stages

        stem = f"{sample_id}{sp.get('name_suffix', '_lowCountThreshold')}"
        gz = ".gz" if bool(sp.get("gzip_outputs", True)) else ""
        split_sentinel = split_dir / f"{stem}_counts.mtx{gz}"

        if not split_sentinel.exists() or force_rerun:
            if split_is_transient:
                log("[rctd_prep] split_prep bundle absent and split_prep "
                    "not in requested stages — producing it as a "
                    "TRANSIENT prerequisite (will be deleted after "
                    "rctd_prep writes the RDS).")
            if not raw_h5ad.exists():
                raise SystemExit(
                    f"rctd_prep stage requested but the raw AnnData "
                    f"does not exist: {raw_h5ad}. Run proseg_to_anndata "
                    f"first."
                )
            from xenium_preprocess.stages.split_prep import run_split_prep
            split_layer = rp.get("source_layer") or sp.get(
                "layer", "maxpost_counts"
            )
            run_split_prep(
                sample_id=sample_id,
                run_id=run_id,
                raw_h5ad=raw_h5ad,
                output_root=output_root,
                layer=split_layer,
                name_suffix=sp.get("name_suffix", "_lowCountThreshold"),
                gzip_outputs=bool(sp.get("gzip_outputs", True)),
                min_counts_cell=sp.get("min_counts_cell", 10),
                force_rerun=force_rerun,
            )

        run_rctd_prep(
            sample_id=sample_id,
            run_id=run_id,
            split_dir=split_dir,
            output_root=output_root,
            name_suffix=sp.get("name_suffix", "_lowCountThreshold"),
            rscript_bin=rp.get("rscript_bin", "Rscript"),
            assay_name=rp.get("assay_name", "Proseg"),
            spatial_key=rp.get("spatial_key", "ST_"),
            force_rerun=force_rerun,
        )

        if split_is_transient and split_dir.exists():
            log(f"[rctd_prep] cleaning up transient split_prep bundle: "
                f"{split_dir}")
            shutil.rmtree(split_dir)

        banner(f"stage {idx}/{n_stages}: rctd_prep — complete in {time.time()-t0:.1f}s")

    banner(f"xenium-preprocess done in {time.time()-pipeline_t0:.1f}s "
           f"-> {the_run_dir}/")
    return 0
