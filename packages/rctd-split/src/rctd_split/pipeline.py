"""Stage-orchestration loop for rctd-split (step 4).

Runs the requested subset of the pipeline stages, each guarded by its
own sentinel-existence resume check (nuke with ``force_rerun``).

Post-2026-08-11 refactor: outputs live under
``<output_root>/<sample_id>/<sample_id>_<run_id>/{spatial_adata,rctd,intermediate}/``,
and ``config.yaml`` at that run folder is written via the shared
``_merge_config`` helper so step 1 / step 3 (spawned as separate sbatch
jobs in the driver chain) don't clobber each other's config sections.

After the terminal stage (``qc_report``) succeeds the entire
``intermediate/`` subfolder is DROPPED (mirrors ref-build; the
final artifacts under ``spatial_adata/``, ``rctd/``, and ``summary/``
survive). Pass ``--keep-intermediate`` / ``keep_intermediate: true``
to retain it for debugging / QC or for a follow-up
``--stages writeback_to_step1_raw`` re-run
(settylab/TracyY123-nexus#26 comment 5322401335).
"""
from __future__ import annotations

import os
import shutil
import socket
import sys
import time
from pathlib import Path

from rctd_split._internal.layout import (
    intermediate_dir,
    resolved_config_path,
    rctd_path,
    run_dir,
)
from rctd_split._internal.logging import banner, log
from rctd_split._internal.merge_config import merge_config


def _drop_intermediate_outputs(
    output_root: Path, sample_id: str, run_id: str,
) -> None:
    """Drop the entire `<run_dir>/intermediate/` directory after the
    terminal stage (`qc_report`) succeeds.

    Removes the heavy R-side `.rds` working files under `split/`, the
    10X-style mtx bundles under `mtx/`, and the intermediate adata
    h5ads / csvs / sentinels under `adata/`. The persisted final
    artifacts under `spatial_adata/`, `rctd/`, and `summary/` are
    untouched.

    Users who need to re-run `--stages writeback_to_step1_raw` against
    an already-completed run without redoing SPLIT/mtx/adata should
    pass `--keep-intermediate` / `keep_intermediate: true` at the
    original invocation to opt out of cleanup (matches Tracy's ask on
    settylab/TracyY123-nexus#26 comment 5322401335). Idempotent — no-op
    when the folder is already gone."""
    p = intermediate_dir(output_root, sample_id, run_id)
    if not p.exists():
        return
    log(f"[pipeline] dropping intermediate/: {p}")
    shutil.rmtree(p)
    log("[pipeline] intermediate/ pruned; pass --keep-intermediate to preserve "
        "for `--stages writeback_to_step1_raw` re-runs.")


def _raise_missing_input(
    *,
    path: Path,
    source: str,
    layout_hint: str,
    producer: str,
    override_flag: str,
) -> None:
    """Uniform fail-loud for a missing test_object / reference_rds.

    Message shape adapts to how the path got resolved: `source="config"`
    means the caller passed it explicitly; `source="layout"` means we
    derived it from the run-folder convention. Points to both remedies
    either way — produce the file via the upstream step, or override
    with the CLI flag — so the operator can pick.
    """
    if source == "layout":
        raise SystemExit(
            f"input not found: {path}\n"
            f"  source: layout-derived ({layout_hint} — {producer} output).\n"
            f"  Either produce it via {producer}, "
            f"or pass {override_flag} explicitly."
        )
    raise SystemExit(
        f"input not found: {path}\n"
        f"  source: explicit ({override_flag} / config).\n"
        f"  Check the path, or omit {override_flag} to let the pipeline "
        f"derive it from the run-folder layout ({layout_hint})."
    )


def log_invocation_banner(argv: list[str], stages: list[str], cfg: dict) -> None:
    """Print everything anyone would want to grep out of the log."""
    banner("rctd-split starting")
    log(f"command:       {' '.join(argv)}")
    log(f"host:          {socket.gethostname()}")
    log(f"cwd:           {os.getcwd()}")
    log(f"pid:           {os.getpid()}")
    log(f"python:        {sys.executable}")
    log(f"python ver:    {sys.version.splitlines()[0]}")
    log(f"sample_id:     {cfg.get('sample_id')}")
    log(f"run_id:        {cfg.get('run_id')}")
    log(f"test_object:   {cfg.get('test_object')}")
    log(f"reference_rds: {cfg.get('reference_rds')}")
    log(f"output_root:   {cfg.get('output_root')}")
    log(f"stages:        {stages}")
    log(f"force_rerun:   {cfg.get('force_rerun', False)}")
    for name in ("numpy", "yaml", "scipy", "pandas", "scanpy", "anndata"):
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
    output_root = Path(cfg["output_root"]).resolve()
    force_rerun = bool(cfg.get("force_rerun", False))
    keep_intermediate = bool(cfg.get("keep_intermediate", False))

    # `test_object` and `reference_rds` auto-derive from the run-folder
    # layout convention when not explicitly set — step 1 (rctd_prep)
    # writes to <run-dir>/rctd/<S>_test_object.rds and step 3
    # (rctd_reference_build) writes to <run-dir>/rctd/<S>_reference.rds,
    # so the workflow-driver chain wires up with no explicit paths.
    # Explicit config values still win (backward compat + one-off
    # experiments with a foreign test-object or reference).
    if cfg.get("test_object"):
        test_object = Path(cfg["test_object"]).resolve()
        test_object_source = "config"
    else:
        test_object = rctd_path(output_root, sample_id, run_id, "test_object")
        test_object_source = "layout"
    if cfg.get("reference_rds"):
        reference_rds = Path(cfg["reference_rds"]).resolve()
        reference_rds_source = "config"
    else:
        reference_rds = rctd_path(output_root, sample_id, run_id, "reference")
        reference_rds_source = "layout"

    # Record the resolved paths back into cfg so the invocation banner
    # and merged resolved_config snapshot show what was actually used.
    cfg = dict(cfg)
    cfg["test_object"] = str(test_object)
    cfg["reference_rds"] = str(reference_rds)

    log_invocation_banner(argv, stages, cfg)
    log(f"[pipeline] test_object   source: {test_object_source}")
    log(f"[pipeline] reference_rds source: {reference_rds_source}")

    if not test_object.exists():
        _raise_missing_input(
            path=test_object,
            source=test_object_source,
            layout_hint="<output_root>/<sample>/<sample>_<run_id>/rctd/<sample>_test_object.rds",
            producer="step 1's rctd_prep",
            override_flag="--test-object",
        )
    if not reference_rds.exists():
        _raise_missing_input(
            path=reference_rds,
            source=reference_rds_source,
            layout_hint="<output_root>/<sample>/<sample>_<run_id>/rctd/<sample>_reference.rds",
            producer="step 3's rctd_reference_build",
            override_flag="--reference-rds",
        )

    # Materialise the run folder up front.
    the_run_dir = run_dir(output_root, sample_id, run_id)
    the_run_dir.mkdir(parents=True, exist_ok=True)
    log(f"[pipeline] run folder: {the_run_dir}")

    # Merged resolved_config: write ONLY the step4 key, preserve any
    # step1/step3/driver sections a sibling pipeline already wrote.
    snap = resolved_config_path(output_root, sample_id, run_id)
    step4_cfg = dict(cfg)
    merge_config(snap, step_key="step4", step_cfg=step4_cfg)
    log(f"[pipeline] merged resolved_config -> {snap}")

    n_stages = len(stages)
    pipeline_t0 = time.time()
    ran_terminal = False

    rr_cfg = cfg.get("rctd_run", {}) or {}
    sp_cfg = cfg.get("split_purify", {}) or {}
    em_cfg = cfg.get("export_mtx", {}) or {}
    mh_cfg = cfg.get("mtx_to_h5ad", {}) or {}
    fs_cfg = cfg.get("filter_status", {}) or {}
    pp_cfg = cfg.get("postprocess", {}) or {}
    wb_cfg = cfg.get("writeback_to_step1_raw", {}) or {}
    ct_cfg = cfg.get("celltype_writeback", {}) or {}
    qc_cfg = cfg.get("qc_report", {}) or {}

    # --- Stage 1: rctd_run ------------------------------------------
    if "rctd_run" in stages:
        idx = stages.index("rctd_run") + 1
        banner(f"stage {idx}/{n_stages}: rctd_run — starting")
        t0 = time.time()
        from rctd_split.stages.rctd_run import run_rctd_run
        run_rctd_run(
            sample_id=sample_id,
            run_id=run_id,
            test_object=test_object,
            reference_rds=reference_rds,
            output_root=output_root,
            rscript_bin=rr_cfg.get("rscript_bin", "Rscript"),
            assay_name=rr_cfg.get("assay_name", "Proseg"),
            spatial_reduction=rr_cfg.get("spatial_reduction", "spatial"),
            UMI_min=rr_cfg.get("UMI_min", 10),
            counts_MIN=rr_cfg.get("counts_MIN", 10),
            UMI_min_sigma=rr_cfg.get("UMI_min_sigma", 100),
            max_cores=rr_cfg.get("max_cores", 4),
            CELL_MIN_INSTANCE=rr_cfg.get("CELL_MIN_INSTANCE", 20),
            doublet_mode=rr_cfg.get("doublet_mode", "doublet"),
            label_slash_replacement=rr_cfg.get("label_slash_replacement", "_"),
            force_rerun=force_rerun,
            r_lib_paths=cfg.get("r_lib_paths", []) or [],
        )
        banner(f"stage {idx}/{n_stages}: rctd_run — complete in {time.time()-t0:.1f}s")

    rctd_results_rds = rctd_path(output_root, sample_id, run_id, "rctd_results")

    # --- Stage 2: split_purify --------------------------------------
    if "split_purify" in stages:
        idx = stages.index("split_purify") + 1
        banner(f"stage {idx}/{n_stages}: split_purify — starting")
        t0 = time.time()
        from rctd_split.stages.split_purify import run_split_purify
        if not rctd_results_rds.exists():
            raise SystemExit(
                f"split_purify stage requested but rctd_results.rds does not exist: "
                f"{rctd_results_rds}. Run the rctd_run stage first."
            )
        unpurified_rds, purified_rds = run_split_purify(
            sample_id=sample_id,
            run_id=run_id,
            test_object=test_object,
            rctd_results_rds=rctd_results_rds,
            output_root=output_root,
            rscript_bin=sp_cfg.get("rscript_bin", "Rscript"),
            assay_name=sp_cfg.get("assay_name", "Proseg"),
            DO_purify_singlets=bool(sp_cfg.get("DO_purify_singlets", True)),
            force_rerun=force_rerun,
            r_lib_paths=cfg.get("r_lib_paths", []) or [],
        )
        banner(f"stage {idx}/{n_stages}: split_purify — complete in {time.time()-t0:.1f}s")
    else:
        from rctd_split._internal.layout import intermediate_path
        unpurified_rds = intermediate_path(output_root, sample_id, run_id, "unpurified_rds")
        purified_rds = intermediate_path(output_root, sample_id, run_id, "purified_rds")

    # --- Stage 3: export_mtx ----------------------------------------
    if "export_mtx" in stages:
        idx = stages.index("export_mtx") + 1
        banner(f"stage {idx}/{n_stages}: export_mtx — starting")
        t0 = time.time()
        from rctd_split.stages.export_mtx import run_export_mtx
        for name, path in (("unpurified.rds", unpurified_rds),
                           ("purified.rds", purified_rds)):
            if not path.exists():
                raise SystemExit(
                    f"export_mtx stage requested but {name} does not exist: "
                    f"{path}. Run the split_purify stage first."
                )
        unpurified_mtx_dir, purified_mtx_dir = run_export_mtx(
            sample_id=sample_id,
            run_id=run_id,
            unpurified_rds=unpurified_rds,
            purified_rds=purified_rds,
            output_root=output_root,
            rscript_bin=em_cfg.get("rscript_bin", "Rscript"),
            assay_name=em_cfg.get("assay_name", "Proseg"),
            gzip_outputs=bool(em_cfg.get("gzip_outputs", True)),
            force_rerun=force_rerun,
            r_lib_paths=cfg.get("r_lib_paths", []) or [],
        )
        banner(f"stage {idx}/{n_stages}: export_mtx — complete in {time.time()-t0:.1f}s")
    else:
        from rctd_split._internal.layout import mtx_bundle_dir
        unpurified_mtx_dir = mtx_bundle_dir(output_root, sample_id, run_id, "unpurified")
        purified_mtx_dir = mtx_bundle_dir(output_root, sample_id, run_id, "purified")

    # --- Stage 4: mtx_to_h5ad ---------------------------------------
    if "mtx_to_h5ad" in stages:
        idx = stages.index("mtx_to_h5ad") + 1
        banner(f"stage {idx}/{n_stages}: mtx_to_h5ad — starting")
        t0 = time.time()
        from rctd_split.stages.mtx_to_h5ad import run_mtx_to_h5ad
        for name, path in (("unpurified mtx bundle", unpurified_mtx_dir),
                           ("purified mtx bundle", purified_mtx_dir)):
            if not path.exists():
                raise SystemExit(
                    f"mtx_to_h5ad stage requested but the {name} directory "
                    f"does not exist: {path}. Run the export_mtx stage first."
                )
        run_mtx_to_h5ad(
            sample_id=sample_id,
            run_id=run_id,
            unpurified_mtx_dir=unpurified_mtx_dir,
            purified_mtx_dir=purified_mtx_dir,
            output_root=output_root,
            h5ad_compression=mh_cfg.get("h5ad_compression", "gzip"),
            force_rerun=force_rerun,
        )
        banner(f"stage {idx}/{n_stages}: mtx_to_h5ad — complete in {time.time()-t0:.1f}s")

    # --- Stage 5: filter_status -------------------------------------
    if "filter_status" in stages:
        idx = stages.index("filter_status") + 1
        banner(f"stage {idx}/{n_stages}: filter_status — starting")
        t0 = time.time()
        from rctd_split.stages.filter_status import (
            DEFAULT_INPLACE_COLUMNS,
            run_filter_status,
        )
        run_filter_status(
            sample_id=sample_id,
            run_id=run_id,
            output_root=output_root,
            write_inplace=bool(fs_cfg.get("write_inplace", True)),
            inplace_columns=tuple(
                fs_cfg.get("inplace_columns", DEFAULT_INPLACE_COLUMNS)
            ),
            first_type_column=fs_cfg.get("first_type_column", "first_type"),
            h5ad_compression=fs_cfg.get(
                "h5ad_compression", mh_cfg.get("h5ad_compression", "gzip")
            ),
            force_rerun=force_rerun,
        )
        banner(f"stage {idx}/{n_stages}: filter_status — complete in {time.time()-t0:.1f}s")

    # --- Stage 6: postprocess ---------------------------------------
    if "postprocess" in stages:
        idx = stages.index("postprocess") + 1
        banner(f"stage {idx}/{n_stages}: postprocess — starting")
        t0 = time.time()
        from rctd_split.stages.postprocess import (
            DEFAULT_LEIDEN_RESOLUTIONS,
            run_postprocess,
        )
        pp_qc = pp_cfg.get("qc", {}) or {}
        pp_norm = pp_cfg.get("normalize", {}) or {}
        pp_nn = pp_cfg.get("neighbors", {}) or {}
        pp_ld = pp_cfg.get("leiden", {}) or {}
        run_postprocess(
            sample_id=sample_id,
            run_id=run_id,
            output_root=output_root,
            min_counts=int(pp_qc.get("min_counts", 50)),
            minprop=float(pp_norm.get("minprop", 1.0e-3)),
            positive_x=bool(pp_norm.get("positive_x", True)),
            n_neighbors=int(pp_nn.get("n_neighbors", 15)),
            random_state=int(pp_nn.get("random_state", 0)),
            leiden_resolutions=pp_ld.get(
                "resolutions", list(DEFAULT_LEIDEN_RESOLUTIONS),
            ),
            process_unpurified_layer=bool(
                pp_cfg.get("process_unpurified_layer", True)
            ),
            h5ad_compression=pp_cfg.get(
                "h5ad_compression", mh_cfg.get("h5ad_compression", "gzip"),
            ),
            force_rerun=force_rerun,
        )
        banner(f"stage {idx}/{n_stages}: postprocess — complete in {time.time()-t0:.1f}s")

    # --- Stage 7: writeback_to_step1_raw ---------------------------
    if "writeback_to_step1_raw" in stages:
        idx = stages.index("writeback_to_step1_raw") + 1
        banner(f"stage {idx}/{n_stages}: writeback_to_step1_raw — starting")
        t0 = time.time()
        from rctd_split.stages.writeback_to_step1_raw import (
            run_writeback_to_step1_raw,
        )
        raw_h5ad_cfg = wb_cfg.get("raw_h5ad")
        run_writeback_to_step1_raw(
            sample_id=sample_id,
            run_id=run_id,
            output_root=output_root,
            raw_h5ad=Path(raw_h5ad_cfg) if raw_h5ad_cfg else None,
            qc_filtered_col=wb_cfg.get("qc_filtered_col", "qc_filtered"),
            exclude_obs_cols=wb_cfg.get("exclude_obs_cols", []) or [],
            h5ad_compression=wb_cfg.get(
                "h5ad_compression", mh_cfg.get("h5ad_compression", "gzip"),
            ),
            force_rerun=force_rerun,
        )
        banner(f"stage {idx}/{n_stages}: writeback_to_step1_raw — "
               f"complete in {time.time()-t0:.1f}s")

    # --- Stage 8: celltype_writeback -------------------------------
    if "celltype_writeback" in stages:
        idx = stages.index("celltype_writeback") + 1
        banner(f"stage {idx}/{n_stages}: celltype_writeback — starting")
        t0 = time.time()
        from rctd_split.stages.celltype_writeback import run_celltype_writeback
        xen_h5ad_cfg = ct_cfg.get("xenium_ranger_h5ad")
        run_celltype_writeback(
            sample_id=sample_id,
            run_id=run_id,
            output_root=output_root,
            xenium_ranger_h5ad=Path(xen_h5ad_cfg) if xen_h5ad_cfg else None,
            reference_label_col=ct_cfg.get("reference_label_col", "first_type"),
            reference_x_col=ct_cfg.get("reference_x_col", "centroid_x"),
            reference_y_col=ct_cfg.get("reference_y_col", "centroid_y"),
            k=int(ct_cfg.get("k", 1)),
            algorithm=ct_cfg.get("algorithm", "auto"),
            metric=ct_cfg.get("metric", "euclidean"),
            tiebreak=ct_cfg.get("tiebreak", "min_dist"),
            distance_threshold=ct_cfg.get("distance_threshold"),
            unmatched_policy=ct_cfg.get("unmatched_policy", "nearest_label"),
            celltype_col=ct_cfg.get("celltype_col", "celltype"),
            distance_col=ct_cfg.get("distance_col", "celltype_source_distance"),
            h5ad_compression=ct_cfg.get(
                "h5ad_compression", mh_cfg.get("h5ad_compression", "gzip"),
            ),
            force_rerun=force_rerun,
        )
        banner(f"stage {idx}/{n_stages}: celltype_writeback — "
               f"complete in {time.time()-t0:.1f}s")

    # --- Stage 9: qc_report ----------------------------------------
    if "qc_report" in stages:
        idx = stages.index("qc_report") + 1
        banner(f"stage {idx}/{n_stages}: qc_report — starting")
        t0 = time.time()
        from rctd_split.stages.qc_report import run_qc_report
        run_qc_report(
            sample_id=sample_id,
            run_id=run_id,
            output_root=output_root,
            purification_status_column=qc_cfg.get(
                "purification_status_column", "purification_status",
            ),
            raw_layer=qc_cfg.get("raw_layer", "maxpost_counts"),
            force_rerun=force_rerun,
            invoking_argv=argv,
            extra_reports=qc_cfg.get("extra_reports", []) or [],
        )
        banner(f"stage {idx}/{n_stages}: qc_report — "
               f"complete in {time.time()-t0:.1f}s")
        ran_terminal = True

    # --- Drop intermediates -----------------------------------------
    # Only if the terminal stage (qc_report) actually ran in THIS
    # invocation. A partial-stages run (e.g. `--stages postprocess`)
    # leaves `intermediate/` in place so a follow-up invocation can
    # pick up from the sentinels. `keep_intermediate=true` opts out
    # entirely (debugging / QC access to the intermediate mtx / h5ad).
    if ran_terminal and not keep_intermediate:
        _drop_intermediate_outputs(output_root, sample_id, run_id)
    elif ran_terminal and keep_intermediate:
        log("[pipeline] keep_intermediate=true — leaving "
            f"{intermediate_dir(output_root, sample_id, run_id)}/ in place")

    banner(f"rctd-split done in {time.time()-pipeline_t0:.1f}s "
           f"-> {the_run_dir}/")
    return 0
