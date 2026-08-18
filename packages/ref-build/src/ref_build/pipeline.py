"""Sub-stage orchestration loop for ref-build (the middle stage of the
three-stage pipeline).

Runs the requested subset of
    load_primary_and_donors → census → assemble → export_mtx → rctd_reference_build

All outputs — intermediate and final — land under the run-scoped
directory
    <output_root>/<sample_id>/<sample_id>_<run_id>/
(internal issue review).

After `rctd_reference_build` succeeds the intermediate outputs
(`loaded/`, `census/`, `mtx_bundle/`) are DROPPED. The two final
artifacts under `rctd/` — `<sample_id>_reference_post_rules.h5ad` (from
`assemble`) and `<sample_id>_reference.rds` (from `rctd_reference_build`)
— survive; they are the inputs rctd-split consumes.

Merged `config.yaml` is written via the shared
`_internal.merge_config` helper under the `ref_build:` top-level key so
xenium-preprocess / rctd-split co-writing the same file don't clobber
each other's sections.
"""
from __future__ import annotations

import os
import shutil
import socket
import sys
import time
from pathlib import Path

from ref_build._internal.layout import (
    rctd_dir,
    rctd_path,
    resolved_config_path,
    run_dir,
)
from ref_build._internal.logging import log, banner
from ref_build._internal.merge_config import merge_config


# Intermediate subdirs under <run_dir>/ that are DROPPED after
# rctd_reference_build succeeds (internal issue review).
_INTERMEDIATE_SUBDIRS = ("loaded", "census", "mtx_bundle")


def log_invocation_banner(argv: list[str], stages: list[str], cfg: dict) -> None:
    """Print everything anyone would want to grep out of the log if the
    job ever silently does nothing again."""
    banner("ref-build starting")
    log(f"command:              {' '.join(argv)}")
    log(f"host:                 {socket.gethostname()}")
    log(f"cwd:                  {os.getcwd()}")
    log(f"pid:                  {os.getpid()}")
    log(f"python:               {sys.executable}")
    log(f"python ver:           {sys.version.splitlines()[0]}")
    log(f"sample_id:            {cfg.get('sample_id')}")
    log(f"run_id:               {cfg.get('run_id')}")
    log(f"primary_h5ad:         {cfg.get('primary_h5ad')}")
    log(f"flex_h5ad_path:       {cfg.get('flex_h5ad_path')}")
    log(f"donor_h5ads:          {cfg.get('donor_h5ads')}")
    log(f"fallback_donor_h5ads: {cfg.get('fallback_donor_h5ads')}")
    log(f"celltype_marker_json: {cfg.get('celltype_marker_json')}")
    cen_dbg = cfg.get("census", {}) or {}
    log(f"celltype_target_list: {cen_dbg.get('celltype_target_list')}")
    log(f"celltype_target_key:  {cen_dbg.get('celltype_target_key')}")
    log(f"tumor_type:           {cfg.get('tumor_type')}")
    log(f"output_root:          {cfg.get('output_root')}")
    log(f"stages:               {stages}")
    log(f"force_rerun:          {cfg.get('force_rerun', False)}")
    cen = cfg.get("census", {}) or {}
    log(f"donor_borrow_cap:     {cen.get('donor_borrow_cap')}")
    log(f"cell_min_instance:    {cen.get('cell_min_instance')}")
    log(f"primary_only_ct:      {cen.get('primary_only_celltypes')}")
    log(f"random_seed:          {cen.get('random_seed')}")
    rr = cfg.get("rctd_reference_build", {}) or {}
    log(f"min_UMI:              {rr.get('min_UMI')}")
    log(f"require_int:          {rr.get('require_int')}")
    # Library versions — wrap each in its own try so one missing library
    # doesn't suppress the rest.
    for name in ("numpy", "yaml", "scipy", "scanpy", "anndata", "pandas"):
        try:
            mod = __import__(name)
            ver = getattr(mod, "__version__", "no __version__ attr")
            log(f"  {name:14s}{ver}")
        except Exception as e:
            log(f"  {name:14s}NOT IMPORTABLE ({type(e).__name__})")
    log("=" * 64)


def _drop_intermediate_outputs(rundir: Path) -> None:
    """Drop the three intermediate output subdirs after
    `rctd_reference_build` succeeds. Idempotent — no-op when a subdir is
    already gone (e.g. from a prior sentinel-skipped run)."""
    for sub in _INTERMEDIATE_SUBDIRS:
        p = rundir / sub
        if p.exists():
            log(f"[pipeline] dropping intermediate outputs: {p}")
            shutil.rmtree(p)


def _write_merged_resolved_config(
    output_root: Path, sample_id: str, run_id: str, cfg: dict,
) -> Path:
    """Write only the `ref_build:` top-level key to the shared merged
    `config.yaml` at `<run_dir>/config.yaml`.
    Preserves any sibling `xenium_preprocess:` / `rctd_split:` / `driver:`
    sections written by adjacent pipelines in the same run.

    The `ref_build:` payload we record is the FULL resolved config for this
    invocation — including `flex_h5ad_path` (reference-in-place — no
    copy, no symlink), the resolved `run_id`, and every sub-stage's knobs.
    """
    resolved = resolved_config_path(output_root, sample_id, run_id)
    ref_build_cfg = dict(cfg)  # shallow copy is enough — we don't mutate below
    merged = merge_config(resolved, "ref_build", ref_build_cfg)
    log(f"[pipeline] resolved config -> {resolved} "
        f"(top-level keys after merge: {sorted(merged.keys())})")
    return resolved


def run(cfg: dict, stages: list[str], argv: list[str]) -> int:
    """Execute the pipeline. Assumes `cfg` has already been validated."""
    sample_id = cfg["sample_id"]
    run_id = cfg["run_id"]
    if not run_id:
        raise SystemExit(
            "[pipeline] cfg['run_id'] must be a non-empty string — the CLI's "
            "_resolve_run_id should have bound it (--run-id > $SLURM_JOB_ID > "
            "timestamp)."
        )
    primary_h5ad = Path(cfg["primary_h5ad"]).resolve()
    donor_h5ads = [Path(p).resolve() for p in (cfg.get("donor_h5ads") or [])]
    fallback_donor_h5ads = [
        Path(p).resolve() for p in (cfg.get("fallback_donor_h5ads") or [])
    ]
    marker_json = Path(cfg["celltype_marker_json"]).resolve()
    output_root = Path(cfg["output_root"]).resolve()
    force_rerun = bool(cfg.get("force_rerun", False))
    tumor_type = cfg.get("tumor_type")

    log_invocation_banner(argv, stages, cfg)

    if not primary_h5ad.exists():
        raise SystemExit(f"primary_h5ad not found: {primary_h5ad}")
    for p in donor_h5ads:
        if not p.exists():
            raise SystemExit(f"donor_h5ad not found: {p}")
    for p in fallback_donor_h5ads:
        if not p.exists():
            raise SystemExit(f"fallback_donor_h5ad not found: {p}")
    if not marker_json.exists():
        raise SystemExit(f"celltype_marker_json not found: {marker_json}")

    # Optional target-list JSON — separate from the marker JSON. When
    # set, its keys override `celltype_marker_json`'s keys as the
    # "required celltype set" the census iterates over (Rule 5 target).
    cen_cfg = cfg.get("census", {}) or {}
    target_list_val = cen_cfg.get("celltype_target_list")
    target_list_path = Path(target_list_val).resolve() if target_list_val else None
    if target_list_path is not None and not target_list_path.exists():
        raise SystemExit(f"celltype_target_list not found: {target_list_path}")
    target_list_key = cen_cfg.get("celltype_target_key")

    # Materialise the run-scoped output directory up front so every
    # stage + the merged resolved config find their parent already.
    rundir = run_dir(output_root, sample_id, run_id)
    rundir.mkdir(parents=True, exist_ok=True)
    rctd_dir(output_root, sample_id, run_id).mkdir(parents=True, exist_ok=True)

    # Merged resolved-config write BEFORE any sub-stage runs. This lets a
    # downstream consumer read `ref_build.flex_h5ad_path` even after a
    # partial-substages invocation.
    _write_merged_resolved_config(output_root, sample_id, run_id, cfg)

    n_stages = len(stages)
    pipeline_t0 = time.time()

    # Per-stage sentinels — declared up front so downstream stages can
    # existence-check their inputs.
    concat_h5ad = rundir / "loaded" / "concat.h5ad"
    census_csv = rundir / "census" / "census.csv"
    reference_h5ad = rctd_path(output_root, sample_id, run_id, "reference_post_rules")
    mtx_dir = rundir / "mtx_bundle"

    ran_rctd_reference_build = False

    # --- Stage 1: load_primary_and_donors ---------------------------
    if "load_primary_and_donors" in stages:
        idx = stages.index("load_primary_and_donors") + 1
        banner(f"stage {idx}/{n_stages}: load_primary_and_donors — starting")
        t0 = time.time()
        from ref_build.stages.load_primary_and_donors import run_load_primary_and_donors
        lpd = cfg.get("load_primary_and_donors", {}) or {}
        run_load_primary_and_donors(
            sample_id=sample_id,
            primary_h5ad=primary_h5ad,
            donor_h5ads=donor_h5ads,
            output_root=output_root,
            celltype_col=lpd.get("celltype_col", "Final_level1_celltype_annotation"),
            force_rerun=force_rerun,
            run_id=run_id,
            tumor_type=tumor_type,
            celltype_col_fallbacks=lpd.get(
                "celltype_col_fallbacks",
                ["global_level1_celltype_annotation"],
            ),
            fallback_donor_h5ads=fallback_donor_h5ads,
        )
        banner(f"stage {idx}/{n_stages}: load_primary_and_donors — complete in {time.time()-t0:.1f}s")

    # --- Stage 2: census --------------------------------------------
    if "census" in stages:
        idx = stages.index("census") + 1
        cen = cfg.get("census", {}) or {}
        banner(f"stage {idx}/{n_stages}: census — starting")
        t0 = time.time()
        if not concat_h5ad.exists():
            raise SystemExit(
                f"census stage requested but the concat AnnData does not exist: "
                f"{concat_h5ad}. Run the load_primary_and_donors stage first."
            )
        # Hybrid rule set (internal issue review).
        # `rule1_threshold` is a deprecated alias for `donor_borrow_cap`;
        # honor it if a legacy config sets ONLY the old key.
        donor_borrow_cap = int(cen.get(
            "donor_borrow_cap",
            cen.get("rule1_threshold", 100),
        ))
        cell_min_instance = int(cen.get("cell_min_instance", 20))
        from ref_build.stages.census import run_census
        run_census(
            sample_id=sample_id,
            concat_h5ad=concat_h5ad,
            output_root=output_root,
            celltype_marker_json=marker_json,
            celltype_col=cen.get("celltype_col", "Final_level1_celltype_annotation"),
            donor_borrow_cap=donor_borrow_cap,
            cell_min_instance=cell_min_instance,
            primary_only_celltypes=cen.get("primary_only_celltypes", ["tumor", "liver"]),
            force_rerun=force_rerun,
            run_id=run_id,
            fuzzy_matching=bool(cen.get("fuzzy_matching", True)),
            include_unknown_maybe=bool(cen.get("include_unknown_maybe", True)),
            celltype_target_list=target_list_path,
            celltype_target_key=target_list_key,
        )
        banner(f"stage {idx}/{n_stages}: census — complete "
               f"in {time.time()-t0:.1f}s")

    # --- Stage 3: assemble ------------------------------------------
    if "assemble" in stages:
        idx = stages.index("assemble") + 1
        banner(f"stage {idx}/{n_stages}: assemble — starting")
        t0 = time.time()
        from ref_build.stages.assemble import run_assemble
        asm = cfg.get("assemble", {}) or {}
        cen = cfg.get("census", {}) or {}
        if not concat_h5ad.exists():
            raise SystemExit(
                f"assemble stage requested but the concat AnnData does not exist: "
                f"{concat_h5ad}. Run the load_primary_and_donors stage first."
            )
        if not census_csv.exists():
            raise SystemExit(
                f"assemble stage requested but census.csv does not exist: "
                f"{census_csv}. Run the census stage first."
            )
        donor_borrow_cap = int(cen.get(
            "donor_borrow_cap",
            cen.get("rule1_threshold", 100),
        ))
        random_seed = int(cen.get("random_seed", cen.get("random_state", 42)))
        run_assemble(
            sample_id=sample_id,
            concat_h5ad=concat_h5ad,
            census_csv=census_csv,
            output_root=output_root,
            celltype_col=asm.get("celltype_col", "Final_level1_celltype_annotation"),
            donor_borrow_cap=donor_borrow_cap,
            random_state=random_seed,
            force_rerun=force_rerun,
            run_id=run_id,
            tumor_type=tumor_type,
            # Mirror the census's fuzzy flags. Assemble MUST see the same
            # rule set that produced the census counts — a mismatch would
            # let a cell be counted (census) then skipped (assemble), so
            # we read the same `census` block, not `assemble`.
            fuzzy_matching=bool(cen.get("fuzzy_matching", True)),
            include_unknown_maybe=bool(cen.get("include_unknown_maybe", True)),
        )
        banner(f"stage {idx}/{n_stages}: assemble — complete in {time.time()-t0:.1f}s")

    # --- Stage 4: export_mtx ----------------------------------------
    if "export_mtx" in stages:
        idx = stages.index("export_mtx") + 1
        banner(f"stage {idx}/{n_stages}: export_mtx — starting")
        t0 = time.time()
        from ref_build.stages.export_mtx import run_export_mtx
        ex = cfg.get("export_mtx", {}) or {}
        if not reference_h5ad.exists():
            raise SystemExit(
                f"export_mtx stage requested but the assembled reference does not exist: "
                f"{reference_h5ad}. Run the assemble stage first."
            )
        run_export_mtx(
            sample_id=sample_id,
            reference_h5ad=reference_h5ad,
            output_root=output_root,
            layer=ex.get("layer", "counts"),
            gzip_outputs=bool(ex.get("gzip_outputs", True)),
            force_rerun=force_rerun,
            run_id=run_id,
        )
        banner(f"stage {idx}/{n_stages}: export_mtx — complete in {time.time()-t0:.1f}s")

    # --- Stage 5: rctd_reference_build ------------------------------
    if "rctd_reference_build" in stages:
        idx = stages.index("rctd_reference_build") + 1
        banner(f"stage {idx}/{n_stages}: rctd_reference_build — starting")
        t0 = time.time()
        from ref_build.stages.rctd_reference_build import run_rctd_reference_build
        rr = cfg.get("rctd_reference_build", {}) or {}
        cen = cfg.get("census", {}) or {}
        if not mtx_dir.exists():
            raise SystemExit(
                f"rctd_reference_build stage requested but the mtx bundle dir does "
                f"not exist: {mtx_dir}. Run the export_mtx stage first."
            )
        run_rctd_reference_build(
            sample_id=sample_id,
            mtx_dir=mtx_dir,
            output_root=output_root,
            rscript_bin=rr.get("rscript_bin", "Rscript"),
            min_umi=rr.get("min_UMI", 10),
            require_int=bool(rr.get("require_int", True)),
            celltype_col=cen.get("celltype_col", "Final_level1_celltype_annotation"),
            label_slash_replacement=rr.get("label_slash_replacement", "_"),
            force_rerun=force_rerun,
            run_id=run_id,
        )
        ran_rctd_reference_build = True
        banner(f"stage {idx}/{n_stages}: rctd_reference_build — complete in {time.time()-t0:.1f}s")

    # --- Drop intermediates -----------------------------------------
    # Only if rctd_reference_build actually ran in THIS invocation — a
    # partial-stages run (e.g. `--stages census`) does not remove the
    # intermediates the next partial invocation would depend on.
    if ran_rctd_reference_build:
        _drop_intermediate_outputs(rundir)

    banner(f"ref-build done in {time.time()-pipeline_t0:.1f}s -> {rundir}/")
    return 0
