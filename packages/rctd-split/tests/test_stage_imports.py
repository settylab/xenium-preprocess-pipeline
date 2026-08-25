"""Smoke test: each stage module imports without ImportError.

Heavy imports (anndata, scanpy, scipy) are LAZY inside `run_<stage>`
bodies, so module-level imports should stay light.
"""
from __future__ import annotations

import importlib

import pytest


STAGES = (
    "rctd_split.stages.rctd_run",
    "rctd_split.stages.split_purify",
    "rctd_split.stages.export_mtx",
    "rctd_split.stages.mtx_to_h5ad",
    "rctd_split.stages.filter_status",
    "rctd_split.stages.postprocess",
    "rctd_split.stages.writeback_to_raw",
    "rctd_split.stages.celltype_writeback",
    "rctd_split.stages.qc_report",
)


@pytest.mark.parametrize("modname", STAGES)
def test_stage_module_imports(modname: str):
    try:
        importlib.import_module(modname)
    except ModuleNotFoundError as e:
        pytest.skip(f"missing dependency for {modname}: {e.name}")


def test_stages_export_run_functions():
    entrypoints = {
        "rctd_split.stages.rctd_run": "run_rctd_run",
        "rctd_split.stages.split_purify": "run_split_purify",
        "rctd_split.stages.export_mtx": "run_export_mtx",
        "rctd_split.stages.mtx_to_h5ad": "run_mtx_to_h5ad",
        "rctd_split.stages.filter_status": "run_filter_status",
        "rctd_split.stages.postprocess": "run_postprocess",
        "rctd_split.stages.writeback_to_raw": "run_writeback_to_raw",
        "rctd_split.stages.celltype_writeback": "run_celltype_writeback",
        "rctd_split.stages.qc_report": "run_qc_report",
    }
    for modname, fname in entrypoints.items():
        try:
            mod = importlib.import_module(modname)
        except ModuleNotFoundError as e:
            pytest.skip(f"missing dependency for {modname}: {e.name}")
        assert hasattr(mod, fname), f"{modname} is missing entry point {fname}"


def test_r_scripts_ship_with_package():
    from rctd_split.stages.rctd_run import _package_r_script as _rctd_run_r
    from rctd_split.stages.split_purify import _package_r_script as _split_r
    from rctd_split.stages.export_mtx import _package_r_script as _export_r

    for get_path, expected_substrs in (
        (_rctd_run_r, ("create.RCTD", "run.RCTD")),
        (_split_r, ("SPLIT::run_post_process_RCTD", "SPLIT::purify")),
        (_export_r, ("writeMM", "features.tsv")),
    ):
        r = get_path()
        assert r.exists(), f"R script missing: {r}"
        txt = r.read_text()
        assert txt.startswith("#!/usr/bin/env Rscript"), f"{r} missing shebang"
        for sub in expected_substrs:
            assert sub in txt, f"{r} missing expected substring: {sub!r}"
