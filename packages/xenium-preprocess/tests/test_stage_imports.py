"""Smoke test: each stage module imports without ImportError.

The heavy stage imports (`anndata`, `scanpy`, `umap`, `leidenalg`,
`scipy`) are LAZILY imported inside the `run_<stage>` bodies, so this
smoke test only exercises the module-level imports (numpy + pandas +
xenium_preprocess._internal.*).

Missing science stack? The test is skipped with a clear reason — do
NOT block CI on the full xeniumPreprocess env being present.
"""
from __future__ import annotations

import importlib

import pytest


STAGES = (
    "xenium_preprocess.stages.proseg_to_anndata",
    "xenium_preprocess.stages.enrich_xenium_id",
    "xenium_preprocess.stages.qc_filter",
    "xenium_preprocess.stages.xenium_ranger_to_anndata",
    "xenium_preprocess.stages.preprocess",
    "xenium_preprocess.stages.split_prep",
    "xenium_preprocess.stages.rctd_prep",
)


@pytest.mark.parametrize("modname", STAGES)
def test_stage_module_imports(modname: str):
    try:
        importlib.import_module(modname)
    except ModuleNotFoundError as e:
        pytest.skip(f"missing dependency for {modname}: {e.name}")


def test_stages_export_run_functions():
    """Each stage exposes a `run_<stage>` entry point that pipeline.py
    calls into."""
    entrypoints = {
        "xenium_preprocess.stages.proseg_to_anndata": "run_proseg_to_anndata",
        "xenium_preprocess.stages.enrich_xenium_id": "run_enrich_xenium_id",
        "xenium_preprocess.stages.qc_filter": "run_qc_filter",
        "xenium_preprocess.stages.xenium_ranger_to_anndata": "run_xenium_ranger_to_anndata",
        "xenium_preprocess.stages.preprocess": "run_preprocess",
        "xenium_preprocess.stages.split_prep": "run_split_prep",
        "xenium_preprocess.stages.rctd_prep": "run_rctd_prep",
    }
    for modname, fname in entrypoints.items():
        try:
            mod = importlib.import_module(modname)
        except ModuleNotFoundError as e:
            pytest.skip(f"missing dependency for {modname}: {e.name}")
        assert hasattr(mod, fname), f"{modname} is missing entry point {fname}"


def test_r_script_ships_with_package():
    """The R script that stage 4 shells out to must live alongside the
    package."""
    from xenium_preprocess.stages.rctd_prep import _package_r_script

    r = _package_r_script()
    assert r.exists(), f"R script missing: {r}"
    # Sanity check: shebang line + it's the actual RCTD script (not empty).
    txt = r.read_text()
    assert txt.startswith("#!/usr/bin/env Rscript")
    assert "SpatialExperiment" in txt
    assert "saveRDS" in txt
