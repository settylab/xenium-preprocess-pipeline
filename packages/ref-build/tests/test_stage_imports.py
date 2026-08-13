"""Smoke test: each stage module imports without ImportError.

The heavy stage imports (`anndata`, `scanpy`, `scipy`) are LAZILY
imported inside the `run_<stage>` bodies, so this smoke test only
exercises the module-level imports (numpy + pandas + ref_build._internal.*).

Missing science stack? The test is skipped with a clear reason — do
NOT block CI on the full refBuild env being present.
"""
from __future__ import annotations

import importlib

import pytest


STAGES = (
    "ref_build.stages.load_primary_and_donors",
    "ref_build.stages.census",
    "ref_build.stages.assemble",
    "ref_build.stages.export_mtx",
    "ref_build.stages.rctd_reference_build",
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
        "ref_build.stages.load_primary_and_donors": "run_load_primary_and_donors",
        "ref_build.stages.census": "run_census",
        "ref_build.stages.assemble": "run_assemble",
        "ref_build.stages.export_mtx": "run_export_mtx",
        "ref_build.stages.rctd_reference_build": "run_rctd_reference_build",
    }
    for modname, fname in entrypoints.items():
        try:
            mod = importlib.import_module(modname)
        except ModuleNotFoundError as e:
            pytest.skip(f"missing dependency for {modname}: {e.name}")
        assert hasattr(mod, fname), f"{modname} is missing entry point {fname}"


def test_r_script_ships_with_package():
    """The R script that stage 5 shells out to must live alongside the
    package."""
    from ref_build.stages.rctd_reference_build import _package_r_script

    r = _package_r_script()
    assert r.exists(), f"R script missing: {r}"
    # Sanity check: shebang + it's the actual Stage C script (not empty).
    txt = r.read_text()
    assert txt.startswith("#!/usr/bin/env Rscript")
    assert "spacexr" in txt
    assert "Reference(" in txt
    assert "saveRDS" in txt
    # spacexr slash-sanitisation is load-bearing (summary line 262).
    assert 'gsub("/"' in txt or "gsub('/'" in txt
