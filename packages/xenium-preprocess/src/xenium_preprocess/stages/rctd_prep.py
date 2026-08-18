"""Stage 4: SPLIT triple → RCTD test object (RDS).

Adapted from
    the internal spatial-RCTD preparation reference

The reference Rmd reads the mtx/features/barcodes triple + the two
sidecars written by stage 3, builds a SpatialExperiment, converts it
to a Seurat object with a "Proseg" assay + a spatial DimReduc keyed
"ST_", and saves as an RDS. That RDS is the RCTD "test object" —
downstream `create.RCTD(spatial_seurat, reference)` calls take this
as the query side of the deconvolution.

The Rmd's second half (building a spacexr `Reference` from a scRNA
mtx triple) is OUT OF SCOPE for xenium-preprocess — the scRNA reference is user
data, not derived from proseg. That part will land in a later step.

Which anndata layer feeds RCTD: `rctd_prep.source_layer` in the config
(default `maxpost_counts` — user request 2026-07-10, (internal issue review).
The choice is applied UPSTREAM at split_prep (pipeline.py threads
`rctd_prep.source_layer` into `run_split_prep(layer=...)`), so this
module itself does not touch the h5ad — it just consumes the mtx
already exported by split_prep. RCTD's `create.RCTD(...,
require_int=TRUE)` requires integer counts, which is why maxpost is
the default (the alternative `expected_counts` layer is continuous).

This module is a thin Python wrapper: it locates the R script that
ships alongside the package (`xenium_preprocess/r/rctd_prep.R`),
shells out to `Rscript`, threads through the resolved paths + args,
and fails loud on non-zero exit.
"""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

from xenium_preprocess._internal.compat import sentinel_exists
from xenium_preprocess._internal.layout import rctd_path
from xenium_preprocess._internal.logging import log


def _package_r_script() -> Path:
    """Absolute path to the R script shipped in the package."""
    return Path(__file__).resolve().parent.parent / "r" / "rctd_prep.R"


def run_rctd_prep(
    sample_id: str,
    run_id: str,
    split_dir: Path,
    output_root: Path,
    name_suffix: str,
    rscript_bin: str,
    assay_name: str,
    spatial_key: str,
    force_rerun: bool,
) -> Path:
    """Invoke `Rscript rctd_prep.R` on the SPLIT triple.

    Returns the path to the output RDS.
    """
    # Final locked layout: rctd/<S>_test_object.rds under the run folder.
    out_rds = rctd_path(output_root, sample_id, run_id, "test_object")
    out_rds.parent.mkdir(parents=True, exist_ok=True)

    if sentinel_exists(out_rds, force_rerun):
        log(f"[rctd_prep] sentinel exists: {out_rds} — skipping "
            f"(pass --force-rerun to re-run).")
        return out_rds

    # Locate Rscript. shutil.which fails loud vs. a silent NotFound
    # deep in subprocess.
    if shutil.which(rscript_bin) is None:
        raise SystemExit(
            f"[rctd_prep] Rscript binary not found on PATH: {rscript_bin!r}. "
            "Options: (a) `module load fhR/4.4.1-foss-2023b` before invoking "
            "the pipeline; (b) pass --rscript-bin /absolute/path/to/Rscript; "
            "(c) set config.rctd_prep.rscript_bin to the explicit binary path."
        )

    r_script = _package_r_script()
    if not r_script.exists():
        raise SystemExit(
            f"[rctd_prep] R script missing from the package: {r_script}. "
            "Something's off with the install; expected r/rctd_prep.R "
            "under the xenium_preprocess package tree."
        )

    stem = f"{sample_id}{name_suffix}"
    args = [
        rscript_bin, "--vanilla",
        str(r_script),
        f"--split-dir={split_dir}",
        f"--stem={stem}",
        f"--out-rds={out_rds}",
        f"--assay-name={assay_name}",
        f"--spatial-key={spatial_key}",
    ]
    log(f"[rctd_prep] launching: {' '.join(args)}")
    # Stream both stdout + stderr into the pipeline log (subprocess
    # inherits our stdout/stderr fds by default). check=True raises on
    # non-zero exit; the resulting CalledProcessError carries the code.
    subprocess.run(args, check=True)

    if not out_rds.exists():
        raise SystemExit(
            f"[rctd_prep] Rscript exited 0 but the expected output "
            f"was not written: {out_rds}."
        )
    log(f"[rctd_prep] wrote {out_rds}")
    return out_rds
