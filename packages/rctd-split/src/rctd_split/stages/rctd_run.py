"""Stage 1: spatial test object + scRNA reference → RCTD results RDS.

Adapted verbatim from
    the internal SPLIT/Proseg workflow
lines 305-363 — the create.RCTD + run.RCTD calls.

Thin Python wrapper: locates the R script that ships alongside the
package (`rctd_split/r/rctd_run.R`), shells out to `Rscript`, threads
through the resolved paths + args, and fails loud on non-zero exit.
"""
from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

from rctd_split._internal.compat import sentinel_exists
from rctd_split._internal.layout import rctd_path
from rctd_split._internal.logging import log


def _package_r_script() -> Path:
    """Absolute path to the R script shipped in the package."""
    return Path(__file__).resolve().parent.parent / "r" / "rctd_run.R"


def run_rctd_run(
    sample_id: str,
    run_id: str,
    test_object: Path,
    reference_rds: Path,
    output_root: Path,
    rscript_bin: str,
    assay_name: str,
    spatial_reduction: str,
    UMI_min: int,
    counts_MIN: int,
    UMI_min_sigma: int,
    max_cores: int,
    CELL_MIN_INSTANCE: int,
    doublet_mode: str,
    label_slash_replacement: str,
    force_rerun: bool,
    r_lib_paths: list[str] | list[Path] | None = None,
) -> Path:
    """Invoke `Rscript rctd_run.R` and return the path to the RDS."""
    out_rds = rctd_path(output_root, sample_id, run_id, "rctd_results")
    out_rds.parent.mkdir(parents=True, exist_ok=True)

    if sentinel_exists(out_rds, force_rerun):
        log(f"[rctd_run] sentinel exists: {out_rds} — skipping "
            f"(pass --force-rerun to re-run).")
        return out_rds

    if shutil.which(rscript_bin) is None:
        raise SystemExit(
            f"[rctd_run] Rscript binary not found on PATH: {rscript_bin!r}. "
            "Options: (a) `module load fhR/4.4.1-foss-2023b` before invoking "
            "the pipeline; (b) pass --rscript-bin /absolute/path/to/Rscript; "
            "(c) set config.rctd_run.rscript_bin to the explicit binary path."
        )

    r_script = _package_r_script()
    if not r_script.exists():
        raise SystemExit(
            f"[rctd_run] R script missing from the package: {r_script}."
        )

    args = [
        rscript_bin, "--vanilla",
        str(r_script),
        f"--test-object={test_object}",
        f"--reference-rds={reference_rds}",
        f"--out-rds={out_rds}",
        f"--assay-name={assay_name}",
        f"--spatial-reduction={spatial_reduction}",
        f"--umi-min={UMI_min}",
        f"--counts-min={counts_MIN}",
        f"--umi-min-sigma={UMI_min_sigma}",
        f"--max-cores={max_cores}",
        f"--cell-min-instance={CELL_MIN_INSTANCE}",
        f"--doublet-mode={doublet_mode}",
        f"--label-slash-replacement={label_slash_replacement}",
    ]
    env = os.environ.copy()
    if r_lib_paths:
        prepend = ":".join(str(p) for p in r_lib_paths)
        existing = env.get("R_LIBS_USER", "")
        env["R_LIBS_USER"] = f"{prepend}:{existing}" if existing else prepend
        log(f"[rctd_run] R_LIBS_USER prepended with: {prepend}")

    log(f"[rctd_run] launching: {' '.join(args)}")
    subprocess.run(args, check=True, env=env)

    if not out_rds.exists():
        raise SystemExit(
            f"[rctd_run] Rscript exited 0 but the expected output "
            f"was not written: {out_rds}."
        )
    log(f"[rctd_run] wrote {out_rds}")
    return out_rds
