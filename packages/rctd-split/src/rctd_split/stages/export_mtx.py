"""Stage 3: unpurified + purified Seurat RDS → 10X-style mtx bundles.

Thin Python wrapper: locates the R script that ships alongside the
package (`rctd_split/r/export_mtx.R`), shells out to `Rscript`, threads
through the resolved paths + args, and fails loud on non-zero exit.

Both bundles are INTERMEDIATE — they live under `intermediate/mtx/`
and are not part of the persisted final output set (dropped per user request
on (internal issue review).
"""
from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

from rctd_split._internal.compat import sentinel_exists
from rctd_split._internal.layout import mtx_bundle_dir
from rctd_split._internal.logging import log


def _package_r_script() -> Path:
    return Path(__file__).resolve().parent.parent / "r" / "export_mtx.R"


def run_export_mtx(
    sample_id: str,
    run_id: str,
    unpurified_rds: Path,
    purified_rds: Path,
    output_root: Path,
    rscript_bin: str,
    assay_name: str,
    gzip_outputs: bool,
    force_rerun: bool,
    r_lib_paths: list[str] | list[Path] | None = None,
) -> tuple[Path, Path]:
    """Invoke `Rscript export_mtx.R`. Returns (unpurified_dir, purified_dir)."""
    unpurified_dir = mtx_bundle_dir(output_root, sample_id, run_id, "unpurified")
    purified_dir = mtx_bundle_dir(output_root, sample_id, run_id, "purified")
    unpurified_dir.mkdir(parents=True, exist_ok=True)
    purified_dir.mkdir(parents=True, exist_ok=True)

    gz = ".gz" if gzip_outputs else ""
    sentinel = purified_dir / f"{sample_id}_purified_counts.mtx{gz}"
    if sentinel_exists(sentinel, force_rerun):
        log(f"[export_mtx] sentinel exists: {sentinel} — skipping "
            f"(pass --force-rerun to re-run).")
        return unpurified_dir, purified_dir

    if shutil.which(rscript_bin) is None:
        raise SystemExit(
            f"[export_mtx] Rscript binary not found on PATH: {rscript_bin!r}."
        )

    r_script = _package_r_script()
    if not r_script.exists():
        raise SystemExit(
            f"[export_mtx] R script missing from the package: {r_script}."
        )

    args = [
        rscript_bin, "--vanilla",
        str(r_script),
        f"--sample-id={sample_id}",
        f"--unpurified-rds={unpurified_rds}",
        f"--purified-rds={purified_rds}",
        f"--unpurified-out-dir={unpurified_dir}",
        f"--purified-out-dir={purified_dir}",
        f"--assay-name={assay_name}",
        f"--gzip-outputs={'TRUE' if gzip_outputs else 'FALSE'}",
    ]
    env = os.environ.copy()
    if r_lib_paths:
        prepend = ":".join(str(p) for p in r_lib_paths)
        existing = env.get("R_LIBS_USER", "")
        env["R_LIBS_USER"] = f"{prepend}:{existing}" if existing else prepend
        log(f"[export_mtx] R_LIBS_USER prepended with: {prepend}")

    log(f"[export_mtx] launching: {' '.join(args)}")
    subprocess.run(args, check=True, env=env)

    if not sentinel.exists():
        raise SystemExit(
            f"[export_mtx] Rscript exited 0 but the expected sentinel "
            f"was not written: {sentinel}."
        )
    log(f"[export_mtx] wrote bundles: {unpurified_dir}, {purified_dir}")
    return unpurified_dir, purified_dir
