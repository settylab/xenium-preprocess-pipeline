"""Stage 5: 10X-style mtx bundle → spacexr Reference .rds (Stage C).

Port of the ref-build-summary Stage C recipe (v3, lines 252-281).
Thin Python wrapper: locates the R script that ships alongside the
package (`ref_build/r/rctd_reference_build.R`), shells out to
`Rscript`, threads through the paths + args, and fails loud on
non-zero exit.

The R side does the actual `spacexr::Reference` call with
`min_UMI = 10, require_int = TRUE`. See the R script for the exact
recipe.

`spacexr` module availability: the default Slurm submit wrapper does
`ml fhR/4.4.1-foss-2023b` which carries Seurat + Matrix + readr, but
NOT spacexr. See `docs/install.md` for the one-time
`~/.claude/r_libs/4.4.1` install. The R script sets .libPaths()
accordingly.
"""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

from ref_build._internal.compat import sentinel_exists
from ref_build._internal.layout import rctd_path
from ref_build._internal.logging import log


def _package_r_script() -> Path:
    """Absolute path to the R script shipped in the package."""
    return Path(__file__).resolve().parent.parent / "r" / "rctd_reference_build.R"


def run_rctd_reference_build(
    sample_id: str,
    mtx_dir: Path,
    output_root: Path,
    rscript_bin: str,
    min_umi: int,
    require_int: bool,
    celltype_col: str,
    label_slash_replacement: str,
    force_rerun: bool,
    run_id: str,
) -> Path:
    """Invoke `Rscript rctd_reference_build.R` on the 10X mtx bundle.
    Writes the .rds DIRECTLY to the final locked-layout path
    `<run_dir>/rctd/<sample_id>_reference.rds`.

    Returns the path to the output RDS.
    """
    out_rds = rctd_path(output_root, sample_id, run_id, "reference")
    out_rds.parent.mkdir(parents=True, exist_ok=True)

    if sentinel_exists(out_rds, force_rerun):
        log(f"[rctd_reference_build] sentinel exists: {out_rds} — skipping "
            f"(pass --force-rerun to re-run).")
        return out_rds

    # Locate Rscript. shutil.which fails loud vs. a silent NotFound
    # deep in subprocess.
    if shutil.which(rscript_bin) is None:
        raise SystemExit(
            f"[rctd_reference_build] Rscript binary not found on PATH: "
            f"{rscript_bin!r}. Options: (a) `module load fhR/4.4.1-foss-2023b` "
            "before invoking the pipeline; (b) pass --rscript-bin "
            "/absolute/path/to/Rscript; (c) set "
            "config.rctd_reference_build.rscript_bin to the explicit binary path."
        )

    r_script = _package_r_script()
    if not r_script.exists():
        raise SystemExit(
            f"[rctd_reference_build] R script missing from the package: {r_script}. "
            "Something's off with the install; expected r/rctd_reference_build.R "
            "under the ref_build package tree."
        )

    args = [
        rscript_bin, "--vanilla",
        str(r_script),
        f"--mtx-dir={mtx_dir}",
        f"--sample-id={sample_id}",
        f"--out-rds={out_rds}",
        f"--min-umi={min_umi}",
        f"--require-int={'TRUE' if require_int else 'FALSE'}",
        f"--celltype-col={celltype_col}",
        f"--label-slash-replacement={label_slash_replacement}",
    ]
    log(f"[rctd_reference_build] launching: {' '.join(args)}")
    # Stream both stdout + stderr into the pipeline log (subprocess
    # inherits our stdout/stderr fds by default). check=True raises on
    # non-zero exit; the resulting CalledProcessError carries the code.
    subprocess.run(args, check=True)

    if not out_rds.exists():
        raise SystemExit(
            f"[rctd_reference_build] Rscript exited 0 but the expected output "
            f"was not written: {out_rds}."
        )
    log(f"[rctd_reference_build] wrote {out_rds}")
    return out_rds
