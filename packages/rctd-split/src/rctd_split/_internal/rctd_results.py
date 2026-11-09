"""Shared helper: extract per-cell RCTD results from a spacexr RCTD
.rds file into a `pandas.DataFrame` indexed by cell_id.

Used by both the pipeline-integrated `celltype_writeback` stage (v0.2.9)
and the standalone `scripts/hotfix_rctd_writeback_to_xenium_ranger.py`.

The R-side extraction lives in ``r/rctd_results_to_csv.R`` so both call
sites share exactly one parser of the RCTD `@results$results_df` shape
(including the SPLIT::run_post_process_RCTD columns).
"""
from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from pathlib import Path

from rctd_split._internal.logging import log


def _r_script_path() -> Path:
    """Path to the bundled R extractor script."""
    return Path(__file__).resolve().parent.parent / "r" / "rctd_results_to_csv.R"


def extract_rctd_results_df(
    *,
    rctd_results_rds: Path,
    rscript_bin: str = "Rscript",
    r_lib_paths: list[str] | list[Path] | None = None,
    tmp_csv: Path | None = None,
):
    """Shell out to `Rscript rctd_results_to_csv.R`, read the CSV,
    return a ``pandas.DataFrame`` indexed by `cell_id`.

    Parameters
    ----------
    rctd_results_rds : Path
        The spacexr RCTD .rds file (as emitted by `rctd_run`).
    rscript_bin : str
        `Rscript` binary name/path (default `"Rscript"`).
    r_lib_paths : list of path-like, optional
        Prepended to `R_LIBS_USER` for the Rscript subprocess, same
        convention as the other R-side stages.
    tmp_csv : Path, optional
        Where to write the CSV. If None, a tempfile is used and
        cleaned up after the DataFrame is read.

    Returns
    -------
    pandas.DataFrame
        One row per cell; `index.name == "cell_id"`. Columns are the
        subset of
        ``{spot_class, first_type, second_type, first_class,
        second_class, min_score, singlet_score, conv_all,
        conv_doublet, purification_status, w1_larger_w2, same_class}``
        actually present in `results_df`.
    """
    import pandas as pd

    rctd_results_rds = Path(rctd_results_rds)
    if not rctd_results_rds.exists():
        raise SystemExit(
            f"[rctd_results] rctd_results_rds not found: {rctd_results_rds}"
        )

    if shutil.which(rscript_bin) is None:
        raise SystemExit(
            f"[rctd_results] Rscript binary not found on PATH: "
            f"{rscript_bin!r}. Load the R module (e.g. "
            "`module load R/4.4.0-gfbf-2023b` on Gizmo) or pass an "
            "explicit path."
        )

    r_script = _r_script_path()
    if not r_script.exists():
        raise SystemExit(
            f"[rctd_results] R extractor script missing from the package: "
            f"{r_script}."
        )

    cleanup_tmp = False
    if tmp_csv is None:
        fd, tmp_name = tempfile.mkstemp(
            prefix="rctd_results_", suffix=".csv",
        )
        os.close(fd)
        tmp_csv = Path(tmp_name)
        cleanup_tmp = True
    else:
        tmp_csv = Path(tmp_csv)
        tmp_csv.parent.mkdir(parents=True, exist_ok=True)

    args = [
        rscript_bin, "--vanilla",
        str(r_script),
        f"--rctd-results-rds={rctd_results_rds}",
        f"--out-csv={tmp_csv}",
    ]
    env = os.environ.copy()
    if r_lib_paths:
        prepend = ":".join(str(p) for p in r_lib_paths)
        args.append(f"--r-lib-paths={prepend}")
        existing = env.get("R_LIBS_USER", "")
        env["R_LIBS_USER"] = f"{prepend}:{existing}" if existing else prepend

    log(f"[rctd_results] launching: {' '.join(args)}")
    try:
        subprocess.run(args, check=True, env=env)
        if not tmp_csv.exists():
            raise SystemExit(
                f"[rctd_results] Rscript exited 0 but no CSV was written: "
                f"{tmp_csv}."
            )
        df = pd.read_csv(tmp_csv, dtype={"cell_id": str})
        df = df.set_index("cell_id")
        log(f"[rctd_results] loaded {len(df)} rows × {df.shape[1]} cols "
            f"({', '.join(df.columns)})")
        return df
    finally:
        if cleanup_tmp and tmp_csv.exists():
            tmp_csv.unlink()
