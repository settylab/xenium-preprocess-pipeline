"""Stage 2: RCTD results + spatial test object → unpurified + purified RDS.

Adapted verbatim from
    the internal SPLIT/Proseg workflow
lines 494-523 — the two SPLIT invocations.

The two RDS files are INTERMEDIATE — they live under the run folder's
`intermediate/split/` and are not part of the persisted final output
set (internal issue review).
"""
from __future__ import annotations

import shutil
import os
import subprocess
from pathlib import Path

from rctd_split._internal.compat import sentinel_exists
from rctd_split._internal.layout import intermediate_path
from rctd_split._internal.logging import log


def _package_r_script() -> Path:
    return Path(__file__).resolve().parent.parent / "r" / "split_purify.R"


def run_split_purify(
    sample_id: str,
    run_id: str,
    test_object: Path,
    rctd_results_rds: Path,
    output_root: Path,
    rscript_bin: str,
    assay_name: str,
    DO_purify_singlets: bool,
    force_rerun: bool,
    r_lib_paths: list[str] | list[Path] | None = None,
) -> tuple[Path, Path]:
    """Invoke `Rscript split_purify.R`. Returns (unpurified_rds, purified_rds)."""
    unpurified_rds = intermediate_path(
        output_root, sample_id, run_id, "unpurified_rds",
    )
    purified_rds = intermediate_path(
        output_root, sample_id, run_id, "purified_rds",
    )
    unpurified_rds.parent.mkdir(parents=True, exist_ok=True)

    if sentinel_exists(purified_rds, force_rerun) and unpurified_rds.exists():
        log(f"[split_purify] sentinels exist: {unpurified_rds}, {purified_rds} "
            f"— skipping (pass --force-rerun to re-run).")
        return unpurified_rds, purified_rds

    if shutil.which(rscript_bin) is None:
        raise SystemExit(
            f"[split_purify] Rscript binary not found on PATH: {rscript_bin!r}."
        )

    r_script = _package_r_script()
    if not r_script.exists():
        raise SystemExit(
            f"[split_purify] R script missing from the package: {r_script}."
        )

    args = [
        rscript_bin, "--vanilla",
        str(r_script),
        f"--test-object={test_object}",
        f"--rctd-results-rds={rctd_results_rds}",
        f"--out-unpurified-rds={unpurified_rds}",
        f"--out-purified-rds={purified_rds}",
        f"--assay-name={assay_name}",
        f"--do-purify-singlets={'TRUE' if DO_purify_singlets else 'FALSE'}",
    ]
    env = os.environ.copy()
    if r_lib_paths:
        prepend = ":".join(str(p) for p in r_lib_paths)
        existing = env.get("R_LIBS_USER", "")
        env["R_LIBS_USER"] = f"{prepend}:{existing}" if existing else prepend
        log(f"[split_purify] R_LIBS_USER prepended with: {prepend}")

    log(f"[split_purify] launching: {' '.join(args)}")
    subprocess.run(args, check=True, env=env)

    for label, path in (("unpurified", unpurified_rds), ("purified", purified_rds)):
        if not path.exists():
            raise SystemExit(
                f"[split_purify] Rscript exited 0 but the {label} RDS was not "
                f"written: {path}."
            )
    log(f"[split_purify] wrote {unpurified_rds}, {purified_rds}")
    return unpurified_rds, purified_rds
