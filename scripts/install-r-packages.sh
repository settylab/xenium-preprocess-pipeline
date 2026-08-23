#!/usr/bin/env bash
# scripts/install-r-packages.sh — install spacexr + SPLIT into an
# explicit R library, self-verifying that they actually landed there.
#
# `remotes::install_github()` checks the installed SHA across ALL
# `.libPaths()` entries, not just the first — so prepending a fresh
# library directory (the pattern this doc used to show) does not
# protect an install from being silently SKIPPED if a contaminated
# default library ($R_LIBS_USER, e.g.
# ~/R/x86_64-pc-linux-gnu-library/4.4) already has a matching SHA. On
# a shared login node where multiple people's R sessions write to a
# common $R_LIBS_USER default, this is the common case, not the
# exception.
#
# This script:
#   1. Sets R_LIBS_USER explicitly to the target dir (--vanilla implies
#      --no-environ but still honours R_LIBS_USER, since that's a
#      process env var, not something read from an .Renviron file).
#   2. Strips any existing $HOME/R/* entry from the resulting
#      .libPaths() so a contaminated default can't be found at all.
#   3. Passes explicit lib= and force=TRUE to remotes::install_github()
#      so the install cannot be silently skipped by a SHA match
#      anywhere else on .libPaths().
#   4. Asserts afterward that spacexr + SPLIT actually landed in the
#      target dir — the failure mode was a silent skip, so
#      absence-of-install is exactly what has to be checked.
#
# Usage: scripts/install-r-packages.sh --r-lib-dir /abs/path/to/R/library
# Run this AFTER loading the R module (`ml fhR/4.4.1-foss-2023b`) or
# activating an off-cluster R 4.4+ that already has
# Seurat/Matrix/SpatialExperiment installed — see docs/installation.md
# § R side.

set -euo pipefail

R_LIB_DIR=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --r-lib-dir) R_LIB_DIR="$2"; shift 2 ;;
        -h|--help) sed -n '2,20p' "$0"; exit 0 ;;
        *) echo "error: unknown argument: $1" >&2; exit 2 ;;
    esac
done

: "${R_LIB_DIR:?--r-lib-dir is required (absolute path to the target R library)}"

mkdir -p "$R_LIB_DIR"
R_LIB_DIR=$(cd "$R_LIB_DIR" && pwd)

if ! command -v Rscript >/dev/null 2>&1; then
    echo "error: [install-r-packages] Rscript not found on PATH." >&2
    echo "       Load the R module first (ml fhR/4.4.1-foss-2023b) or activate" >&2
    echo "       an off-cluster R 4.4+." >&2
    exit 3
fi

echo "[install-r-packages] target library: $R_LIB_DIR"

R_LIBS_USER="$R_LIB_DIR" Rscript --vanilla -e '
  lib <- Sys.getenv("R_LIBS_USER")
  stopifnot(nzchar(lib))
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)

  # Strip any contaminated default user library so install_github
  # cannot silently match an already-installed SHA there instead of
  # actually installing into `lib`.
  home <- Sys.getenv("HOME")
  paths <- .libPaths()
  if (nzchar(home)) {
    paths <- paths[!grepl(paste0("^", home, "/R($|/)"), paths)]
  }
  if (!(lib %in% paths)) paths <- c(lib, paths)
  .libPaths(paths)
  cat("[install-r-packages] .libPaths():\n")
  print(.libPaths())

  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", lib = lib, repos = "https://cloud.r-project.org")
  }

  remotes::install_github("dmcable/spacexr", lib = lib, force = TRUE)
  remotes::install_github("bdsc-tds/SPLIT",   lib = lib, force = TRUE)

  missing <- Filter(function(pkg) !dir.exists(file.path(lib, pkg)), c("spacexr", "SPLIT"))
  if (length(missing) > 0) {
    cat("ERROR: not installed into", lib, ":", paste(missing, collapse = ", "), "\n")
    quit(status = 1)
  }
  cat("[install-r-packages] verified: spacexr + SPLIT installed in", lib, "\n")
'
echo "[install-r-packages] done."
