#!/bin/bash -l
# Small end-to-end verify of the R_LIBS_USER propagation fix.
# Simulates Tracy's env by setting R_LIBS_USER to a path that does NOT
# contain spacexr, so only the .libPaths() prepend in the R script can
# make spacexr reachable — and only the fix's Sys.setenv can make it
# reachable to PSOCK workers.
#
#SBATCH --job-name=rctd-verify-parallel
#SBATCH --partition=campus-new
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:20:00

set -euo pipefail

LOG_DIR="/fh/fast/setty_m/user/ryang/nexus/nexus-code/scratch/rctd-split/tests/.verify-logs"
mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/verify_${SLURM_JOB_ID:-local-$(date +%s)}.log"
exec >"$LOG_FILE" 2>&1

echo "===== verify_parallel_fix job ====="
echo "JOBID:    ${SLURM_JOB_ID:-NA}"
echo "NODELIST: ${SLURM_NODELIST:-NA}"
echo "DATE:     $(date)"
echo "LOG_FILE: $LOG_FILE"

# Load fhR module the same way rctd-split's real submit does.
source /app/lmod/lmod/init/bash 2>/dev/null || true
ml fhR/4.4.1-foss-2023b
echo "which Rscript: $(command -v Rscript)"

# Simulate Tracy's env: R_LIBS_USER is a directory that does NOT hold spacexr.
FAKE_LIB=/tmp/rctd_verify_fake_lib_$$
mkdir -p "$FAKE_LIB"
export R_LIBS_USER="$FAKE_LIB"
echo "R_LIBS_USER (pre-Rscript) = $R_LIBS_USER"

# Run the verify script. Uses --vanilla mirroring rctd_run.R's invocation.
Rscript --vanilla /fh/fast/setty_m/user/ryang/nexus/nexus-code/scratch/rctd-split/tests/verify_parallel_fix.R

echo "===== verify DONE at $(date) ====="
touch "$LOG_DIR/verify_${SLURM_JOB_ID:-local}.sentinel"
