#!/usr/bin/env bash
# Run Step 6 (Dataset Integration) outside RStudio to avoid OOM crashes.
# Usage:  bash run_06_integration.sh [--no-log]
#
# Output: output/GSE194315/SCT/RData/checkpoint_06_data_integrated.rds
#         logs/06_integration_<timestamp>.log  (unless --no-log)

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")" && pwd)"
SCRIPT="${REPO_ROOT}/code_claude/06_integration.R"
LOG_DIR="${REPO_ROOT}/logs"
TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
LOG_FILE="${LOG_DIR}/06_integration_${TIMESTAMP}.log"

# R memory flags: expand max vector heap; allow R to grab more RAM than default
R_FLAGS="--no-save --no-restore --quiet"
# Raise the C stack limit (helps with deep Seurat recursion)
export R_MAX_VSIZE="${R_MAX_VSIZE:-32Gb}"

mkdir -p "$LOG_DIR"

if [[ "${1:-}" == "--no-log" ]]; then
  echo "[INFO] Logging disabled — output to stdout only"
  Rscript ${R_FLAGS} "${SCRIPT}"
else
  echo "[INFO] Log: ${LOG_FILE}"
  # tee: show progress live AND write log
  Rscript ${R_FLAGS} "${SCRIPT}" 2>&1 | tee "${LOG_FILE}"
  echo "[INFO] Done. Log saved: ${LOG_FILE}"
fi
