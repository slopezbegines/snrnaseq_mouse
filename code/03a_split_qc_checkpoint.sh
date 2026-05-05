#!/usr/bin/env bash
# 03a_split_qc_checkpoint.sh — One-time split of checkpoint_02 into per-library files
#
# Run from project root:
#   bash code/03a_split_qc_checkpoint.sh
#
# Uses 22 GB MemoryMax to accommodate the full checkpoint_02 in RAM (~12-15 GB).
# This only needs to run ONCE. After it completes, use 03_filter_cells.sh normally.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
R_SCRIPT="${SCRIPT_DIR}/03a_split_qc_checkpoint.R"
LOG_DIR="${PROJECT_ROOT}/logs"
LOG_FILE="${LOG_DIR}/03a_split_qc_$(date +%Y%m%d_%H%M%S).log"

cd "$PROJECT_ROOT"
mkdir -p "$LOG_DIR"

echo "=============================================="
echo " Step 03a — Split QC Checkpoint (one-time)"
echo " $(date '+%Y-%m-%d %H:%M:%S')"
echo " Project root : $PROJECT_ROOT"
echo " Log file     : $LOG_FILE"
echo "=============================================="

systemd-run \
  --user \
  --scope \
  --quiet \
  -p MemoryMax=22G \
  -p MemorySwapMax=4G \
  -- \
  nice -n 10 \
  ionice -c 3 \
  Rscript --vanilla "$R_SCRIPT" \
  2>&1 | tee "$LOG_FILE"

EXIT_CODE=${PIPESTATUS[0]}

echo "=============================================="
if [ "$EXIT_CODE" -eq 0 ]; then
  echo " FINISHED OK — $(date '+%Y-%m-%d %H:%M:%S')"
  echo " Next: bash code/03_filter_cells.sh"
else
  echo " FAILED (exit $EXIT_CODE) — $(date '+%Y-%m-%d %H:%M:%S')"
  echo " Check log: $LOG_FILE"
fi
echo "=============================================="

exit "$EXIT_CODE"
