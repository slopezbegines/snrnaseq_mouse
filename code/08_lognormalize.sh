#!/usr/bin/env bash
# 08_lognormalize.sh — Step 8 (LogNorm path): merge → normalise → integrate → cluster → UMAP
#
# Run from project root:
#   bash code/08_lognormalize.sh
#
# Peak RAM estimate:
#   Merge of N singlet objects (~N × single-sample size, BPCells on disk) +
#   ScaleData for ~2000 HVGs across all cells (in memory, largest spike).
#   Keep MemoryMax at 14 G for a 16 GB machine.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
R_SCRIPT="${SCRIPT_DIR}/08_lognormalize.R"
LOG_DIR="${PROJECT_ROOT}/logs"
LOG_FILE="${LOG_DIR}/08_lognormalize_$(date +%Y%m%d_%H%M%S).log"

cd "$PROJECT_ROOT"
mkdir -p "$LOG_DIR"

echo "=============================================="
echo " Step 08 — LogNormalize Integration Pipeline"
echo " $(date '+%Y-%m-%d %H:%M:%S')"
echo " Project root : $PROJECT_ROOT"
echo " Log file     : $LOG_FILE"
echo "=============================================="

# Option A: systemd-run with 14 GB RAM cap (recommended for 16 GB laptop)
systemd-run \
  --user \
  --scope \
  --quiet \
  -p MemoryMax=14G \
  -p MemorySwapMax=24G \
  -- \
  nice -n 10 \
  ionice -c 3 \
  Rscript --vanilla "$R_SCRIPT" \
  2>&1 | tee "$LOG_FILE"

# Option B: no memory guard (uncomment and comment out Option A above)
# nice -n 10 ionice -c 3 Rscript --vanilla "$R_SCRIPT" 2>&1 | tee "$LOG_FILE"

EXIT_CODE=${PIPESTATUS[0]}

echo "=============================================="
if [ "$EXIT_CODE" -eq 0 ]; then
  echo " FINISHED OK — $(date '+%Y-%m-%d %H:%M:%S')"
  echo " Results:"
  echo "   RPCA   → output/LogNorm/RData/checkpoint_lognorm_02_rpca.rds"
  echo "   Harmony→ output/LogNorm/RData/checkpoint_lognorm_03_harmony.rds"
else
  echo " FAILED (exit $EXIT_CODE) — $(date '+%Y-%m-%d %H:%M:%S')"
  echo " Check log: $LOG_FILE"
fi
echo "=============================================="

exit "$EXIT_CODE"
