#!/usr/bin/env bash
# 05_sct_normalization.sh — Step 5: SCTransform normalization per library
#
# Run from project root:
#   bash code/05_sct_normalization.sh
#
# Processes one library at a time to stay within laptop RAM limits (16 GB).
# Peak usage: checkpoint_04 list (minus nulled slots) + one SCT result (~10-13 GB).
#
# Memory guard options (uncomment ONE block):
#
#   [A] systemd-run — recommended: cgroups v2 enforces a hard RAM cap.
#       Kills R cleanly with OOM if the limit is exceeded.
#
#   [B] No guard — uncomment plain Rscript at the bottom.

set -euo pipefail

# ── Paths ──────────────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
R_SCRIPT="${SCRIPT_DIR}/05_sct_normalization.R"
LOG_DIR="${PROJECT_ROOT}/logs"
LOG_FILE="${LOG_DIR}/05_sct_normalization_$(date +%Y%m%d_%H%M%S).log"

cd "$PROJECT_ROOT"
mkdir -p "$LOG_DIR"

echo "=============================================="
echo " Step 05 — SCTransform Normalization"
echo " $(date '+%Y-%m-%d %H:%M:%S')"
echo " Project root : $PROJECT_ROOT"
echo " Log file     : $LOG_FILE"
echo "=============================================="

# ── Option A: systemd-run with 14 GB RAM cap (recommended for 16 GB laptop) ───
# MemoryMax=14G  → hard RSS ceiling; process is OOM-killed if exceeded
# nice -n 10     → lower CPU priority (keeps desktop responsive)
# ionice -c 3    → idle I/O class (swap pressure won't freeze desktop)
systemd-run \
  --user \
  --scope \
  --quiet \
  -p MemoryMax=14G \
  -p MemorySwapMax=16G \
  -- \
  nice -n 10 \
  ionice -c 3 \
  Rscript --vanilla "$R_SCRIPT" \
  2>&1 | tee "$LOG_FILE"

# ── Option B: no memory guard (uncomment and comment out Option A above) ───────
# nice -n 10 ionice -c 3 Rscript --vanilla "$R_SCRIPT" 2>&1 | tee "$LOG_FILE"

EXIT_CODE=${PIPESTATUS[0]}

echo "=============================================="
if [ "$EXIT_CODE" -eq 0 ]; then
  echo " FINISHED OK — $(date '+%Y-%m-%d %H:%M:%S')"
  echo " Next: bash code/06_integration.sh (server recommended)"
else
  echo " FAILED (exit $EXIT_CODE) — $(date '+%Y-%m-%d %H:%M:%S')"
  echo " Check log: $LOG_FILE"
fi
echo "=============================================="

exit "$EXIT_CODE"
