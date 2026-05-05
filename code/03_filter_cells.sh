#!/usr/bin/env bash
# 03_filter_cells.sh — Step 3: Cell filtering (Demuxlet singlets + QC thresholds)
#
# Run from project root:
#   bash code_claude/03_filter_cells.sh
#
# Memory guard options (uncomment ONE block):
#
#   [A] systemd-run — recommended: cgroups v2 enforces a hard RAM cap.
#       Kills R cleanly with OOM if the limit is exceeded.
#       Requires: systemd user session (already confirmed present).
#
#   [B] No guard — safe for this step (filtering is lightweight).
#       Uncomment the plain Rscript line at the bottom.

set -euo pipefail

# ── Paths ──────────────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
R_SCRIPT="${SCRIPT_DIR}/03_filter_cells.R"
LOG_DIR="${PROJECT_ROOT}/logs"
LOG_FILE="${LOG_DIR}/03_filter_cells_$(date +%Y%m%d_%H%M%S).log"

cd "$PROJECT_ROOT"
mkdir -p "$LOG_DIR"

echo "=============================================="
echo " Step 03 — Cell Filtering"
echo " $(date '+%Y-%m-%d %H:%M:%S')"
echo " Project root : $PROJECT_ROOT"
echo " Log file     : $LOG_FILE"
echo "=============================================="

# ── Option A: systemd-run with 13 GB RAM cap (recommended) ────────────────────
# MemoryMax=13G  → hard RSS ceiling via cgroups v2; process is OOM-killed if hit
# nice -n 10     → slightly lower CPU priority (keeps desktop responsive)
# ionice -c 3    → idle I/O class (swap pressure won't freeze the desktop)
systemd-run \
  --user \
  --scope \
  --quiet \
  -p MemoryMax=13G \
  -p MemorySwapMax=2G \
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
else
  echo " FAILED (exit $EXIT_CODE) — $(date '+%Y-%m-%d %H:%M:%S')"
  echo " Check log: $LOG_FILE"
fi
echo "=============================================="

exit "$EXIT_CODE"
