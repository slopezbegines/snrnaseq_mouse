#!/usr/bin/env bash
# 06_integration.sh — Step 6: Dataset integration (RPCA + SCT)
#
# Run from project root:
#   bash code/06_integration.sh
#
# Integration is memory-heavy (all SCT objects loaded at once).
# Recommended: run on a server with >= 32 GB RAM.
# For a local 16 GB laptop attempt, use Option B (systemd-run with swap headroom).
#
# Memory guard options (uncomment ONE block):
#
#   [A] Plain Rscript — default, suitable for servers (no systemd required).
#
#   [B] systemd-run — for a local laptop attempt with a hard RAM cap.
#       May OOM-kill on 16 GB; increase MemorySwapMax if you have swap space.

set -euo pipefail

# ── Paths ──────────────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
R_SCRIPT="${SCRIPT_DIR}/06_integration.R"
LOG_DIR="${PROJECT_ROOT}/logs"
LOG_FILE="${LOG_DIR}/06_integration_$(date +%Y%m%d_%H%M%S).log"

cd "$PROJECT_ROOT"
mkdir -p "$LOG_DIR"

echo "=============================================="
echo " Step 06 — Dataset Integration (RPCA)"
echo " $(date '+%Y-%m-%d %H:%M:%S')"
echo " Project root : $PROJECT_ROOT"
echo " Log file     : $LOG_FILE"
echo "=============================================="

# ── Option A: plain Rscript — active when running in multi-user.target ─────────
nice -n 10 ionice -c 3 Rscript --vanilla "$R_SCRIPT" 2>&1 | tee "$LOG_FILE"

# ── Option B: systemd-run for local laptop (16 GB + swap) ─────────────────────
# Requires an active user D-Bus session (graphical login). Fails after
# 'systemctl isolate multi-user.target' because $DBUS_SESSION_BUS_ADDRESS
# is gone. Use Option A instead when running headless.
#
# systemd-run \
#   --user \
#   --scope \
#   --quiet \
#   -p MemoryMax=14G \
#   -p MemorySwapMax=24G \
#   -- \
#   nice -n 10 \
#   ionice -c 3 \
#   Rscript --vanilla "$R_SCRIPT" \
#   2>&1 | tee "$LOG_FILE"

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
