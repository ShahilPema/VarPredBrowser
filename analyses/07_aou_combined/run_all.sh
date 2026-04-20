#!/usr/bin/env bash
# Run the full 07_aou_combined pipeline on the workbench.
#
# Usage (typical):
#   mkdir -p output/logs    # ensure log dir exists BEFORE nohup's redirect
#   nohup ./run_all.sh > output/logs/pipeline.log 2>&1 &
#   echo $! > output/logs/pipeline.pid
#
# (The mkdir is also done by run_all.sh once it starts, but the parent shell's
# redirect needs the dir to exist first — there's a race otherwise.)
#
# Monitor:
#   tail -f output/logs/pipeline.log              # stage transitions + summary
#   tail -f output/logs/$(ls -t output/logs/*.stage.log 2>/dev/null | head -1)  # current stage detail
#   ls -lt output/logs/*.stage.log                # newest = currently running
#   cat output/logs/STATUS                        # one-line current state
#
# Options:
#   --skip-pull          skip 00_pull_inputs (use existing inputs/)
#   --skip-build         skip 01_build_per_site (use existing per_site/)
#   --skip-diagnostics   skip 02_diagnostic_plots
#   --no-gp              skip popeve_gp (07) entirely
#   --no-pack            skip 08_pack_outputs (don't push to GCS)
#   --datasets gnomad_only,aou_only,combined   subset (default: all 3)
#   --workers N          override workers passed to fit scripts (default 80)

set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
LOG_DIR="$ROOT/output/logs"
STATUS_FILE="$LOG_DIR/STATUS"
mkdir -p "$LOG_DIR"

DATASETS="gnomad_only,aou_only,combined"
WORKERS=80
SKIP_PULL=0
SKIP_BUILD=0
SKIP_DIAG=0
NO_GP=0
NO_PACK=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --skip-pull) SKIP_PULL=1; shift ;;
        --skip-build) SKIP_BUILD=1; shift ;;
        --skip-diagnostics) SKIP_DIAG=1; shift ;;
        --no-gp) NO_GP=1; shift ;;
        --no-pack) NO_PACK=1; shift ;;
        --datasets) DATASETS="$2"; shift 2 ;;
        --workers) WORKERS="$2"; shift 2 ;;
        *) echo "Unknown arg: $1" >&2; exit 1 ;;
    esac
done

cd "$ROOT"

stamp() { date '+%Y-%m-%d %H:%M:%S'; }

write_status() {
    echo "$(stamp) | $1" > "$STATUS_FILE"
}

# Run a stage. $1 = stage name (used in logs/status), rest = command + args.
stage() {
    local name="$1"; shift
    local log="$LOG_DIR/${name}.stage.log"
    write_status "RUNNING ${name}"
    echo
    echo "[$(stamp)] === ${name} === (detail log: $log)"
    if "$@" >"$log" 2>&1; then
        echo "[$(stamp)] === ${name} done ==="
        write_status "DONE ${name}"
    else
        local rc=$?
        echo "[$(stamp)] === ${name} FAILED (exit $rc, see $log) ===" >&2
        write_status "FAILED ${name} (exit $rc)"
        exit $rc
    fi
}

write_status "STARTED pipeline (datasets=$DATASETS, workers=$WORKERS)"

# Stage 00: pull inputs from GCS
if [[ $SKIP_PULL -eq 0 ]]; then
    stage 00_pull_inputs python3 -u scripts/00_pull_inputs.py
fi

# Stage 01: build per-site parquets (Polars; AoU VAT join on workbench)
if [[ $SKIP_BUILD -eq 0 ]]; then
    stage 01_build_per_site python3 -u scripts/01_build_per_site.py --datasets "$DATASETS"
fi

# Stage 02: diagnostic plots (k summary + sanity scatter + 8-panel mean-E)
if [[ $SKIP_DIAG -eq 0 ]]; then
    stage 02_diagnostic_plots python3 -u scripts/02_diagnostic_plots.py
fi

# Stage 03: cloglog fits (9 fits, ~70 min on 80 vCPU)
stage 03_fit_cloglog python3 -u scripts/03_fit_cloglog.py \
    --datasets "$DATASETS" --workers "$WORKERS"

# Stage 04: quadprog fits (3 fits, ~25 min)
stage 04_fit_quadprog python3 -u scripts/04_fit_quadprog.py \
    --datasets "$DATASETS" --workers "$WORKERS"

# Stage 07: popeve_gp per-dataset so each dataset is restartable if it crashes
# Each dataset is ~7 h on 80 vCPU; 3 datasets = ~20 h total.
if [[ $NO_GP -eq 0 ]]; then
    IFS=',' read -ra DS_ARR <<< "$DATASETS"
    for ds in "${DS_ARR[@]}"; do
        stage "07_fit_popeve_gp_${ds}" python3 -u scripts/07_fit_popeve_gp.py \
            --datasets "$ds" --workers "$WORKERS"
    done
fi

# Stage 08: pack outputs and push to GCS
if [[ $NO_PACK -eq 0 ]]; then
    stage 08_pack_outputs python3 -u scripts/08_pack_outputs.py
fi

echo
echo "[$(stamp)] === pipeline complete ==="
write_status "COMPLETE"
echo "Outputs at gs://zoghbi-lab-data/shahil/07_aou_combined_outputs/"
echo "Pull at BCM: gsutil -m cp gs://zoghbi-lab-data/shahil/07_aou_combined_outputs/*.tar.zst ."
