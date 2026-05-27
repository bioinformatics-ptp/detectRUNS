#!/bin/bash
# Run validation for each (dataset × type × method × params) as its own R process.
# 64 subprocesses total (4 datasets × 2 types × 2 methods × 4 param sets).
# Each process exits cleanly, resetting memory before the next scan starts.

set -euo pipefail

cd "$(dirname "$0")/.."

DATASETS=(pigData SELMOL Innovagen_HD ADAPTmap)
TYPES=(ROHom ROHet)
METHODS=(sliding consecutive)
PARAMS=(very_lenient lenient strict very_strict)
LOGDIR="dev"
RSCRIPT="Rscript"

N_TOTAL=$(( ${#DATASETS[@]} * ${#TYPES[@]} * ${#METHODS[@]} * ${#PARAMS[@]} ))
N_DONE=0

echo "========================================================"
echo " detectRUNS per-scan validation"
echo " Start: $(date '+%Y-%m-%d %H:%M:%S')"
echo " Total subprocesses: $N_TOTAL"
echo " Order: ${DATASETS[*]}"
echo "========================================================"
echo

for DS in "${DATASETS[@]}"; do
    echo "------------------------------------------------------------"
    echo " [$(date '+%H:%M:%S')] Dataset: $DS"
    echo "------------------------------------------------------------"

    # Clear per-dataset CSVs so this run starts fresh (not appended to a prior run).
    for F in master_summary.csv sanity_checks.csv step_timing.csv; do
        rm -f "Ext_Data/results/$DS/$F"
    done

    for TYPE in "${TYPES[@]}"; do
        for METHOD in "${METHODS[@]}"; do
            for PARAM in "${PARAMS[@]}"; do
                N_DONE=$(( N_DONE + 1 ))
                TAG="${TYPE}_${METHOD}_${PARAM}"
                LOG="${LOGDIR}/validation_${DS}_${TAG}.log"

                echo -n " [$(date '+%H:%M:%S')] ($N_DONE/$N_TOTAL) $DS | $TAG ... "

                DATASET_FILTER="$DS" \
                TYPE_FILTER="$TYPE" \
                METHOD_FILTER="$METHOD" \
                PARAMS_FILTER="$PARAM" \
                    $RSCRIPT "$LOGDIR/run_validation.R" > "$LOG" 2>&1

                EC=$?
                if [ $EC -eq 0 ]; then
                    echo "OK"
                else
                    echo "FAILED (exit $EC) — see $LOG"
                    exit $EC
                fi

            done
        done
    done

    echo " [$(date '+%H:%M:%S')] $DS complete."
    echo
done

echo "------------------------------------------------------------"
echo " [$(date '+%H:%M:%S')] All scans done. Merging summaries..."
echo "------------------------------------------------------------"
$RSCRIPT dev/merge_summaries.R >> dev/validation_master.log 2>&1
echo " [$(date '+%H:%M:%S')] Merge done. Final CSVs in Ext_Data/results/"
echo
echo "========================================================"
echo " ALL DONE: $(date '+%Y-%m-%d %H:%M:%S')"
echo "========================================================"
