#!/usr/bin/env bash
# Sweep clustrad values on 3I/ATLAS test data to find recovery floor.
# Runs four sequential Docker containers; safe to start before going to bed.
#
# Usage:  bash bin/sweep_3I_clustrad.sh [/abs/path/to/repo]
# Default repo path: directory containing this script's parent.

set -euo pipefail

REPO="${1:-$(cd "$(dirname "$0")/.." && pwd)}"
IMAGE="reductionista/interstellar"
OUT_DIR="${REPO}/test_data/heliolinx/clustrad_sweep"

mkdir -p "$OUT_DIR"

echo "=== 3I clustrad sweep ==="
echo "Repo:   $REPO"
echo "Output: $OUT_DIR"
echo "Image:  $IMAGE"
echo ""

# For each clustrad, set maxrms = 85% of clustrad (matching the empirical rule).
# maxrms is disabled (50M km) above clustrad=1M to preserve the original behaviour
# and let us see the full totRMS distribution.

run_one() {
    local clustrad="$1"
    local maxrms="$2"
    local label="clustrad${clustrad}km"
    local logfile="${OUT_DIR}/${label}.log"

    echo "--- Starting  clustrad=${clustrad} km  maxrms=${maxrms} km ---"
    docker run --rm \
        -v "${REPO}:/app/interstellar" \
        "$IMAGE" \
        python /app/interstellar/python_scripts/run_3I_heliolinx.py \
            --clustrad "$clustrad" \
            --maxrms   "$maxrms" \
            --output-dir "/app/interstellar/test_data/heliolinx/clustrad_sweep" \
        2>&1 | tee "$logfile"
    echo "--- Finished  clustrad=${clustrad} km  (log: $logfile) ---"
    echo ""
}

# Smallest to largest so we fail fast if something is broken.
run_one  100000   85000
run_one  200000  170000
run_one  500000  425000
run_one 1000000  850000

echo "=== Sweep complete. Results in $OUT_DIR ==="
echo ""
echo "Quick summary — hyperbolic candidates per run:"
for csv in "${OUT_DIR}"/3I_refined_clustrad*K.csv; do
    label=$(basename "$csv" .csv)
    count=$(python3 -c "
import pandas as pd, sys
try:
    df = pd.read_csv('$csv')
    if 'orbit_e' in df.columns:
        print(len(df[df['orbit_e'] > 1.0]))
    else:
        print('no orbit_e col')
except Exception as e:
    print(f'error: {e}')
" 2>/dev/null)
    echo "  $label : $count hyperbolic candidates"
done
