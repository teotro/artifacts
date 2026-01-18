#!/bin/bash

# ==============================================================================
# TraceQ: Full Pipeline Reproduction Script (HPCA 2026)
# ==============================================================================
# This script automates the end-to-end reproduction of Figures 12-17.
# Total execution time depends on the hardware and number of cores provided.
# ==============================================================================

# Exit on any unhandled error
set -e

# --- Configuration ---
NUM_CORES=${1:-30} # Default to 30 cores if not specified
export TRACEQ_ROOT=$(pwd)
export TRACEQ_RESULTS_DIR="$TRACEQ_ROOT/targets"

echo "======================================================================"
echo "🚀 Starting TraceQ Full Pipeline Reproduction"
echo "Using $NUM_CORES CPU cores"
echo "Results will be stored in: $TRACEQ_RESULTS_DIR"
echo "======================================================================"

# 1. Environment Verification
if [[ "$CONDA_DEFAULT_ENV" != "traceq" ]]; then
    echo "⚠️ Warning: 'traceq' conda environment not active."
    echo "Please run 'conda activate traceq' before starting this script."
    exit 1
fi

# 2. Stage 1: Trace Reconstruction & DAG Generation
echo -e "\n--- [Step 1/4] Trace Reconstruction (Figures 12-15) ---"
LAYOUTS=("compact" "intermediate" "sparse")

for LAYOUT in "${LAYOUTS[@]}"; do
    echo "▶️ Running $LAYOUT layout reconstruction..."
    ./run_tracemaker.sh "$LAYOUT" --parallel True --num_cores "$NUM_CORES"
done

# 3. Stage 2: Rerunning Failed Cases (Handling Timeouts)
echo -e "\n--- [Step 2/4] Resolving Timeouts (Statistical Rigor) ---"
for LAYOUT in "${LAYOUTS[@]}"; do
    echo "▶️ Retrying $LAYOUT failed perturbations (up to 10 attempts)..."
    ./rerun_failed.sh "$LAYOUT" --num_cores "$NUM_CORES"
done

# 4. Stage 3: Subroutine Matching (Figures 16-17)
echo -e "\n--- [Step 3/4] Subroutine Matching (VF3 Algorithm) ---"

echo "▶️ Converting DAGs to .grf format..."
python convert_targets_to_simple_grf.py

# We simulate the progressive timeout increases used in the paper
MATCH_TIMEOUTS=(1 10 60 600 3600)
for LAYOUT in "${LAYOUTS[@]}"; do
    echo "▶️ Matching subroutines for $LAYOUT..."
    ./matching_all.sh "$LAYOUT" 1
    
    # Progressive reruns for remaining cases
    ./rerun_timeouts.sh "$LAYOUT" 10 1
    ./rerun_timeouts.sh "$LAYOUT" 60 10
    ./rerun_timeouts.sh "$LAYOUT" 600 60
    ./rerun_timeouts.sh "$LAYOUT" 3600 600
    
    echo "▶️ Aggregating $LAYOUT matching statistics..."
    ./analyze_layouts.sh "$LAYOUT"
done

# 5. Stage 4: Figure Generation
echo -e "\n--- [Step 4/4] Generating Paper Figures ---"
mkdir -p "$TRACEQ_RESULTS_DIR/figures"

# Figures 12-15
python plot_ambiguity_by_layout.py
python dfs_time_boxplot.py
python generate_success_heatmap.py
python plot_endpoints_vs_dfs_time.py

# Figures 16-17
python analyze_success_rate.py
python plot.py

echo "======================================================================"
echo "✅ Reproduction Complete!"
echo "All generated figures can be found in: $TRACEQ_RESULTS_DIR/figures"
echo "Check .err files in $TRACEQ_RESULTS_DIR for any persistent failures."
echo "======================================================================"