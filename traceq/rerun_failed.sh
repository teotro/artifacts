#!/bin/bash

# --- CONFIGURATION ---
MAX_ATTEMPTS=10
TIMEOUT_SEC=3600
DEFAULT_CORES=1 # Set conservatively; you can safely raise this to your core count (e.g., 32)
# Ensure NUM_CORES is set based on your available physical cores for maximum safety.

# --- Argument Parsing ---
if [ -z "$1" ]; then
    echo "Usage: $0 <layout_type> [--num_cores <int>]"
    echo "Example: $0 compact --num_cores 32"
    exit 1
fi
LAYOUT=$1
shift

# Parse optional arguments for global core limit
NUM_CORES=$DEFAULT_CORES
while [ "$#" -gt 0 ]; do
    case "$1" in
        --num_cores)
            if [[ "$2" =~ ^[0-9]+$ ]]; then
                NUM_CORES=$2
            fi
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

echo "Starting SAFE Global Parallel Retry script."
echo "Max concurrent failure checks (single Python process each): $NUM_CORES"

# Declare associative array: benchmark → num_qubits
declare -A qubit_lookup=(
    ["random_1"]=29 ["random_2"]=26 ["random_3"]=34 ["random_4"]=56 ["random_5"]=35
    ["random_6"]=42 ["random_7"]=38 ["random_8"]=68 ["random_9"]=12 ["random_10"]=15
    ["random_11"]=32 ["random_12"]=32 ["random_13"]=25 ["random_14"]=40 ["random_15"]=34
    ["random_16"]=34 ["random_17"]=39 ["random_18"]=50 ["random_19"]=56 ["random_20"]=72
)

# --- Core Function: Runs the 10 attempts SEQUENTIALLY for ONE failed perturbation ---
# This function is executed by one of the xargs worker threads.
process_failed_perturbation_sequentially() {
    local BENCH_NAME=$1
    local PERT_NUM=$2
    local LAYOUT=$3
    local NUM_QUBITS=$4
    local MAX_ATTEMPTS=$5
    local TIMEOUT_SEC=$6

    local PERT="pert_${PERT_NUM}"
    local BENCH_DIR="/home/george/artifacts/traceq/synthetic_benchmarks/$BENCH_NAME"
    local OUTPUT_DIR="targets/$LAYOUT/$BENCH_NAME"
    local QASM_FILE="$BENCH_DIR/${PERT}.qasm"
    local OUTPUT_FILE="$OUTPUT_DIR/${PERT}.txt"
    local ERR_FILE="$OUTPUT_DIR/${PERT}.txt.err"

    # Only proceed if the error file still exists
    if [ ! -f "$ERR_FILE" ]; then
        echo "   [SKIP] $BENCH_NAME/$PERT: Error file already removed by another process."
        return 0
    fi
    
    echo "   [START] $BENCH_NAME/$PERT: Retrying sequentially (up to $MAX_ATTEMPTS times)..."
    
    # Loop through attempts SEQUENTIALLY
    for attempt in $(seq 1 $MAX_ATTEMPTS); do
        echo "   🔄 Attempt $attempt..."

        # Re-initialize the error file before running
        > "$ERR_FILE"
        
        # Execute the job
        timeout "$TIMEOUT_SEC" python tracemaker_theo.py \
            --path="$QASM_FILE" \
            --num_qubits="$NUM_QUBITS" \
            --$LAYOUT \
            --benchmark_name="$BENCH_NAME" \
            --benchmark_pert="$PERT" \
            > "$OUTPUT_FILE" 2> "$ERR_FILE"
        
        EXIT_CODE=$?

        if [ $EXIT_CODE -eq 0 ] && [ ! -s "$ERR_FILE" ]; then
            # SUCCESS: Exit code 0 and empty error file
            echo "   ✅ Success on attempt $attempt. Removing .err file."
            rm -f "$ERR_FILE"
            return 0 # Exit the function on success
        elif [ $EXIT_CODE -eq 124 ]; then
            # TIMEOUT: Keep ERR_FILE (which may be empty) and log the timeout
            echo "   ⏰ Timeout after 1 hour on attempt $attempt." >> "$ERR_FILE"
        else
            # PYTHON/OTHER ERROR: Keep ERR_FILE and log the failure
            echo "   ⚠️ Error on attempt $attempt. Check log: $ERR_FILE"
        fi
    done
    
    # If the function reaches here, all attempts failed
    echo "   ❌ Still failing after $MAX_ATTEMPTS attempts. Final status remains FAILED."
    return 1
}
export -f process_failed_perturbation_sequentially

# --- MAIN JOB COLLECTION ---

FAILED_JOBS=()

echo "🔍 Scanning for failed perturbations across all 20 benchmarks..."

# Loop sequentially to identify ALL failed jobs first
for j in {1..20}; do
    BENCH_NAME="random_${j}"
    OUTPUT_DIR="targets/$LAYOUT/$BENCH_NAME"
    NUM_QUBITS=${qubit_lookup[$BENCH_NAME]}

    mkdir -p "$OUTPUT_DIR"

    for i in {1..30}; do
        ERR_FILE="$OUTPUT_DIR/pert_${i}.txt.err"

        if [ -f "$ERR_FILE" ]; then
            # Store job parameters as a single command string
            FAILED_JOBS+=("$BENCH_NAME $i $LAYOUT $NUM_QUBITS $MAX_ATTEMPTS $TIMEOUT_SEC")
        fi
    done
done

TOTAL_FAILED=${#FAILED_JOBS[@]}

if [ $TOTAL_FAILED -eq 0 ]; then
    echo "🎉 No failed perturbations found. Script finished."
    exit 0
fi

echo "Found $TOTAL_FAILED failed jobs. Running them globally in parallel (max $NUM_CORES concurrent checks)..."

# --- GLOBAL PARALLEL EXECUTION ---

# Pipe all commands into 'xargs' to execute them in parallel
printf "%s\n" "${FAILED_JOBS[@]}" | xargs -I {} -P "$NUM_CORES" bash -c 'process_failed_perturbation_sequentially {}'

echo "---"
echo "🎯 Safe Global Parallel Retry script completed for layout: $LAYOUT"