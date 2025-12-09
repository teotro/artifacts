#!/bin/bash

# Default values
PARALLEL_RUN=false
NUM_CORES=1  # Default to 1 if not parallel

# --- Argument Parsing ---

# Check for minimum usage: layout_type
if [ -z "$1" ]; then
    echo "Usage: $0 <layout_type> [--parallel <true|false>] [--num_cores <int>]"
    echo "Example: $0 compact --parallel True --num_cores 64"
    exit 1
fi

LAYOUT=$1
shift # Shift past the layout argument

# Parse optional arguments
while [ "$#" -gt 0 ]; do
    case "$1" in
        --parallel)
            if [ "$2" == "True" ] || [ "$2" == "true" ]; then
                PARALLEL_RUN=true
            fi
            shift 2
            ;;
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

if $PARALLEL_RUN; then
    echo "Parallel execution ENABLED with $NUM_CORES cores."
else
    echo "Parallel execution DISABLED. Running sequentially."
fi

# Declare associative array: benchmark → num_qubits
declare -A qubit_lookup=(
    ["random_1"]=29 ["random_2"]=26 ["random_3"]=34 ["random_4"]=56 ["random_5"]=35
    ["random_6"]=42 ["random_7"]=38 ["random_8"]=68 ["random_9"]=12 ["random_10"]=15
    ["random_11"]=32 ["random_12"]=32 ["random_13"]=25 ["random_14"]=40 ["random_15"]=34
    ["random_16"]=34 ["random_17"]=39 ["random_18"]=50 ["random_19"]=56 ["random_20"]=72
)

# --- Function to run a single perturbation ---
# This function is executed by each worker thread/process
run_perturbation() {
    local BENCH_NAME=$1
    local BENCH_DIR=$2
    local LAYOUT=$3
    local NUM_QUBITS=$4
    local OUTPUT_DIR=$5
    local i=$6

    # Specific file names for this perturbation
    local QASM_FILE="$BENCH_DIR/pert_${i}.qasm"
    local PERT="pert_${i}"
    local OUTPUT_FILE="$OUTPUT_DIR/${PERT}.txt"
    # Create a unique temporary error file name for this job
    local ERR_FILE="$OUTPUT_FILE.err" 

    echo "  ▶️  Running $PERT ($BENCH_NAME)..."

    # Execute the python script with a timeout
    timeout 3600 python tracemaker_theo.py \
        --path="$QASM_FILE" \
        --num_qubits="$NUM_QUBITS" \
        --$LAYOUT \
        --benchmark_name="$BENCH_NAME" \
        --benchmark_pert="$PERT" \
        > "$OUTPUT_FILE" 2> "$ERR_FILE"

    local EXIT_CODE=$?

    if [ $EXIT_CODE -ne 0 ]; then
        if [ $EXIT_CODE -eq 124 ]; then
            echo "  ❌ Timeout ($BENCH_NAME - $PERT). Check log: $ERR_FILE"
            echo "Timeout (3600s)" > "$ERR_FILE"
        else
            echo "  ❌ Error in $PERT ($BENCH_NAME): Check log: $ERR_FILE"
            # The error details are already in $ERR_FILE
        fi
    else
        # If successful, remove the temporary error file
        rm -f "$ERR_FILE"
        echo "  ✅ Completed $PERT ($BENCH_NAME)"
    fi
}
export -f run_perturbation # Export the function for use in subshells

# --- Main Logic ---

# Loop through each random benchmark
for j in {18..20}; do
    BENCH_NAME="random_${j}"
    BENCH_DIR="/home/george/artifacts/traceq/synthetic_benchmarks/$BENCH_NAME" # Replace with your own path!
    OUTPUT_DIR="targets/$LAYOUT/$BENCH_NAME"

    NUM_QUBITS=${qubit_lookup[$BENCH_NAME]}
    
    mkdir -p "$OUTPUT_DIR"

    echo "---"
    echo "🔄 Starting $BENCH_NAME with $NUM_QUBITS qubits in layout $LAYOUT..."

    # Array to hold all job commands for this benchmark
    JOB_COMMANDS=()
    
    # Loop through pert_1.qasm to pert_30.qasm
    for i in {1..30}; do
        # Build the command string: Function name followed by its arguments
        COMMAND="run_perturbation $BENCH_NAME $BENCH_DIR $LAYOUT $NUM_QUBITS $OUTPUT_DIR $i"
        JOB_COMMANDS+=("$COMMAND")
    done

    # Run the jobs
    if $PARALLEL_RUN; then
        # Pipe all job commands into 'xargs' to execute them in parallel
        printf "%s\n" "${JOB_COMMANDS[@]}" | xargs -I {} -P "$NUM_CORES" bash -c '{}'
    else
        # Run sequentially
        for COMMAND in "${JOB_COMMANDS[@]}"; do
            bash -c "$COMMAND"
        done
    fi

    echo "✅ Completed $BENCH_NAME"
done

echo "---"
echo "🎉 All benchmarks completed for layout: $LAYOUT"

