#!/bin/bash

LAYOUT=$1
TIMEOUT=$2  # e.g., 1, 10, 60, 600
VF3="./bin/vf3"
QUERIES_DIR="/home/george/artifacts/traceq/queries" # Replace this with your own path!
TARGETS_DIR="/home/george/artifacts/traceq/targets/$LAYOUT" # Replace this with your own path!
RESULTS_DIR="results_${TIMEOUT}s"
TODO_DIR="todo"
mkdir -p "$RESULTS_DIR" "$TODO_DIR"

# === Mapping of synthetic benchmarks to queries ===
declare -A BENCH_INSTANCES
BENCH_INSTANCES[random_1]="qft_5 t_npe_4 add_3"
BENCH_INSTANCES[random_2]="qft_6 t_npe_3 add_4"
BENCH_INSTANCES[random_3]="qft_4 t_npe_4 add_5"
BENCH_INSTANCES[random_4]="qft_11 t_npe_5 add_7"
BENCH_INSTANCES[random_5]="t_npe_3 add_9"
BENCH_INSTANCES[random_6]="qft_13 add_10"
BENCH_INSTANCES[random_7]="qft_13 t_npe_5"
BENCH_INSTANCES[random_8]="add_11 t_npe_6"
BENCH_INSTANCES[random_9]="add_3 qft_4"
BENCH_INSTANCES[random_10]="add_4 qft_4"
BENCH_INSTANCES[random_11]="prod_3 outofplace_add_3 qft_4"
BENCH_INSTANCES[random_12]="prod_3 cuccaro_4 qft_4"
BENCH_INSTANCES[random_13]="t_npe_3 cuccaro_7"
BENCH_INSTANCES[random_14]="outofplace_add_5 qft_10 add_5"
BENCH_INSTANCES[random_15]="add_4 outofplace_add_4 cuccaro_4"
BENCH_INSTANCES[random_16]="t_npe_3 qft_5 cuccaro_9"
BENCH_INSTANCES[random_17]="cuccaro_12 qft_4 t_npe_3"
BENCH_INSTANCES[random_18]="add_12 t_npe_3 qft_6"
BENCH_INSTANCES[random_19]="t_npe_4 add_12 qft_5"
BENCH_INSTANCES[random_20]="qft_9 cuccaro_13 adder_12"

# TIMEOUT_FILE="$TODO_DIR/timeouts_${TIMEOUT}s.txt"
TIMEOUT_FILE="$TODO_DIR/${LAYOUT}_timeouts_${TIMEOUT}s.txt"
> "$TIMEOUT_FILE"

for bench in "${!BENCH_INSTANCES[@]}"; do
  queries=(${BENCH_INSTANCES[$bench]})
  log="${RESULTS_DIR}/${bench}_${LAYOUT}.log"
  echo "=== $bench - timeout ${TIMEOUT}s ===" > "$log"

  for i in {1..30}; do
    target="${TARGETS_DIR}/${bench}/pert_${i}.grf"
    if [[ ! -f "$target" ]]; then
      echo "⏩ Missing target: $target" >> "$log"
      continue
    fi

    echo "--- Perturbation: pert_${i} ---" >> "$log"

    for query_name in "${queries[@]}"; do
      query="${QUERIES_DIR}/${query_name}.grf"
      if [[ ! -f "$query" ]]; then
        echo "  ⏩ Missing query: $query" >> "$log"
        continue
      fi

      echo "  Matching: $query_name → pert_${i}" >> "$log"
      result=$(timeout ${TIMEOUT}s "$VF3" "$query" "$target" 2>&1)
      status=$?

      if [[ $status -eq 124 ]]; then
        echo "    ⏱ Timeout (${TIMEOUT}s)" >> "$log"
        echo "$bench pert_${i} $query_name" >> "$TIMEOUT_FILE"
      elif echo "$result" | grep -Eq '(^|\s)-1($|\s)'; then
        echo "    ❌ Not a subgraph" >> "$log"
        echo "$result" | sed 's/^/      /' >> "$log"
      elif [[ $status -eq 0 ]]; then
        echo "    ✅ Match found" >> "$log"
        echo "$result" | sed 's/^/      /' >> "$log"
      else
        echo "    ❌ Runtime error (exit code $status)" >> "$log"
        echo "$result" | sed 's/^/      /' >> "$log"
      fi
    done
    echo "" >> "$log"
  done
done
