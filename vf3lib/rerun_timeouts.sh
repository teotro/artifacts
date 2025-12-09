#!/bin/bash

LAYOUT=$1
TIMEOUT=$2
PREV_TIMEOUT=$3

VF3="./bin/vf3"
QUERIES_DIR="/home/george/artifacts/traceq/queries" # Replace this with your own path!
TARGETS_DIR="/home/george/artifacts/traceq/targets/$LAYOUT" # Replace this with your own path!
RESULTS_DIR="results_${TIMEOUT}s"
TODO_DIR="todo"
mkdir -p "$RESULTS_DIR"

# INPUT_TODO="$TODO_DIR/timeouts_${PREV_TIMEOUT}s.txt"
# OUTPUT_TODO="$TODO_DIR/timeouts_${TIMEOUT}s.txt"
INPUT_TODO="$TODO_DIR/${LAYOUT}_timeouts_${PREV_TIMEOUT}s.txt"
OUTPUT_TODO="$TODO_DIR/${LAYOUT}_timeouts_${TIMEOUT}s.txt"
> "$OUTPUT_TODO"

while read -r bench pert query_name; do
  target="${TARGETS_DIR}/${bench}/${pert}.grf"
  query="${QUERIES_DIR}/${query_name}.grf"
  log="${RESULTS_DIR}/${bench}_${LAYOUT}.log"

  if [[ ! -f "$query" || ! -f "$target" ]]; then
    echo "⏩ Skipping missing: $query or $target" >> "$log"
    continue
  fi

  echo "  Matching: $query_name → $pert (timeout ${TIMEOUT}s)" >> "$log"
  result=$(timeout ${TIMEOUT}s "$VF3" "$query" "$target" 2>&1)
  status=$?

  if [[ $status -eq 124 ]]; then
    echo "    ⏱ Timeout (${TIMEOUT}s)" >> "$log"
    echo "$bench $pert $query_name" >> "$OUTPUT_TODO"
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
done < "$INPUT_TODO"
