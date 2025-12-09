#!/bin/bash

# === INPUT CHECK ===
if [[ $# -ne 1 ]]; then
  echo "Usage: $0 <layout>   # e.g., compact, sparse, intermediate"
  exit 1
fi

LAYOUT=$1
TIMEOUTS=(1 10 60 600 3600)
BENCHMARKS=(random_{1..20})
OUTDIR="layout_results"
mkdir -p "$OUTDIR"
OUTFILE="$OUTDIR/${LAYOUT}.txt"
> "$OUTFILE"  # Clear existing

for BENCH in "${BENCHMARKS[@]}"; do
  declare -A SUCC=()    # query → count of success
  declare -A FAIL=()    # query → count of fail
  declare -A TIME=()    # query → count of timeout
  declare -A TOTAL=()   # query → total attempts
  declare -A SEEN=()    # query_pert → avoid duplicate recording

  for TO in "${TIMEOUTS[@]}"; do
    LOG="results_${TO}s/${BENCH}_${LAYOUT}.log"
    [[ ! -f "$LOG" ]] && continue

    PERT=""
    QUERY=""
    CURRENT_KEY=""

    while IFS= read -r line; do
      line=$(echo "$line" | sed 's/^[[:space:]]*//')

      if [[ "$line" =~ ^---\ Perturbation:\ pert_([0-9]+) ]]; then
        PERT="pert_${BASH_REMATCH[1]}"
      elif [[ "$line" =~ ^Matching:\ ([a-zA-Z0-9_]+)[[:space:]]*→[[:space:]]*pert_ ]]; then
        QUERY="${BASH_REMATCH[1]}"
        CURRENT_KEY="${QUERY}_${PERT}"
      elif [[ "$line" =~ ✅\ Match\ found ]]; then
        if [[ -z "${SEEN[$CURRENT_KEY]}" ]]; then
          SUCC["$QUERY"]=$((SUCC["$QUERY"] + 1))
          TOTAL["$QUERY"]=$((TOTAL["$QUERY"] + 1))
          SEEN["$CURRENT_KEY"]=1
        fi
      elif [[ "$line" =~ ❌\ Not\ a\ subgraph ]]; then
        if [[ -z "${SEEN[$CURRENT_KEY]}" ]]; then
          FAIL["$QUERY"]=$((FAIL["$QUERY"] + 1))
          TOTAL["$QUERY"]=$((TOTAL["$QUERY"] + 1))
          SEEN["$CURRENT_KEY"]=1
        fi
      elif [[ "$line" =~ ⏱\ Timeout ]]; then
        if [[ -z "${SEEN[$CURRENT_KEY]}" ]]; then
          TIME["$QUERY"]=$((TIME["$QUERY"] + 1))
          TOTAL["$QUERY"]=$((TOTAL["$QUERY"] + 1))
          SEEN["$CURRENT_KEY"]=1
        fi
      fi
    done < "$LOG"
  done

  # Print summary per benchmark
  if [[ ${#TOTAL[@]} -gt 0 ]]; then
    queries=("${!TOTAL[@]}")
    IFS=$'\n' sorted=($(sort <<<"${queries[*]}"))
    
    echo -n "$BENCH consists of "
    (
      first=1
      for q in "${sorted[@]}"; do
        [[ $first -eq 0 ]] && echo -n " + "
        echo -n "$q"
        first=0
      done
    )
    echo ","
    
    echo -n "I found: "
    (
      first=1
      for q in "${sorted[@]}"; do
        [[ $first -eq 0 ]] && echo -n " + "
        echo -n "${SUCC[$q]:-0}/${TOTAL[$q]} $q"
        first=0
      done
    )
    echo

    echo -n "Failed:  "
    (
      first=1
      for q in "${sorted[@]}"; do
        [[ $first -eq 0 ]] && echo -n " + "
        echo -n "${FAIL[$q]:-0}/${TOTAL[$q]} $q"
        first=0
      done
    )
    echo

    echo -n "Timeout: "
    (
      first=1
      for q in "${sorted[@]}"; do
        [[ $first -eq 0 ]] && echo -n " + "
        echo -n "${TIME[$q]:-0}/${TOTAL[$q]} $q"
        first=0
      done
    )
    echo -e "\n"

    # === Save to .txt file ===
    for q in "${sorted[@]}"; do
      echo "$BENCH $q ${SUCC[$q]:-0} ${FAIL[$q]:-0} ${TIME[$q]:-0} ${TOTAL[$q]:-0}" >> "$OUTFILE"
    done
  fi
done
