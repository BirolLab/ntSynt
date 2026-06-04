#!/bin/bash

#!/usr/bin/env bash
set -euo pipefail

# Usage:
# ./kmc_find_common_kmers.sh genome_list.txt k output_prefix

if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <genome_list.txt> <k> <output_prefix> <tmp dir>"
    exit 1
fi

GENOME_LIST=$(realpath $1)
K="$2"
OUT_PREFIX="$3"
TMP_DIR=$4

CONFIG_FILE="kmc_complex.conf"
DB_LIST=()

BASE_DIR=$(pwd)
mkdir -p "$TMP_DIR" && cd "$TMP_DIR"

echo "[INFO] Counting k-mers..." >&2

# Step 1: Build KMC DBs
while IFS= read -r genome; do
    base=$(basename "$genome")
    db="${base}_kmc"

    echo "[INFO] Processing $genome" >&2

    kmc -v -k"$K" -ci1 -cs255 -t48 -m100 -fm "$genome" "$db" ./

    DB_LIST+=("$db")
done < "$GENOME_LIST"

echo "[INFO] Building KMC complex config..." >&2

# Step 2: Build config file
(
    echo "INPUT:"

    # Assign short aliases (d1, d2, ...)
    for i in "${!DB_LIST[@]}"; do
        idx=$((i + 1))
        echo "d${idx} = ${DB_LIST[$i]}"
    done

    echo ""
    echo "OUTPUT:"

    # Build intersection expression: d1 * d2 * d3 * ...
    expr="d1"
    for ((i=2; i<=${#DB_LIST[@]}; i++)); do
        expr="${expr} * d${i}"
    done

    echo "result = ${expr}"
 ) > "$CONFIG_FILE"

echo "[INFO] Running kmc_tools complex..." >&2

# Step 3: Run intersection
kmc_tools complex "$CONFIG_FILE"

echo "[INFO] Dumping final k-mers..." >&2

# Step 4: Dump result
kmc_tools transform result dump "${OUT_PREFIX}.kmers.txt"
cp "${OUT_PREFIX}.kmers.txt" "${BASE_DIR}"

echo "[INFO] Done." >&2

cd ../