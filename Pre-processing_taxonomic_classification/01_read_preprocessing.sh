#!/bin/bash

# Exit immediately if a command fails, if an unset variable is used, 
# or if any command in a pipeline fails.
set -euo pipefail

# Prevent the loop from running if no matching directories are found.
shopt -s nullglob

# Input/output directories
INPUT_DIR="pass"
OUTPUT_DIR="reads_clean"

# Length filter
MIN_LEN=1300
MAX_LEN=1800

# Trimming
TRIM_5=21
TRIM_3=21

# Process reads
mkdir -p "${OUTPUT_DIR}"

for dir in "${INPUT_DIR}"/barcode[0-9][0-9]; do

    # Extract barcode number
    bname=$(basename "$dir")
    x=${bname#barcode}

    # Collect FASTQ files
    files=("${dir}"/*.fastq.gz)

    if [ ${#files[@]} -eq 0 ]; then
        continue
    fi

    # Define output filename with .gz extension for auto-compression
    output_file="${OUTPUT_DIR}/bc${x}.fastq.gz"

    echo "Processing ${dir} -> ${output_file}"

    # Step 1: low complexity filtering
    # Step 2: length filtering
    # Step 3: trimming

    seqkit grep -vsip \
    "aaaaaaaaaaaaaaaaaaaa,cccccccccccccccccccc,acacacacacacacacacac,agagagagagagagagagag,atatatatatatatatatat,cgcgcgcgcgcgcgcgcgcg" \
    "${files[@]}" |
    seqkit seq -g -m "${MIN_LEN}" -M "${MAX_LEN}" |
    seqkit subseq -r "${TRIM_5}:-${TRIM_3}" \
    -o "${output_file}"

done
