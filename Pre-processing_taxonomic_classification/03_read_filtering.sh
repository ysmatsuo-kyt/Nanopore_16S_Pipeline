#!/bin/bash

set -euo pipefail
shopt -s nullglob

# Input/output directories
INPUT_DIR="reads_clean"
OUTPUT_DIR="reads_filt"

# User-defined quality cutoff
Q_CUTOFF=17

mkdir -p "${OUTPUT_DIR}"

echo "Quality cutoff: Q${Q_CUTOFF}"
echo

for file in "${INPUT_DIR}"/*.fastq.gz; do

    sample=$(basename "${file}")

    output_file="${OUTPUT_DIR}/${sample}"

    echo "Processing ${file} -> ${output_file}"

    seqkit seq \
        -Q "${Q_CUTOFF}" \
        "${file}" \
        -o "${output_file}"

done
