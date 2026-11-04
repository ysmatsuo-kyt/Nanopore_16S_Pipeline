#!/bin/bash

set -euo pipefail
shopt -s nullglob

# Input directory
INPUT_DIR="reads_clean"

# Quality thresholds to evaluate
Q_LIST=(15 16 17 18 19 20)

# Output CSV
OUTPUT_CSV="qual_prev.csv"

header="Sample,Input"

for q in "${Q_LIST[@]}"; do
    header="${header},Q${q}"
done

echo "${header}" > "${OUTPUT_CSV}"

printf "%-8s %10s" "Sample" "Input"

for q in "${Q_LIST[@]}"; do
    printf " %10s" "Q${q}"
done

echo
printf '%*s\n' 90 '' | tr ' ' '-'

for file in "${INPUT_DIR}"/*.fastq.gz; do

    sample=$(basename "${file}" .fastq.gz)

    input_reads=$(
        seqkit stats -T "${file}" 2>/dev/null |
        awk 'NR==2 {print $4}'
    )

    row="${sample},${input_reads}"

    printf "%-8s %10s" \
        "${sample}" \
        "${input_reads}"

    for q in "${Q_LIST[@]}"; do

        n=$(
            seqkit seq -n -Q "${q}" "${file}" |
            wc -l |
            tr -d ' '
        )

        printf " %10s" "${n}"

        row="${row},${n}"
    done

    echo

    echo "${row}" >> "${OUTPUT_CSV}"

done