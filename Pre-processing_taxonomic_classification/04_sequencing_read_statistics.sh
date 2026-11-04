#!/bin/bash

set -euo pipefail
shopt -s nullglob

seqkit stats -T -a *.fastq.gz |
{
    IFS=$'\t' read -r header
    printf "%s,AvgReadQ\n" "${header//$'\t'/,}"

    while IFS=$'\t' read -r line; do

        file=$(printf '%s\n' "$line" | cut -f1)

        sample=$(basename "$file" .fastq.gz)

        mean_qual=$(
            seqkit fx2tab -q "$file" |
            awk '
            {
                sum += $NF
                n++
            }
            END {
                if (n > 0)
                    printf "%.2f", sum/n
                else
                    print "NA"
            }'
        )

        rest=$(printf '%s\n' "$line" | cut -f2-)

        printf "%s,%s,%s\n" \
            "$sample" \
            "${rest//$'\t'/,}" \
            "$mean_qual"
    done
} > summ_stats.csv
