#!/bin/bash

set -euo pipefail
shopt -s nullglob

# Define the output parent directory for EPI2ME analysis
OUT_DIR="epi2me"

# Create the output directory if it does not exist
mkdir -p "$OUT_DIR"

# Process each FASTQ file in the current directory
for f in *.fastq.gz; do
    # Extract the filename without the .fastq.gz extension (e.g., N01.fastq.gz -> N01)
    base="${f%.fastq.gz}"
    
    # Create a subdirectory for each sample
    target_dir="$OUT_DIR/$base"
    mkdir -p "$target_dir"
    
    # Move the FASTQ file into the corresponding subdirectory
    mv -n "$f" "$target_dir/"

    echo "Moved $f to $target_dir/"

done

# Running the EPI2ME wf-16s workflow
nextflow run epi2me-labs/wf-16s \
  --fastq epi2me \
  --taxonomic_rank S \
  --out_dir results
