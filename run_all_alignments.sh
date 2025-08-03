#!/bin/bash

# Define input and output folder
INPUT_DIR="."
OUTPUT_DIR="my_alignments"
SCRIPT="python msa_mafft.py"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through all .fasta files
for file in "$INPUT_DIR"/*.fasta; do
    base=$(basename "$file" .fasta)
    outfile="${base}.al"
    echo "Aligning $file..."
    $SCRIPT --infile "$file" --outdir "$OUTPUT_DIR" --outfile "$outfile"
done

