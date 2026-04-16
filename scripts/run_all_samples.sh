#!/bin/bash

# Author: Bruno Jurakic
# Runs the cli on all deer FASTQ files.

BINARY="./build/gene_variant_finder"
DATA_DIR="./data/fastq"
OUTPUT_DIR="./data/results"

mkdir -p "$OUTPUT_DIR"

echo "=== Running cli on all deer samples ==="
echo ""

for fastq in "$DATA_DIR"/J*.fastq; do
    sample=$(basename "$fastq" .fastq)
    echo "--- $sample ---"
    "$BINARY" \
        --input "$fastq" \
        --output "$OUTPUT_DIR/${sample}_alleles.fasta"
    echo ""
done

echo "=== Done. Results in $OUTPUT_DIR ==="
