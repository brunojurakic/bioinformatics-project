#!/bin/bash

# Author: Martin Saincevic
# Simple demo for the example data.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BINARY="${ROOT_DIR}/build/gene_variant_finder"
EXAMPLE_DIR="${ROOT_DIR}/data/example"
OUTPUT_DIR="${ROOT_DIR}/data/example/out"

mkdir -p "${ROOT_DIR}/build" "${OUTPUT_DIR}"

if [[ ! -f "${ROOT_DIR}/build/Makefile" ]]; then
  cmake -S "${ROOT_DIR}" -B "${ROOT_DIR}/build"
fi

cmake --build "${ROOT_DIR}/build"

echo "=== Running demo on bundled example data ==="
echo ""

"${BINARY}" \
   --input "${EXAMPLE_DIR}/J29_example.fastq" \
   --output "${OUTPUT_DIR}/J29_alleles.fasta" \
   --expected "${EXAMPLE_DIR}/expected.fasta" \
   --eval-output "${OUTPUT_DIR}/J29_eval.tsv"

echo ""

"${BINARY}" \
   --input "${EXAMPLE_DIR}/J30_example.fastq" \
   --output "${OUTPUT_DIR}/J30_alleles.fasta" \
   --expected "${EXAMPLE_DIR}/expected.fasta" \
   --eval-output "${OUTPUT_DIR}/J30_eval.tsv"

echo ""
echo "=== Done. Results in ${OUTPUT_DIR} ==="
