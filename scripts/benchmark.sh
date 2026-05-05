#!/bin/bash

# Author: Martin Saincevic
# Runs the pipeline with multiple clustering thresholds and logs results.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BINARY="${ROOT_DIR}/build/gene_variant_finder"
INPUT_FILE="${ROOT_DIR}/data/example/J29_example.fastq"
OUTPUT_DIR="${ROOT_DIR}/data/example/bench"
REPORT_FILE="${OUTPUT_DIR}/benchmark.tsv"

mkdir -p "${ROOT_DIR}/build" "${OUTPUT_DIR}"

if [[ ! -f "${ROOT_DIR}/build/Makefile" ]]; then
  cmake -S "${ROOT_DIR}" -B "${ROOT_DIR}/build"
fi

cmake --build "${ROOT_DIR}/build"

echo -e "threshold\talleles\tparse_ms\tfilter_ms\ttrim_ms\tcluster_ms\tconsensus_ms\ttotal_ms\tpeak_mb" > "${REPORT_FILE}"

for threshold in 10 12 15 18 20; do
  echo "--- threshold=${threshold} ---"
  run_output="$(${BINARY} --input "${INPUT_FILE}" --output "${OUTPUT_DIR}/alleles_${threshold}.fasta" --cluster-threshold "${threshold}")"

  alleles=$(echo "${run_output}" | awk '/Found [0-9]+ alleles/ {print $2; exit}')
  parse_ms=$(echo "${run_output}" | awk '/Parsing:/ {print $2; exit}')
  filter_ms=$(echo "${run_output}" | awk '/Filtering:/ {print $2; exit}')
  trim_ms=$(echo "${run_output}" | awk '/Trimming:/ {print $2; exit}')
  cluster_ms=$(echo "${run_output}" | awk '/Clustering:/ {print $2; exit}')
  consensus_ms=$(echo "${run_output}" | awk '/Consensus:/ {print $2; exit}')
  total_ms=$(echo "${run_output}" | awk '/Total:/ {print $2; exit}')
  peak_mb=$(echo "${run_output}" | awk '/Peak memory usage:/ {print $4; exit}')

  echo -e "${threshold}\t${alleles}\t${parse_ms}\t${filter_ms}\t${trim_ms}\t${cluster_ms}\t${consensus_ms}\t${total_ms}\t${peak_mb}" >> "${REPORT_FILE}"
  echo ""
done

echo "Benchmark report written to ${REPORT_FILE}"

