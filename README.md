# Gene Variant Finder (Bioinformatics 1 Project)

Course project for **Bioinformatics 1** (FER): this tool processes sequencing
reads to detect gene variants and supports evaluation on known samples.

[Course website](https://www.fer.unizg.hr/predmet/bio1)

## What this project currently does

- Parses FASTQ reads from deer samples.
- Filters reads by dominant length and trims known adapters.
- Extracts gene region (target length 249 bp).
- Clusters sequences with greedy centroid clustering (Hamming distance).
- Builds one majority-vote consensus sequence per cluster.
- Writes discovered alleles to FASTA output.
- Evaluates discovered alleles against expected FASTA (`--expected`).
- Optionally writes evaluation results to TSV (`--eval-output`).
- Integrates minimap2 mapping and reports read counts by reference
  (`--reference`, `--mapping-output`).

## Requirements

- Linux
- CMake (>= 3.14)
- C++17 compiler
- `minimap2` (optional, only for mapping mode)

## Build

```zsh
cd build
cmake ..
make
```

## Run examples

### 1) Single-sample allele discovery

```zsh
./build/gene_variant_finder --input data/fastq/J29_B_CE_IonXpress_005.fastq --output data/J29_alleles.fasta
```

### 2) Single-sample + evaluation against expected alleles

```zsh
./build/gene_variant_finder --input data/fastq/J29_B_CE_IonXpress_005.fastq --output data/J29_alleles.fasta --expected data/J29B_expected.fasta --eval-output data/J29_eval.tsv
```

### 3) Minimap2 mapping counts by reference

```zsh
./build/gene_variant_finder --input data/fastq/J29_B_CE_IonXpress_005.fastq --reference data/J29B_expected.fasta --mapping-output data/J29_mapping.tsv --minimap-min-mapq 20
```

### 4) Multi-sample mode (all `J*.fastq` in a directory)

```zsh
./build/gene_variant_finder --input-dir data/fastq --output data/all_unique_alleles.fasta
```

## Useful options

- `--cluster-threshold <n>`: max Hamming distance for cluster assignment
  (default `15`)
- `--length-tolerance <n>`: allowed deviation from dominant read length
  (default `5`)
- `--min-cluster-size <n>`: ignore small clusters (default `3`)
- `--merge-threshold <n>`: merge distance across samples in multi-sample mode
  (default `3`)
- `--verbose`: print more detailed output

## Notes

- Minimap2 mapping is currently enabled in `--input` mode.
- In `--input-dir` mode, clustering and cross-sample allele merge are used.
- See `scripts/run_all_samples.sh` for a batch run helper.

