// Author: Bruno Jurakic

#ifndef SEQUENCE_FILTER_H
#define SEQUENCE_FILTER_H

#include <map>
#include <vector>

#include "fastq_parser.h"

// Counts how many reads have each length.
std::map<int, int> ComputeLengthHistogram(const std::vector<FastqRead>& reads);

// Returns the length that appears most frequently.
int FindMostCommonLength(const std::map<int, int>& histogram);

// Keeps only reads whose length is within [target - tolerance, target +
// tolerance].
std::vector<FastqRead> FilterByLength(const std::vector<FastqRead>& reads,
                                      int target_length, int tolerance);

// Prints a summary of the length distribution to stdout.
void PrintLengthStats(const std::map<int, int>& histogram);

// Trims adapter sequences from a read, returning only the gene region.
// Returns an empty string if the prefix adapter is not found.
std::string TrimAdapters(const std::string& sequence,
                         const std::string& prefix_adapter,
                         const std::string& suffix_adapter);

// Trims adapters from all reads and returns only the ones that produce
// exactly target_gene_length bases after trimming.
std::vector<std::string> TrimAndExtractGenes(
    const std::vector<FastqRead>& reads, const std::string& prefix_adapter,
    const std::string& suffix_adapter, int target_gene_length);

// A trimmed gene sequence with the name of its original FASTQ read.
struct NamedGene {
  std::string read_name;
  std::string sequence;
};

// Trims adapters and keeps the original read name alongside each sequence.
// Only includes reads that produce exactly target_gene_length bases.
std::vector<NamedGene> TrimAndExtractNamedGenes(
    const std::vector<FastqRead>& reads, const std::string& prefix_adapter,
    const std::string& suffix_adapter, int target_gene_length);

#endif
