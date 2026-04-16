// Author: Bruno Jurakic

#ifndef MULTI_SAMPLE_H
#define MULTI_SAMPLE_H

#include <string>
#include <vector>

// Alleles and cluster sizes for a single sample.
struct SampleResult {
  std::string sample_name;
  std::vector<std::string> alleles;
  std::vector<int> cluster_sizes;
};

// Unique allele from all samples and which samples have it.
struct UniqueAllele {
  std::string sequence;
  std::vector<std::pair<int, int>> occurrences;
};

// Merges same alleles from samples and counts how many samples have each
// allele.
std::vector<UniqueAllele> MergeAlleles(const std::vector<SampleResult>& results,
                                       int merge_threshold);

// Prints summary of unique alleles and which samples have them.
void PrintAlleleSummary(const std::vector<UniqueAllele>& unique_alleles,
                        const std::vector<SampleResult>& results);

#endif
