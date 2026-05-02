// Author: Bruno Jurakic

#include "include/sequence_filter.h"

#include <algorithm>
#include <iostream>

std::map<int, int> ComputeLengthHistogram(const std::vector<FastqRead>& reads) {
  std::map<int, int> histogram;
  for (const auto& read : reads) {
    histogram[read.sequence.size()]++;
  }
  return histogram;
}

int FindMostCommonLength(const std::map<int, int>& histogram) {
  int best_length = 0;
  int best_count = 0;
  for (const auto& [length, count] : histogram) {
    if (count > best_count) {
      best_count = count;
      best_length = length;
    }
  }
  return best_length;
}

std::vector<FastqRead> FilterByLength(const std::vector<FastqRead>& reads,
                                      int target_length, int tolerance) {
  std::vector<FastqRead> filtered;
  for (const auto& read : reads) {
    int len = read.sequence.size();
    if (len >= target_length - tolerance && len <= target_length + tolerance) {
      filtered.push_back(read);
    }
  }
  return filtered;
}

// Counts mismatches between two strings over the given length.
static int CountMismatches(const std::string& a, const std::string& b,
                           int length) {
  int mismatches = 0;
  for (int i = 0; i < length; ++i) {
    if (a[i] != b[i]) mismatches++;
  }
  return mismatches;
}

std::string TrimAdapters(const std::string& sequence,
                         const std::string& prefix_adapter,
                         const std::string& suffix_adapter) {
  int prefix_len = prefix_adapter.size();
  int suffix_len = suffix_adapter.size();

  // Read must be long enough to have both prefix and suffix.
  if ((int)sequence.size() < prefix_len + suffix_len) return "";

  // Check if prefix matches, must 2 or less mismatches.
  int prefix_mismatches = CountMismatches(sequence, prefix_adapter, prefix_len);
  if (prefix_mismatches > 2) return "";

  // Extract the part between prefix and suffix.
  int gene_start = prefix_len;
  int gene_end = sequence.size() - suffix_len;

  if (gene_end <= gene_start) return "";

  return sequence.substr(gene_start, gene_end - gene_start);
}

std::vector<std::string> TrimAndExtractGenes(
    const std::vector<FastqRead>& reads, const std::string& prefix_adapter,
    const std::string& suffix_adapter, int target_gene_length) {
  std::vector<std::string> genes;
  for (const auto& read : reads) {
    std::string gene =
        TrimAdapters(read.sequence, prefix_adapter, suffix_adapter);
    // Only keep sequences that are exactly the expected gene length.
    if ((int)gene.size() == target_gene_length) {
      genes.push_back(gene);
    }
  }
  return genes;
}

// Trims adapters and returns read_name, gene_sequence for reads that
// have the expected gene length after trimming.
std::vector<NamedGene> TrimAndExtractNamedGenes(
    const std::vector<FastqRead>& reads, const std::string& prefix_adapter,
    const std::string& suffix_adapter, int target_gene_length) {
  std::vector<NamedGene> named_genes;
  for (const auto& read : reads) {
    std::string gene =
        TrimAdapters(read.sequence, prefix_adapter, suffix_adapter);
    if ((int)gene.size() == target_gene_length) {
      named_genes.push_back({read.name, gene});
    }
  }
  return named_genes;
}

void PrintLengthStats(const std::map<int, int>& histogram) {
  int total = 0;
  for (const auto& [length, count] : histogram) {
    total += count;
  }

  std::cout << "Length distribution (top 10):\n";
  // Collect entries sorted by count descending.
  std::vector<std::pair<int, int>> sorted_entries(histogram.begin(),
                                                  histogram.end());
  std::sort(sorted_entries.begin(), sorted_entries.end(),
            [](const auto& a, const auto& b) { return a.second > b.second; });

  int shown = 0;
  for (const auto& [length, count] : sorted_entries) {
    if (shown >= 10) break;
    std::cout << "  " << length << " bp: " << count << " reads ("
              << (100.0 * count / total) << "%)\n";
    shown++;
  }
  std::cout << "  Total: " << total << " reads\n";
}
