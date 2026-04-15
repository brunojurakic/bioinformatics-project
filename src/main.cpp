// Author: Bruno Jurakic

#include <iostream>
#include <string>

#include "include/distance.h"
#include "include/fastq_parser.h"
#include "include/sequence_filter.h"

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <input.fastq>\n";
    return 1;
  }

  std::string input_path = argv[1];
  auto reads = ParseFastq(input_path);
  std::cout << "Loaded " << reads.size() << " reads from " << input_path
            << "\n\n";

  auto length_histogram = ComputeLengthHistogram(reads);
  PrintLengthStats(length_histogram);

  int mode_length = FindMostCommonLength(length_histogram);
  int tolerance = 5;
  std::cout << "\nFiltering to " << mode_length << " +/- " << tolerance
            << " bp\n";

  auto filtered = FilterByLength(reads, mode_length, tolerance);
  std::cout << "Reads after filtering: " << filtered.size() << "\n";

  const std::string prefix_adapter = "GATCCTCTCTCTGCAGCACATTTCCTG";
  const std::string suffix_adapter = "CAGCGGCGAGGTGACGCGAA";
  int gene_length = 249;

  auto genes = TrimAndExtractGenes(filtered, prefix_adapter, suffix_adapter,
                                   gene_length);
  std::cout << "Sequences after trimming: " << genes.size() << "\n";

  if (genes.empty()) {
    std::cerr << "No valid gene sequences after filtering/trimming.\n";
    return 1;
  }

  const std::string& reference = genes.front();
  const auto distances = DistancesToReference(genes, reference);
  const auto distance_histogram = BuildDistanceHistogram(distances);

  std::cout << "Distance histogram vs first sequence (reference):\n";
  for (const auto& [distance, count] : distance_histogram) {
    std::cout << "d=" << distance << " -> " << count << '\n';
  }

  return 0;
}
