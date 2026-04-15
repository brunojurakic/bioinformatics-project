// Author: Bruno Jurakic

#include <iostream>
#include <string>

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

  auto histogram = ComputeLengthHistogram(reads);
  PrintLengthStats(histogram);

  int mode_length = FindMostCommonLength(histogram);
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

  return 0;
}
