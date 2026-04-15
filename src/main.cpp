// Author: Bruno Jurakic

#include <iostream>
#include <string>

#include "fastq_parser.h"
#include "sequence_filter.h"

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

  return 0;
}
