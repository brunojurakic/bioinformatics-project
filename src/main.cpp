// Author: Bruno Jurakic

#include <iostream>
#include <string>

#include "fastq_parser.h"

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <input.fastq>\n";
    return 1;
  }

  std::string input_path = argv[1];
  auto reads = ParseFastq(input_path);
  std::cout << "Loaded " << reads.size() << " reads from " << input_path
            << "\n";

  return 0;
}
