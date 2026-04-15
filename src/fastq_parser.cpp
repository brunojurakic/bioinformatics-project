// Author: Bruno Jurakic

#include "include/fastq_parser.h"

#include <fstream>
#include <stdexcept>

std::vector<FastqRead> ParseFastq(const std::string& filepath) {
  std::ifstream file(filepath);
  if (!file.is_open()) {
    throw std::runtime_error("Cannot open file: " + filepath);
  }

  std::vector<FastqRead> reads;
  std::string name, sequence, plus_line, quality;

  // Each FASTQ record is 4 lines: @name, sequence, +, quality.
  while (std::getline(file, name)) {
    if (name.empty()) continue;
    if (!name.empty() && name.back() == '\r') name.pop_back();

    if (name[0] != '@') {
      throw std::runtime_error("Expected '@' at start of read name: " + name);
    }

    // Remove the leading '@'.
    name = name.substr(1);

    if (!std::getline(file, sequence)) {
      throw std::runtime_error("Unexpected end of file reading sequence");
    }
    if (!sequence.empty() && sequence.back() == '\r') sequence.pop_back();

    if (!std::getline(file, plus_line)) {
      throw std::runtime_error("Unexpected end of file reading '+' line");
    }

    if (!std::getline(file, quality)) {
      throw std::runtime_error("Unexpected end of file reading quality");
    }
    if (!quality.empty() && quality.back() == '\r') quality.pop_back();

    reads.push_back({name, sequence, quality});
  }

  return reads;
}
