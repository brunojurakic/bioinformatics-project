// Author: Bruno Jurakic

#ifndef FASTQ_PARSER_H
#define FASTQ_PARSER_H

#include <string>
#include <vector>

// Represents a single read from a FASTQ file.
struct FastqRead {
  std::string name;
  std::string sequence;
  std::string quality;
};

// Reads all records from a FASTQ file.
// Each record is 4 lines: @name, sequence, +, quality.
std::vector<FastqRead> ParseFastq(const std::string& filepath);

#endif
