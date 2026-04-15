// Author: Bruno Jurakic

#include "fasta_io.h"

#include <fstream>
#include <stdexcept>

std::vector<FastaRecord> ReadFasta(const std::string& filepath) {
  std::ifstream file(filepath);
  if (!file.is_open()) {
    throw std::runtime_error("Cannot open file: " + filepath);
  }

  std::vector<FastaRecord> records;
  std::string line;
  FastaRecord current;

  while (std::getline(file, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty()) continue;

    if (line[0] == '>') {
      if (!current.name.empty()) {
        records.push_back(current);
      }
      current.name = line.substr(1);
      current.sequence.clear();
    } else {
      current.sequence += line;
    }
  }

  // Save the last record.
  if (!current.name.empty()) {
    records.push_back(current);
  }

  return records;
}

void WriteFasta(const std::string& filepath,
                const std::vector<FastaRecord>& records) {
  std::ofstream file(filepath);
  if (!file.is_open()) {
    throw std::runtime_error("Cannot create file: " + filepath);
  }

  for (const auto& record : records) {
    file << ">" << record.name << "\n" << record.sequence << "\n";
  }
}
