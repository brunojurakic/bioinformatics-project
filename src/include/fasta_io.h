// Author: Bruno Jurakic

#ifndef FASTA_IO_H
#define FASTA_IO_H

#include <string>
#include <vector>

// Represents a single record from a FASTA file.
struct FastaRecord {
  std::string name;
  std::string sequence;
};

// Reads all lines from a FASTA file.
std::vector<FastaRecord> ReadFasta(const std::string& filepath);

// Writes records to a FASTA file by writing one sequence per line.
void WriteFasta(const std::string& filepath,
                const std::vector<FastaRecord>& records);

#endif
