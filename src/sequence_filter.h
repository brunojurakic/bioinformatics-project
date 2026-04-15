// Author: Bruno Jurakic

#ifndef SEQUENCE_FILTER_H
#define SEQUENCE_FILTER_H

#include <map>
#include <vector>

#include "fastq_parser.h"

// Counts how many reads have each length.
std::map<int, int> ComputeLengthHistogram(const std::vector<FastqRead>& reads);

// Returns the length that appears most frequently.
int FindMostCommonLength(const std::map<int, int>& histogram);

// Keeps only reads whose length is within [target - tolerance, target + tolerance].
std::vector<FastqRead> FilterByLength(const std::vector<FastqRead>& reads,
                                      int target_length, int tolerance);

// Prints a summary of the length distribution to stdout.
void PrintLengthStats(const std::map<int, int>& histogram);

#endif
