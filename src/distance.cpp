// Author: Martin Saincevic

#include "include/distance.h"

#include <stdexcept>
#include <string>

int HammingDistance(const std::string& first, const std::string& second) {
  // Hamming distance is defined only for strings of equal length.
  if (first.size() != second.size()) {
    throw std::invalid_argument(
        "HammingDistance requires sequences of equal length.");
  }

  int mismatches = 0;
  for (size_t i = 0; i < first.size(); ++i) {
    if (first[i] != second[i]) {
      ++mismatches;
    }
  }

  return mismatches;
}

std::vector<int> DistancesToReference(const std::vector<std::string>& sequences,
                                      const std::string& reference) {
  // Preserve input order so indices still match original sequences.
  std::vector<int> distances;
  distances.reserve(sequences.size());

  for (const auto& sequence : sequences) {
    distances.push_back(HammingDistance(sequence, reference));
  }
  return distances;
}

std::map<int, int> BuildDistanceHistogram(const std::vector<int>& distances) {
  // Group equal distances and count their frequency.
  std::map<int, int> histogram;
  for (int d : distances) {
    ++histogram[d];
  }
  return histogram;
}
