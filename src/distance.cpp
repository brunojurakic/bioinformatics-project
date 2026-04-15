// Author: Martin Šainčević

#include "include/distance.h"

#include <stdexcept>
#include <string>

int HammingDistance(const std::string& first, const std::string& second) {
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
