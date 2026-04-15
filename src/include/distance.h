// Author: Martin Šainčević

#ifndef DISTANCE_H
#define DISTANCE_H

#include <string>

// Returns Hamming distance between two sequences of equal length.
// Throws std::invalid_argument if lengths differ.
int HammingDistance(const std::string& first, const std::string& second);

#endif