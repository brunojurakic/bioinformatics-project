// Author: Martin Saincevic

#ifndef DISTANCE_H
#define DISTANCE_H

#include <map>
#include <string>
#include <vector>

// Returns Hamming distance between two sequences of equal length.
// Throws std::invalid_argument if lengths differ.
int HammingDistance(const std::string& first, const std::string& second);

// Computes distance from each sequence to one chosen reference sequence.
// Throws std::invalid_argument if any sequence length differs from reference.
std::vector<int> DistancesToReference(const std::vector<std::string>& sequences,
                                      const std::string& reference);

// Builds a histogram where key=distance and value=number of occurrences.
std::map<int, int> BuildDistanceHistogram(const std::vector<int>& distances);

#endif