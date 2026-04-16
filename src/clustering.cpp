// Author: Martin Šainčević

#include "include/clustering.h"

#include <array>
#include <limits>
#include <stdexcept>

#include "include/distance.h"

std::vector<Cluster> GreedyCentroidClustering(
    const std::vector<std::string>& sequences, int max_distance) {
  std::vector<Cluster> clusters;
  // Upper bound optimization: in worst case each sequence forms its own
  // cluster.
  clusters.reserve(sequences.size());

  for (int i = 0; i < static_cast<int>(sequences.size()); ++i) {
    const std::string& sequence = sequences[i];

    // Track the nearest existing centroid for the current sequence.
    int best_cluster = -1;
    int best_distance = std::numeric_limits<int>::max();

    for (int c = 0; c < static_cast<int>(clusters.size()); ++c) {
      int distance = HammingDistance(sequence, clusters[c].centroid);
      if (distance < best_distance) {
        best_distance = distance;
        best_cluster = c;
      }
    }

    // Join the closest cluster only if it is within the allowed threshold.
    if (best_cluster != -1 && best_distance <= max_distance) {
      clusters[best_cluster].member_indices.push_back(i);
      continue;
    }

    // Otherwise start a new cluster; centroid is fixed to first member.
    Cluster cluster;
    cluster.centroid = sequence;
    cluster.member_indices.push_back(i);
    clusters.push_back(cluster);
  }

  return clusters;
}

std::string BuildMajorityConsensus(const std::vector<std::string>& sequences,
                                   const Cluster& cluster) {
  // Empty cluster has no representative sequence.
  if (cluster.member_indices.empty()) {
    return "";
  }

  const int first_index = cluster.member_indices.front();
  if (first_index < 0 || first_index >= static_cast<int>(sequences.size())) {
    throw std::invalid_argument("Cluster member index is out of range.");
  }

  // Assume equal read length inside one cluster and validate as we iterate.
  const size_t sequence_length = sequences[first_index].size();
  std::string consensus(sequence_length, 'N');

  for (size_t pos = 0; pos < sequence_length; ++pos) {
    // Count bases at this position: A, C, G, T, other.
    std::array<int, 5> counts = {0, 0, 0, 0, 0};

    for (int member_index : cluster.member_indices) {
      if (member_index < 0 ||
          member_index >= static_cast<int>(sequences.size())) {
        throw std::invalid_argument("Cluster member index is out of range.");
      }
      const std::string& seq = sequences[member_index];
      if (seq.size() != sequence_length) {
        throw std::invalid_argument(
            "All sequences inside a cluster must have equal length.");
      }

      switch (seq[pos]) {
        case 'A':
          ++counts[0];
          break;
        case 'C':
          ++counts[1];
          break;
        case 'G':
          ++counts[2];
          break;
        case 'T':
          ++counts[3];
          break;
        default:
          ++counts[4];
          break;
      }
    }

    const char bases[5] = {'A', 'C', 'G', 'T', 'N'};
    // Deterministic tie-break comes from fixed base order above.
    int best_base = 0;
    for (int i = 1; i < 5; ++i) {
      if (counts[i] > counts[best_base]) {
        best_base = i;
      }
    }
    consensus[pos] = bases[best_base];
  }

  return consensus;
}
