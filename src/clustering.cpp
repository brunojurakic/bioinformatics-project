// Author: Martin Šainčević

#include "include/clustering.h"

#include <limits>

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
