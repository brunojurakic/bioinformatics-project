// Author: Martin Šainčević

#ifndef CLUSTERING_H
#define CLUSTERING_H

#include <string>
#include <vector>

// Greedy cluster with one fixed centroid per cluster.
struct Cluster {
  std::string centroid;
  std::vector<int> member_indices;
};

// Assign each sequence to the nearest centroid if distance <= max_distance.
// Otherwise create a new cluster with this sequence as centroid.
std::vector<Cluster> GreedyCentroidClustering(
    const std::vector<std::string>& sequences, int max_distance);

#endif
