// Author: Bruno Jurakic

#include "include/multi_sample.h"

#include <iostream>

#include "include/distance.h"

std::vector<UniqueAllele> MergeAlleles(const std::vector<SampleResult>& results,
                                       int merge_threshold) {
  std::vector<UniqueAllele> unique_alleles;

  for (int s = 0; s < (int)results.size(); ++s) {
    const auto& sample = results[s];

    for (int a = 0; a < (int)sample.alleles.size(); ++a) {
      const auto& allele = sample.alleles[a];
      int cluster_size = sample.cluster_sizes[a];

      // Match against existing unique alleles.
      int best_match = -1;
      int best_dist = merge_threshold + 1;

      for (int u = 0; u < (int)unique_alleles.size(); ++u) {
        if (unique_alleles[u].sequence.size() != allele.size()) continue;
        int dist = HammingDistance(unique_alleles[u].sequence, allele);
        if (dist < best_dist) {
          best_dist = dist;
          best_match = u;
        }
      }

      if (best_match >= 0 && best_dist <= merge_threshold) {
        // Add this sample to the existing unique allele.
        unique_alleles[best_match].occurrences.push_back({s, cluster_size});
      } else {
        UniqueAllele ua;
        ua.sequence = allele;
        ua.occurrences.push_back({s, cluster_size});
        unique_alleles.push_back(ua);
      }
    }
  }

  return unique_alleles;
}

void PrintAlleleSummary(const std::vector<UniqueAllele>& unique_alleles,
                        const std::vector<SampleResult>& results) {
  std::cout << "\n=== Allele Summary Across All Samples ===\n";
  std::cout << "Total unique alleles: " << unique_alleles.size() << "\n";
  std::cout << "Total samples: " << results.size() << "\n\n";

  for (int u = 0; u < (int)unique_alleles.size(); ++u) {
    const auto& ua = unique_alleles[u];
    std::cout << "Unique allele " << (u + 1) << ": found in "
              << ua.occurrences.size() << " sample(s)\n";
    std::cout << "  Sequence: " << ua.sequence.substr(0, 30) << "...\n";
    std::cout << "  Samples: ";
    for (const auto& [sample_idx, size] : ua.occurrences) {
      std::cout << results[sample_idx].sample_name << "(" << size << ") ";
    }
    std::cout << "\n";
  }
}
