// Author: Bruno Jurakic

#include "sequence_filter.h"

#include <iostream>

std::map<int, int> ComputeLengthHistogram(const std::vector<FastqRead>& reads) {
  std::map<int, int> histogram;
  for (const auto& read : reads) {
    histogram[read.sequence.size()]++;
  }
  return histogram;
}

int FindMostCommonLength(const std::map<int, int>& histogram) {
  int best_length = 0;
  int best_count = 0;
  for (const auto& [length, count] : histogram) {
    if (count > best_count) {
      best_count = count;
      best_length = length;
    }
  }
  return best_length;
}

std::vector<FastqRead> FilterByLength(const std::vector<FastqRead>& reads,
                                      int target_length, int tolerance) {
  std::vector<FastqRead> filtered;
  for (const auto& read : reads) {
    int len = read.sequence.size();
    if (len >= target_length - tolerance && len <= target_length + tolerance) {
      filtered.push_back(read);
    }
  }
  return filtered;
}

void PrintLengthStats(const std::map<int, int>& histogram) {
  int total = 0;
  for (const auto& [length, count] : histogram) {
    total += count;
  }

  std::cout << "Length distribution (top 10):\n";
  // Collect entries sorted by count descending.
  std::vector<std::pair<int, int>> sorted_entries(histogram.begin(),
                                                  histogram.end());
  std::sort(sorted_entries.begin(), sorted_entries.end(),
            [](const auto& a, const auto& b) { return a.second > b.second; });

  int shown = 0;
  for (const auto& [length, count] : sorted_entries) {
    if (shown >= 10) break;
    std::cout << "  " << length << " bp: " << count << " reads ("
              << (100.0 * count / total) << "%)\n";
    shown++;
  }
  std::cout << "  Total: " << total << " reads\n";
}
