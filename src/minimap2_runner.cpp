// Author: Martin Saincevic

#include "include/minimap2_runner.h"

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace {

// PAF format is tab-separated.
std::vector<std::string> SplitTab(const std::string& line) {
  std::vector<std::string> fields;
  std::stringstream ss(line);
  std::string field;
  while (std::getline(ss, field, '\t')) {
    fields.push_back(field);
  }
  return fields;
}

std::string Quote(const std::string& path) { return "\"" + path + "\""; }

}  // namespace

std::map<std::string, int> RunMinimap2ReadCounts(
    const std::string& minimap2_path, const std::string& reference_fasta,
    const std::string& reads_fastq, int min_mapq) {
  namespace fs = std::filesystem;

  const auto now_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(
                          std::chrono::steady_clock::now().time_since_epoch())
                          .count();
  const std::string paf_path =
      (fs::temp_directory_path() /
       ("minimap2_mapping_" + std::to_string(now_ns) + ".paf"))
          .string();

  // Use PAF output to keep parsing simple and deterministic.
  const std::string command = Quote(minimap2_path) +
                              " -x sr --secondary=no -c " +
                              Quote(reference_fasta) + " " +
                              Quote(reads_fastq) + " > " + Quote(paf_path);
  const int exit_code = std::system(command.c_str());
  if (exit_code != 0) {
    throw std::runtime_error("minimap2 failed with exit code " +
                             std::to_string(exit_code));
  }

  std::ifstream paf_file(paf_path);
  if (!paf_file.is_open()) {
    throw std::runtime_error("Cannot open minimap2 output: " + paf_path);
  }

  struct BestHit {
    std::string target;
    int mapq = -1;
    int aligned_block = -1;
  };

  std::unordered_map<std::string, BestHit> best_by_read;
  std::string line;
  while (std::getline(paf_file, line)) {
    if (line.empty()) continue;
    auto fields = SplitTab(line);
    if (fields.size() < 12) continue;

    const std::string& query_name = fields[0];
    const std::string& target_name = fields[5];
    const int aligned_block = std::stoi(fields[10]);
    const int mapq = std::stoi(fields[11]);
    if (mapq < min_mapq) continue;

    // Keep one best alignment per read: MAPQ first, aligned block as tie-break.
    auto it = best_by_read.find(query_name);
    if (it == best_by_read.end() || mapq > it->second.mapq ||
        (mapq == it->second.mapq && aligned_block > it->second.aligned_block)) {
      best_by_read[query_name] = {target_name, mapq, aligned_block};
    }
  }

  // Convert best-hit map to final counts per reference sequence.
  std::map<std::string, int> counts;
  for (const auto& [read_name, hit] : best_by_read) {
    (void)read_name;
    counts[hit.target]++;
  }

  std::error_code ec;
  fs::remove(paf_path, ec);
  return counts;
}

void WriteMappingCountsTsv(const std::string& filepath,
                           const std::map<std::string, int>& counts) {
  std::ofstream output(filepath);
  if (!output.is_open()) {
    throw std::runtime_error("Cannot create file: " + filepath);
  }

  output << "reference_name\tread_count\n";
  for (const auto& [reference_name, read_count] : counts) {
    output << reference_name << "\t" << read_count << "\n";
  }
}
