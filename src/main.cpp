// Author: Bruno Jurakic, Martin Saincevic

#include <algorithm>
#include <chrono>
#include <iostream>
#include <string>

#include "include/clustering.h"
#include "include/distance.h"
#include "include/fasta_io.h"
#include "include/fastq_parser.h"
#include "include/sequence_filter.h"

// Config object for cli arguments.
struct Config {
  std::string input_path;
  std::string output_path;
  std::string expected_path;
  int cluster_threshold = 15;
  int length_tolerance = 5;
  int min_cluster_size = 3;
  bool verbose = false;
};

// Parses cli arguments into Config struct, returns false if something is
// missing.
bool ParseArgs(int argc, char* argv[], Config& config) {
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--input" && i + 1 < argc) {
      config.input_path = argv[++i];
    } else if (arg == "--output" && i + 1 < argc) {
      config.output_path = argv[++i];
    } else if (arg == "--expected" && i + 1 < argc) {
      config.expected_path = argv[++i];
    } else if (arg == "--cluster-threshold" && i + 1 < argc) {
      config.cluster_threshold = std::stoi(argv[++i]);
    } else if (arg == "--length-tolerance" && i + 1 < argc) {
      config.length_tolerance = std::stoi(argv[++i]);
    } else if (arg == "--min-cluster-size" && i + 1 < argc) {
      config.min_cluster_size = std::stoi(argv[++i]);
    } else if (arg == "--verbose") {
      config.verbose = true;
    } else if (arg == "--help" || arg == "-h") {
      return false;
    } else {
      std::cerr << "Unknown argument: " << arg << "\n";
      return false;
    }
  }
  return !config.input_path.empty();
}

void PrintUsage(const std::string& program_name) {
  std::cerr
      << "Usage: " << program_name << " --input <file> [options]\n"
      << "\n"
      << "Options:\n"
      << "  --input <file>            Input FASTQ file (required)\n"
      << "  --output <file>           Output FASTA file with discovered "
         "alleles\n"
      << "  --expected <file>         Expected alleles FASTA for evaluation\n"
      << "  --cluster-threshold <n>   Max Hamming distance within cluster "
         "(default: 15)\n"
      << "  --length-tolerance <n>    Length filter tolerance in bp "
         "(default: 5)\n"
      << "  --min-cluster-size <n>    Minimum reads per cluster (default: 3)\n"
      << "  --verbose                 Print detailed progress info\n"
      << "  --help                    Show this help message\n";
}

int main(int argc, char* argv[]) {
  Config config;
  if (!ParseArgs(argc, argv, config)) {
    PrintUsage(argv[0]);
    return 1;
  }

  using Clock = std::chrono::high_resolution_clock;
  auto pipeline_start = Clock::now();

  // Parse FASTQ file.
  auto t0 = Clock::now();
  auto reads = ParseFastq(config.input_path);
  auto t1 = Clock::now();
  std::cout << "Loaded " << reads.size() << " reads from " << config.input_path
            << "\n";

  // Filter by length.
  auto histogram = ComputeLengthHistogram(reads);
  if (config.verbose) {
    PrintLengthStats(histogram);
  }

  int mode_length = FindMostCommonLength(histogram);
  auto filtered = FilterByLength(reads, mode_length, config.length_tolerance);
  auto t2 = Clock::now();
  std::cout << "Filtered to " << filtered.size() << " reads (" << mode_length
            << " +/- " << config.length_tolerance << " bp)\n";

  // Trim adapters.
  const std::string prefix_adapter = "GATCCTCTCTCTGCAGCACATTTCCTG";
  const std::string suffix_adapter = "CAGCGGCGAGGTGACGCGAA";
  int gene_length = 249;

  auto genes = TrimAndExtractGenes(filtered, prefix_adapter, suffix_adapter,
                                   gene_length);
  auto t3 = Clock::now();
  std::cout << "Trimmed to " << genes.size() << " gene sequences ("
            << gene_length << " bp)\n";

  if (genes.empty()) {
    std::cerr << "No valid gene sequences after filtering/trimming.\n";
    return 1;
  }

  // Cluster sequences.
  auto clusters = GreedyCentroidClustering(genes, config.cluster_threshold);
  auto t4 = Clock::now();

  // Keep only the clusters with enough reads.
  std::vector<Cluster> significant;
  for (auto& cluster : clusters) {
    if ((int)cluster.member_indices.size() >= config.min_cluster_size) {
      significant.push_back(std::move(cluster));
    }
  }

  // Sort by size descending.
  std::sort(significant.begin(), significant.end(),
            [](const Cluster& a, const Cluster& b) {
              return a.member_indices.size() > b.member_indices.size();
            });

  std::cout << "Found " << significant.size() << " alleles (min cluster size "
            << config.min_cluster_size << ")\n";

  // Build consensus for each cluster.
  auto consensi = BuildClusterConsensi(genes, significant);
  auto t5 = Clock::now();

  for (int i = 0; i < (int)consensi.size(); ++i) {
    std::cout << "  Allele " << (i + 1) << ": "
              << significant[i].member_indices.size() << " reads\n";
    if (config.verbose) {
      std::cout << "    " << consensi[i] << "\n";
    }
  }

  // Write output FASTA.
  if (!config.output_path.empty()) {
    std::vector<FastaRecord> records;
    for (int i = 0; i < (int)consensi.size(); ++i) {
      FastaRecord record;
      record.name = "allele_" + std::to_string(i + 1) + " reads=" +
                    std::to_string(significant[i].member_indices.size());
      record.sequence = consensi[i];
      records.push_back(record);
    }
    WriteFasta(config.output_path, records);
    std::cout << "Results written to " << config.output_path << "\n";
  }

  // Compare with expected alleles if provided.
  if (!config.expected_path.empty()) {
    auto expected = ReadFasta(config.expected_path);
    std::cout << "\nEvaluation against " << expected.size()
              << " expected alleles:\n";

    for (const auto& exp : expected) {
      int best_dist = -1;
      int best_allele = -1;
      for (int i = 0; i < (int)consensi.size(); ++i) {
        if (consensi[i].size() != exp.sequence.size()) continue;
        int dist = HammingDistance(consensi[i], exp.sequence);
        if (best_dist < 0 || dist < best_dist) {
          best_dist = dist;
          best_allele = i + 1;
        }
      }
      if (best_allele >= 0) {
        std::cout << "  " << exp.name << " -> allele_" << best_allele
                  << " (distance=" << best_dist << ")\n";
      } else {
        std::cout << "  " << exp.name << " -> NO MATCH\n";
      }
    }
  }

  // Print timing results.
  auto pipeline_end = Clock::now();
  auto ms = [](auto start, auto end) {
    return std::chrono::duration<double, std::milli>(end - start).count();
  };
  std::cout << "\nTiming:\n";
  std::cout << "  Parsing:    " << ms(t0, t1) << " ms\n";
  std::cout << "  Filtering:  " << ms(t1, t2) << " ms\n";
  std::cout << "  Trimming:   " << ms(t2, t3) << " ms\n";
  std::cout << "  Clustering: " << ms(t3, t4) << " ms\n";
  std::cout << "  Consensus:  " << ms(t4, t5) << " ms\n";
  std::cout << "  Total:      " << ms(pipeline_start, pipeline_end) << " ms\n";

  return 0;
}
