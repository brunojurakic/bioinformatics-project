// Author: Bruno Jurakic

#include "include/pipeline.h"

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include "include/clustering.h"
#include "include/distance.h"
#include "include/fasta_io.h"
#include "include/fastq_parser.h"
#include "include/sequence_filter.h"

SampleResult ProcessSample(const std::string& input_path,
                           const Config& config) {
  const std::string prefix_adapter = "GATCCTCTCTCTGCAGCACATTTCCTG";
  const std::string suffix_adapter = "CAGCGGCGAGGTGACGCGAA";
  const int gene_length = 249;

  auto reads = ParseFastq(input_path);
  auto histogram = ComputeLengthHistogram(reads);
  int mode_length = FindMostCommonLength(histogram);
  auto filtered = FilterByLength(reads, mode_length, config.length_tolerance);
  auto genes = TrimAndExtractGenes(filtered, prefix_adapter, suffix_adapter,
                                   gene_length);

  // Extract sample name from file path.
  std::string sample_name = std::filesystem::path(input_path).stem().string();

  SampleResult result;
  result.sample_name = sample_name;

  if (genes.empty()) return result;

  auto clusters = GreedyCentroidClustering(genes, config.cluster_threshold);

  // Keep only clusters with enough members and sort by size descending.
  std::vector<Cluster> significant;
  for (auto& cluster : clusters) {
    if ((int)cluster.member_indices.size() >= config.min_cluster_size) {
      significant.push_back(std::move(cluster));
    }
  }
  std::sort(significant.begin(), significant.end(),
            [](const Cluster& a, const Cluster& b) {
              return a.member_indices.size() > b.member_indices.size();
            });

  auto consensi = BuildClusterConsensi(genes, significant);

  for (int i = 0; i < (int)consensi.size(); ++i) {
    result.alleles.push_back(consensi[i]);
    result.cluster_sizes.push_back(significant[i].member_indices.size());
  }

  return result;
}

int RunSingleSample(const Config& config) {
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

int RunMultiSample(const Config& config) {
  namespace fs = std::filesystem;

  // Collect all deer fastq files sorted by name.
  std::vector<std::string> fastq_files;
  for (const auto& entry : fs::directory_iterator(config.input_dir)) {
    std::string filename = entry.path().filename().string();
    if (filename[0] == 'J' && filename.size() > 6 &&
        filename.substr(filename.size() - 6) == ".fastq") {
      fastq_files.push_back(entry.path().string());
    }
  }
  std::sort(fastq_files.begin(), fastq_files.end());

  std::cout << "Found " << fastq_files.size() << " deer samples in "
            << config.input_dir << "\n\n";

  // Process every sample.
  std::vector<SampleResult> results;
  int skipped = 0;
  for (const auto& file : fastq_files) {
    auto result = ProcessSample(file, config);
    if (result.alleles.empty()) {
      if (config.verbose) {
        std::cout << "  " << result.sample_name << ": skipped (no alleles)\n";
      }
      skipped++;
      continue;
    }
    std::cout << "  " << result.sample_name << ": " << result.alleles.size()
              << " alleles\n";
    results.push_back(std::move(result));
  }

  std::cout << "\nProcessed " << results.size() << " samples (" << skipped
            << " skipped)\n";

  // Merge alleles across samples.
  auto unique_alleles = MergeAlleles(results, config.merge_threshold);
  PrintAlleleSummary(unique_alleles, results);

  // Write all unique alleles to output file.
  if (!config.output_path.empty()) {
    std::vector<FastaRecord> records;
    for (int i = 0; i < (int)unique_alleles.size(); ++i) {
      FastaRecord record;
      record.name = "unique_allele_" + std::to_string(i + 1) + " samples=" +
                    std::to_string(unique_alleles[i].occurrences.size());
      record.sequence = unique_alleles[i].sequence;
      records.push_back(record);
    }
    WriteFasta(config.output_path, records);
    std::cout << "\nAll unique alleles written to " << config.output_path
              << "\n";
  }

  return 0;
}
