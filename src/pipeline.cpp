// Author: Bruno Jurakic, Martin Saincevic

#include "include/pipeline.h"

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "include/clustering.h"
#include "include/distance.h"
#include "include/fasta_io.h"
#include "include/fastq_parser.h"
#include "include/memory_usage.h"
#include "include/minimap2_runner.h"
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

  // Optional minimap2-based read assignment against provided references.
  if (!config.reference_path.empty()) {
    auto mapping_counts =
        RunMinimap2ReadCounts(config.minimap2_path, config.reference_path,
                              config.input_path, config.minimap_min_mapq);

    int total_mapped = 0;
    std::cout << "Minimap2 mapping counts by reference:\n";
    for (const auto& [reference_name, count] : mapping_counts) {
      total_mapped += count;
      std::cout << "  " << reference_name << ": " << count << "\n";
    }
    std::cout << "  Total mapped reads: " << total_mapped << "\n";

    if (!config.mapping_output_path.empty()) {
      WriteMappingCountsTsv(config.mapping_output_path, mapping_counts);
      std::cout << "Mapping report written to " << config.mapping_output_path
                << "\n";
    }
    std::cout << "\n";
  }

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

  auto named_genes = TrimAndExtractNamedGenes(filtered, prefix_adapter,
                                              suffix_adapter, gene_length);

  // Extract the gene sequences for clustering.
  std::vector<std::string> genes;
  genes.reserve(named_genes.size());
  for (const auto& ng : named_genes) genes.push_back(ng.sequence);
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

  // Write cluster assignments for each read.
  if (!config.cluster_assignments_path.empty()) {
    std::ofstream assignments_file(config.cluster_assignments_path);
    if (!assignments_file.is_open()) {
      throw std::runtime_error("Cannot create file: " +
                               config.cluster_assignments_path);
    }
    assignments_file << "read_name\tcluster_id\n";

    // Mark each read by which significant cluster it belongs to.
    std::vector<int> read_to_cluster(named_genes.size(), -1);
    for (int i = 0; i < (int)significant.size(); ++i) {
      for (int idx : significant[i].member_indices) {
        read_to_cluster[idx] = i + 1;
      }
    }

    for (int i = 0; i < (int)named_genes.size(); ++i) {
      if (read_to_cluster[i] >= 0) {
        assignments_file << named_genes[i].read_name << "\tallele_"
                         << read_to_cluster[i] << "\n";
      } else {
        assignments_file << named_genes[i].read_name << "\tnoise\n";
      }
    }
    std::cout << "Cluster assignments written to "
              << config.cluster_assignments_path << "\n";
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
    std::vector<std::string> eval_rows;
    // Collect rows for optional evaluation report.
    eval_rows.emplace_back("expected_name\tbest_allele\tdistance");

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
        eval_rows.push_back(exp.name + "\tallele_" +
                            std::to_string(best_allele) + "\t" +
                            std::to_string(best_dist));
      } else {
        std::cout << "  " << exp.name << " -> NO MATCH\n";
        eval_rows.push_back(exp.name + "\tNO_MATCH\tNA");
      }
    }

    if (!config.evaluation_output_path.empty()) {
      // Write evaluation rows in TSV format.
      std::ofstream eval_file(config.evaluation_output_path);
      if (!eval_file.is_open()) {
        throw std::runtime_error("Cannot create file: " +
                                 config.evaluation_output_path);
      }
      for (const auto& row : eval_rows) {
        eval_file << row << "\n";
      }
      std::cout << "Evaluation report written to "
                << config.evaluation_output_path << "\n";
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

  std::cout << "\nPeak memory usage: " << GetPeakMemoryMB() << " MB\n";

  return 0;
}

int RunMultiSample(const Config& config) {
  namespace fs = std::filesystem;

  if (!config.reference_path.empty()) {
    std::cerr << "Warning: minimap2 mapping is currently available only in "
                 "--input mode; skipping in --input-dir mode.\n";
  }

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

  std::cout << "\nPeak memory usage: " << GetPeakMemoryMB() << " MB\n";

  return 0;
}
