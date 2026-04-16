// Author: Bruno Jurakic

#include "include/config.h"

#include <iostream>
#include <string>

bool ParseArgs(int argc, char* argv[], Config& config) {
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--input" && i + 1 < argc) {
      config.input_path = argv[++i];
    } else if (arg == "--input-dir" && i + 1 < argc) {
      config.input_dir = argv[++i];
    } else if (arg == "--output" && i + 1 < argc) {
      config.output_path = argv[++i];
    } else if (arg == "--expected" && i + 1 < argc) {
      config.expected_path = argv[++i];
    } else if (arg == "--eval-output" && i + 1 < argc) {
      config.evaluation_output_path = argv[++i];
    } else if (arg == "--cluster-threshold" && i + 1 < argc) {
      config.cluster_threshold = std::stoi(argv[++i]);
    } else if (arg == "--length-tolerance" && i + 1 < argc) {
      config.length_tolerance = std::stoi(argv[++i]);
    } else if (arg == "--min-cluster-size" && i + 1 < argc) {
      config.min_cluster_size = std::stoi(argv[++i]);
    } else if (arg == "--merge-threshold" && i + 1 < argc) {
      config.merge_threshold = std::stoi(argv[++i]);
    } else if (arg == "--verbose") {
      config.verbose = true;
    } else if (arg == "--help" || arg == "-h") {
      return false;
    } else {
      std::cerr << "Unknown argument: " << arg << "\n";
      return false;
    }
  }
  return !config.input_path.empty() || !config.input_dir.empty();
}

void PrintUsage(const std::string& program_name) {
  std::cerr
      << "Usage: " << program_name << " --input <file> [options]\n"
      << "\n"
      << "Options:\n"
      << "  --input <file>            Input FASTQ file\n"
      << "  --input-dir <dir>         Process all J*.fastq files in directory\n"
      << "  --output <file>           Output FASTA file with discovered "
         "alleles\n"
      << "  --expected <file>         Expected alleles FASTA for evaluation\n"
      << "  --eval-output <file>      Output TSV file for evaluation results\n"
      << "  --cluster-threshold <n>   Max Hamming distance within cluster "
         "(default: 15)\n"
      << "  --length-tolerance <n>    Length filter tolerance in bp "
         "(default: 5)\n"
      << "  --min-cluster-size <n>    Minimum reads per cluster (default: 3)\n"
      << "  --merge-threshold <n>     Max distance to merge alleles across "
         "samples (default: 3)\n"
      << "  --verbose                 Print detailed progress info\n"
      << "  --help                    Show this help message\n";
}
