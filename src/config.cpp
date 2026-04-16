// Author: Bruno Jurakic, Martin Saincevic

#include "include/config.h"

#include <iostream>
#include <stdexcept>
#include <string>

namespace {

bool ParseNonNegativeInt(const char* value, const std::string& name,
                         int& out_value) {
  try {
    out_value = std::stoi(value);
  } catch (const std::exception&) {
    std::cerr << "Invalid integer for " << name << ": " << value << "\n";
    return false;
  }

  if (out_value < 0) {
    std::cerr << name << " must be >= 0\n";
    return false;
  }
  return true;
}

}  // namespace

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
      if (!ParseNonNegativeInt(argv[++i], "--cluster-threshold",
                               config.cluster_threshold)) {
        return false;
      }
    } else if (arg == "--length-tolerance" && i + 1 < argc) {
      if (!ParseNonNegativeInt(argv[++i], "--length-tolerance",
                               config.length_tolerance)) {
        return false;
      }
    } else if (arg == "--min-cluster-size" && i + 1 < argc) {
      if (!ParseNonNegativeInt(argv[++i], "--min-cluster-size",
                               config.min_cluster_size)) {
        return false;
      }
    } else if (arg == "--merge-threshold" && i + 1 < argc) {
      if (!ParseNonNegativeInt(argv[++i], "--merge-threshold",
                               config.merge_threshold)) {
        return false;
      }
    } else if (arg == "--verbose") {
      config.verbose = true;
    } else if (arg == "--help" || arg == "-h") {
      return false;
    } else {
      std::cerr << "Unknown argument: " << arg << "\n";
      return false;
    }
  }

  if (!config.input_path.empty() && !config.input_dir.empty()) {
    std::cerr << "Use either --input or --input-dir, not both.\n";
    return false;
  }
  if (config.input_path.empty() && config.input_dir.empty()) {
    std::cerr << "Missing required input: use --input or --input-dir.\n";
    return false;
  }
  if (!config.evaluation_output_path.empty() && config.expected_path.empty()) {
    std::cerr << "--eval-output requires --expected.\n";
    return false;
  }
  if (config.min_cluster_size <= 0) {
    std::cerr << "--min-cluster-size must be > 0\n";
    return false;
  }
  return true;
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
  std::cerr << "\nRules:\n"
            << "  Exactly one of --input or --input-dir must be provided.\n"
            << "  --eval-output can be used only with --expected.\n"
            << "  Threshold arguments must be non-negative; --min-cluster-size "
               "must be > 0.\n";
}
