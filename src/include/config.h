// Author: Bruno Jurakic

#ifndef CONFIG_H
#define CONFIG_H

#include <string>

// Holds all configuration parsed from cli arguments.
struct Config {
  std::string input_path;
  std::string input_dir;
  std::string output_path;
  std::string expected_path;
  int cluster_threshold = 15;
  int length_tolerance = 5;
  int min_cluster_size = 3;
  int merge_threshold = 3;
  bool verbose = false;
};

// Parses cli arguments into Config struct, returns false if something is
// missing.
bool ParseArgs(int argc, char* argv[], Config& config);

// Prints usage instructions to stderr.
void PrintUsage(const std::string& program_name);

#endif
