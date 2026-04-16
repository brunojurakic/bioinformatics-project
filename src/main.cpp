// Author: Bruno Jurakic, Martin Saincevic

#include <exception>
#include <iostream>

#include "include/config.h"
#include "include/pipeline.h"

int main(int argc, char* argv[]) {
  Config config;
  if (!ParseArgs(argc, argv, config)) {
    PrintUsage(argv[0]);
    return 1;
  }

  try {
    if (!config.input_dir.empty()) {
      return RunMultiSample(config);
    }
    return RunSingleSample(config);
  } catch (const std::exception& ex) {
    std::cerr << "Fatal error: " << ex.what() << "\n";
  } catch (...) {
    std::cerr << "Fatal error: unknown exception\n";
  }

  return 1;
}
