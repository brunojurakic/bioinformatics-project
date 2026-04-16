// Author: Bruno Jurakic, Martin Saincevic

#include "include/config.h"
#include "include/pipeline.h"

int main(int argc, char* argv[]) {
  Config config;
  if (!ParseArgs(argc, argv, config)) {
    PrintUsage(argv[0]);
    return 1;
  }

  if (!config.input_dir.empty()) {
    return RunMultiSample(config);
  }
  return RunSingleSample(config);
}
