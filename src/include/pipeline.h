// Author: Bruno Jurakic

#ifndef PIPELINE_H
#define PIPELINE_H

#include "config.h"
#include "multi_sample.h"

// Processes a fastq file and returns its discovered alleles.
SampleResult ProcessSample(const std::string& input_path, const Config& config);

// Runs the single sample pipeline with timing and evaluation.
int RunSingleSample(const Config& config);

// Runs multi sample mode: processes all J*.fastq in a directory,
// merges alleles across samples and prints a summary.
int RunMultiSample(const Config& config);

#endif
