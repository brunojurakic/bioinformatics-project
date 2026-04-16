// Author: Martin Saincevic

#ifndef MINIMAP2_RUNNER_H
#define MINIMAP2_RUNNER_H

#include <map>
#include <string>

// Runs minimap2 on reads against a reference FASTA and returns read counts by
// reference sequence name. Keeps one best hit per read and filters by MAPQ.
std::map<std::string, int> RunMinimap2ReadCounts(
    const std::string& minimap2_path, const std::string& reference_fasta,
    const std::string& reads_fastq, int min_mapq);

// Writes mapping counts to a TSV file with columns: reference_name, read_count.
void WriteMappingCountsTsv(const std::string& filepath,
                           const std::map<std::string, int>& counts);

#endif
