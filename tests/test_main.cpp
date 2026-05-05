// Author: Martin Saincevic

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "include/clustering.h"
#include "include/distance.h"
#include "include/sequence_filter.h"

// Simple test state holder for pass/fail tracking.
struct TestState {
  int passed = 0;
  int failed = 0;
};

// Records a boolean expectation and prints a message on failure.
void Expect(bool condition, const std::string& message, TestState& state) {
  if (condition) {
    state.passed++;
  } else {
    state.failed++;
    std::cerr << "[FAIL] " << message << "\n";
  }
}

// Expects std::invalid_argument to be thrown by the provided callable.
template <typename Func>
void ExpectInvalidArgument(Func func, const std::string& message,
                           TestState& state) {
  try {
    func();
    state.failed++;
    std::cerr << "[FAIL] " << message << " (no exception)\n";
  } catch (const std::invalid_argument&) {
    state.passed++;
  } catch (...) {
    state.failed++;
    std::cerr << "[FAIL] " << message << " (wrong exception type)\n";
  }
}

// Tests Hamming distance behavior for equal/unequal length inputs.
void TestHammingDistance(TestState& state) {
  Expect(HammingDistance("AAAA", "AAAA") == 0,
         "HammingDistance equal strings", state);
  Expect(HammingDistance("ACGT", "AGGT") == 1,
         "HammingDistance counts mismatches", state);
  ExpectInvalidArgument(
      []() { HammingDistance("AAA", "AA"); },
      "HammingDistance throws on unequal length", state);
}

// Tests consensus generation including tie-break behavior.
void TestMajorityConsensus(TestState& state) {
  std::vector<std::string> sequences = {"ACG", "ATG", "ATG"};
  Cluster cluster;
  cluster.member_indices = {0, 1, 2};
  Expect(BuildMajorityConsensus(sequences, cluster) == "ATG",
         "BuildMajorityConsensus majority vote", state);

  std::vector<std::string> tie_sequences = {"A", "C"};
  Cluster tie_cluster;
  tie_cluster.member_indices = {0, 1};
  Expect(BuildMajorityConsensus(tie_sequences, tie_cluster) == "A",
         "BuildMajorityConsensus tie-break order", state);

  Cluster empty_cluster;
  Expect(BuildMajorityConsensus(sequences, empty_cluster).empty(),
         "BuildMajorityConsensus empty cluster", state);

  std::vector<std::string> bad_lengths = {"AAA", "AA"};
  Cluster bad_cluster;
  bad_cluster.member_indices = {0, 1};
  ExpectInvalidArgument(
      [&]() { BuildMajorityConsensus(bad_lengths, bad_cluster); },
      "BuildMajorityConsensus throws on unequal lengths", state);
}

// Tests adapter trimming for matching and mismatching prefixes.
void TestTrimAdapters(TestState& state) {
  const std::string prefix = "AAA";
  const std::string suffix = "TT";

  Expect(TrimAdapters("AAACCCCTT", prefix, suffix) == "CCCC",
         "TrimAdapters extracts gene", state);
  Expect(TrimAdapters("AATCCCCTT", prefix, suffix) == "CCCC",
         "TrimAdapters allows <=2 mismatches", state);
  Expect(TrimAdapters("CCCCCCTT", prefix, suffix).empty(),
         "TrimAdapters rejects mismatched prefix", state);
  Expect(TrimAdapters("AAAT", prefix, suffix).empty(),
         "TrimAdapters rejects too short read", state);
}

int main() {
  TestState state;
  TestHammingDistance(state);
  TestMajorityConsensus(state);
  TestTrimAdapters(state);

  std::cout << "Tests passed: " << state.passed
            << ", failed: " << state.failed << "\n";
  return state.failed == 0 ? 0 : 1;
}

