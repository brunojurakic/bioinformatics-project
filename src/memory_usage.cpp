// Author: Bruno Jurakic

#include "include/memory_usage.h"

#include <sys/resource.h>

double GetPeakMemoryMB() {
  struct rusage usage;
  if (getrusage(RUSAGE_SELF, &usage) != 0) {
    return 0.0;
  }
  // Convert from kilobytes to megabytes because on macOS, ru_maxrss is in
  // bytes, and on Linux its in kilobytes.
#ifdef __APPLE__
  return usage.ru_maxrss / (1024.0 * 1024.0);
#else
  return usage.ru_maxrss / 1024.0;
#endif
}
