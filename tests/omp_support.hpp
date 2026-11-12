#pragma once
// Thread-count control for tests that compare OpenMP runs against a
// single-thread reference. In serial builds set_threads() does nothing and
// the same comparisons run on the calling thread.

#ifdef _OPENMP
#include <omp.h>
#endif

namespace kumo_test {

// Use exactly `threads` OpenMP threads for the following parallel regions
// (dynamic adjustment of the thread count is disabled).
inline void set_threads(int threads) {
#ifdef _OPENMP
    omp_set_dynamic(0);
    omp_set_num_threads(threads);
#else
    (void)threads;
#endif
}

} // namespace kumo_test
