#pragma once

// Helpers shared by the block-parallel evaluation loops. Grid points are split
// into independent blocks; with OpenMP the blocks run on worker threads, and
// in serial builds the same loop runs on the calling thread.

#include <cstddef>
#include <cstdint>
#include <exception>
#include <limits>
#include <stdexcept>

namespace chem::kumo::detail {

/// Number of blocks of at most `block_size` points needed for `count` points.
/// OpenMP's canonical loop needs a signed counter (std::int64_t here); throws
/// std::length_error with `what` if the block count would not fit.
inline std::size_t parallel_block_count(std::size_t count, std::size_t block_size,
                                        const char* what) {
    const std::size_t nblocks = count / block_size + (count % block_size != 0);
    if (nblocks > static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max()))
        throw std::length_error(what);
    return nblocks;
}

/// Carries an exception out of a parallel loop.
///
/// An exception must not escape an OpenMP parallel region (the program would
/// terminate). Each block therefore catches everything and calls capture();
/// the first exception is kept, the others are dropped, and the remaining
/// blocks still run. After the loop, rethrow_if_any() rethrows the kept
/// exception on the calling thread.
///
///   FirstException errors;
///   #pragma omp parallel for
///   for (std::int64_t block = 0; block < n; ++block) {
///       try { ... } catch (...) { errors.capture(); }
///   }
///   errors.rethrow_if_any();
class FirstException {
public:
    /// Record the exception being handled. Call only inside a catch block.
    void capture() noexcept {
#ifdef _OPENMP
#pragma omp critical(kumo_first_exception)
#endif
        {
            if (!error_) error_ = std::current_exception();
        }
    }

    void rethrow_if_any() const {
        if (error_) std::rethrow_exception(error_);
    }

private:
    std::exception_ptr error_;
};

} // namespace chem::kumo::detail
