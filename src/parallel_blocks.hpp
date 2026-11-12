#pragma once

// Helpers shared by the block-parallel evaluation loops. Grid points are split
// into independent blocks; with OpenMP the blocks run on worker threads, and
// in serial builds the same loop runs on the calling thread.

#include <cstddef>
#include <cstdint>
#include <exception>
#include <functional>
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

/// Reports progress of a block loop to an optional user callback
/// (ProgressCallback in <chem/kumo/density.hpp>).
///
/// advance() may be called from any OpenMP worker; the callback calls are
/// serialized, so it is never run concurrently and `done` increases
/// monotonically. An exception thrown by the callback is rethrown by
/// advance() outside the critical section, so that the block loop can
/// capture it like any other error (see FirstException).
class Progress {
public:
    Progress(const std::function<void(std::size_t, std::size_t)>& callback, std::size_t total)
        : callback_(callback), total_(total) {}

    /// Record that `points` more grid points are done.
    void advance(std::size_t points) {
        if (!callback_) return;
        std::exception_ptr error;
#ifdef _OPENMP
#pragma omp critical(kumo_progress)
#endif
        {
            try {
                done_ += points;
                callback_(done_, total_);
            } catch (...) {
                error = std::current_exception();
            }
        }
        if (error) std::rethrow_exception(error);
    }

private:
    const std::function<void(std::size_t, std::size_t)>& callback_;
    std::size_t total_;
    std::size_t done_ = 0;
};

} // namespace chem::kumo::detail
