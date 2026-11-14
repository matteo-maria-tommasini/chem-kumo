#pragma once

// Helpers shared by the block-parallel evaluation loops. Points (of a grid or
// a list) are split into independent blocks; with OpenMP the blocks run on worker threads, and
// in serial builds the same loop runs on the calling thread.

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <functional>
#include <limits>
#include <optional>
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

    /// Record that `points` more points are done and report the new total.
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

/// Call body(block) for block = 0, ..., nblocks - 1, on OpenMP threads
/// (dynamic schedule) when there is more than one block. The first exception
/// thrown by any block is rethrown after all blocks have run; blocks must
/// write disjoint outputs.
template <typename Body> void for_each_block(std::size_t nblocks, const Body& body) {
    FirstException errors;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) if (nblocks > 1)
#endif
    for (std::int64_t block = 0; block < static_cast<std::int64_t>(nblocks); ++block) {
        try {
            body(static_cast<std::size_t>(block));
        } catch (...) {
            errors.capture();
        }
    }
    errors.rethrow_if_any();
}

/// Consecutive points [start, start + size) of a point block.
struct PointBlock {
    std::size_t start;
    std::size_t size;
};

/// Split `count` points into consecutive blocks of at most `block_size` and
/// call body(PointBlock) on each, as for_each_block does. `what` is the
/// parallel_block_count message.
template <typename Body>
void for_each_point_block(std::size_t count, std::size_t block_size, const char* what,
                          const Body& body) {
    for_each_block(parallel_block_count(count, block_size, what), [&](std::size_t block) {
        const std::size_t start = block * block_size;
        body(PointBlock{start, std::min(block_size, count - start)});
    });
}

/// As for_each_point_block, with a workspace per OpenMP worker (the calling
/// thread in serial builds): make_workspace() is called once per worker and
/// body(PointBlock, workspace&) reuses it for all the blocks of that worker.
/// If make_workspace() throws, that worker skips its blocks and the
/// exception is rethrown after the loop.
template <typename MakeWorkspace, typename Body>
void for_each_point_block_with_workspace(std::size_t count, std::size_t block_size,
                                         const char* what, const MakeWorkspace& make_workspace,
                                         const Body& body) {
    const std::size_t nblocks = parallel_block_count(count, block_size, what);
    FirstException errors;
#ifdef _OPENMP
#pragma omp parallel if (nblocks > 1)
#endif
    {
        std::optional<decltype(make_workspace())> workspace;
        try {
            workspace.emplace(make_workspace());
        } catch (...) {
            errors.capture();
        }
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
        for (std::int64_t block = 0; block < static_cast<std::int64_t>(nblocks); ++block) {
            if (!workspace) continue;
            try {
                const std::size_t start = static_cast<std::size_t>(block) * block_size;
                body(PointBlock{start, std::min(block_size, count - start)}, *workspace);
            } catch (...) {
                errors.capture();
            }
        }
    }
    errors.rethrow_if_any();
}

} // namespace chem::kumo::detail
