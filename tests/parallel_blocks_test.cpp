// Block-parallel helpers: block counts, and propagation of an exception
// thrown inside a parallel loop to the calling thread. In OpenMP builds the
// loop really runs on several threads; in serial builds the same code runs
// on one thread.
#include "parallel_blocks.hpp"

#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using chem::kumo::detail::FirstException;
using chem::kumo::detail::parallel_block_count;

namespace {
int failures = 0;

void check(bool ok, const std::string& what) {
    if (!ok) {
        std::cerr << "FAIL: " << what << '\n';
        ++failures;
    }
}

struct BlockError : std::runtime_error {
    explicit BlockError(std::int64_t block)
        : std::runtime_error("block " + std::to_string(block)), block(block) {}

    std::int64_t block;
};

// Run `nblocks` blocks; blocks listed in `failing` throw BlockError. Returns
// which blocks completed, and stores the rethrown exception's block (or -1).
std::vector<char> run_blocks(std::int64_t nblocks, const std::vector<std::int64_t>& failing,
                             std::int64_t& rethrown_block) {
    std::vector<char> completed(static_cast<std::size_t>(nblocks), 0);
    FirstException errors;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::int64_t block = 0; block < nblocks; ++block) {
        try {
            for (auto f : failing)
                if (block == f) throw BlockError(block);
            completed[static_cast<std::size_t>(block)] = 1;
        } catch (...) {
            errors.capture();
        }
    }
    rethrown_block = -1;
    try {
        errors.rethrow_if_any();
    } catch (const BlockError& e) {
        rethrown_block = e.block;
    }
    return completed;
}
} // namespace

int main() {
    check(parallel_block_count(0, 128, "x") == 0, "0 points need 0 blocks");
    check(parallel_block_count(1, 128, "x") == 1, "1 point needs 1 block");
    check(parallel_block_count(256, 128, "x") == 2, "256 points in blocks of 128 need 2 blocks");
    check(parallel_block_count(257, 128, "x") == 3, "257 points in blocks of 128 need 3 blocks");

    std::int64_t rethrown = 0;
    auto completed = run_blocks(64, {}, rethrown);
    check(rethrown == -1, "no exception is rethrown when no block fails");
    for (std::size_t b = 0; b < completed.size(); ++b)
        check(completed[b] == 1, "block " + std::to_string(b) + " completed without failures");

    completed = run_blocks(64, {7, 40}, rethrown);
    check(rethrown == 7 || rethrown == 40,
          "the exception of a failing block (7 or 40) reaches the caller with its type (got " +
              std::to_string(rethrown) + ")");
    for (std::size_t b = 0; b < completed.size(); ++b)
        if (b != 7 && b != 40)
            check(completed[b] == 1,
                  "block " + std::to_string(b) + " still completed when other blocks failed");

    if (failures != 0) {
        std::cerr << failures << " check(s) failed\n";
        return 1;
    }
    std::cout << "parallel block helpers OK\n";
    return 0;
}
