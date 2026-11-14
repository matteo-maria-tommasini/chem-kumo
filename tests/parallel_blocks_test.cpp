// Block-parallel helpers: block counts, propagation of an exception thrown
// inside a parallel loop to the calling thread, for_each_block (every block
// visited exactly once, first exception rethrown, no call for zero blocks),
// and the point-block loops (consecutive ranges covering every point once,
// one workspace per worker, a failing workspace reported). In OpenMP builds the loops really run on
// several threads; in serial builds the same code runs on one thread.
#include "parallel_blocks.hpp"
#include "test_support.hpp"

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using chem::kumo::detail::FirstException;
using chem::kumo::detail::for_each_block;
using chem::kumo::detail::for_each_point_block;
using chem::kumo::detail::for_each_point_block_with_workspace;
using chem::kumo::detail::parallel_block_count;
using chem::kumo::detail::PointBlock;
using kumo_test::check;

namespace {

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

// for_each_block: every block index runs exactly once, the exception of a
// failing block is rethrown after the loop (the first one captured, so one of
// the failing blocks), and no block runs for nblocks 0.
void check_for_each_block() {
    const std::string schedule = "dynamic schedule";
    for (std::size_t nblocks :
         {std::size_t{1}, std::size_t{2}, std::size_t{7}, std::size_t{64}, std::size_t{1000}}) {
        std::vector<std::atomic<int>> visits(nblocks);
        for_each_block(nblocks, [&](std::size_t block) {
            if (block < nblocks) visits[block].fetch_add(1);
        });
        std::size_t wrong = 0;
        for (const auto& v : visits)
            if (v.load() != 1) ++wrong;
        check(wrong == 0, schedule + ": for_each_block(" + std::to_string(nblocks) +
                              ") visits every block exactly once (" + std::to_string(wrong) +
                              " blocks visited a different number of times)");
    }

    bool called = false;
    for_each_block(0, [&](std::size_t) { called = true; });
    check(!called, schedule + ": for_each_block(0) calls nothing");

    constexpr std::size_t nblocks = 64;
    std::vector<std::atomic<int>> visits(nblocks);
    std::int64_t rethrown = -1;
    try {
        for_each_block(nblocks, [&](std::size_t block) {
            visits[block].fetch_add(1);
            if (block == 7 || block == 40) throw BlockError(static_cast<std::int64_t>(block));
        });
    } catch (const BlockError& e) {
        rethrown = e.block;
    }
    check(rethrown == 7 || rethrown == 40,
          schedule +
              ": for_each_block rethrows the exception of a failing block (7 or 40) "
              "after the loop (got " +
              std::to_string(rethrown) + ")");
    std::size_t wrong = 0;
    for (const auto& v : visits)
        if (v.load() != 1) ++wrong;
    check(wrong == 0, schedule +
                          ": for_each_block still runs every other block once when "
                          "blocks fail (" +
                          std::to_string(wrong) + " blocks visited a different number of times)");
}

// Every point of [0, count) belongs to exactly one block, blocks start at
// multiples of block_size and only the last one may be shorter.
void check_point_blocks() {
    for (std::size_t count : {std::size_t{0}, std::size_t{1}, std::size_t{127}, std::size_t{128},
                              std::size_t{129}, std::size_t{1000}}) {
        std::vector<std::atomic<int>> visits(count);
        std::atomic<int> misplaced{0};
        for_each_point_block(count, 128, "x", [&](PointBlock range) {
            if (range.start % 128 != 0 || range.size == 0 || range.size > 128 ||
                (range.size < 128 && range.start + range.size != count))
                misplaced.fetch_add(1);
            for (std::size_t i = range.start; i < range.start + range.size && i < count; ++i)
                visits[i].fetch_add(1);
        });
        std::size_t wrong = 0;
        for (const auto& v : visits)
            if (v.load() != 1) ++wrong;
        check(wrong == 0 && misplaced.load() == 0,
              "for_each_point_block(" + std::to_string(count) +
                  ", 128) covers every point once with aligned blocks (" + std::to_string(wrong) +
                  " points visited a different number of times, " +
                  std::to_string(misplaced.load()) + " misplaced blocks)");
    }
}

// Workspaces: each worker builds one and reuses it for all its blocks, so
// there are at most as many workspaces as blocks and every block sees a built
// workspace; a throwing make_workspace() is rethrown after the loop.
void check_point_blocks_with_workspace() {
    constexpr std::size_t count = 1000;
    std::vector<std::atomic<int>> visits(count);
    std::atomic<int> made{0}, unbuilt{0};

    struct Workspace {
        std::vector<double> buffer;
    };

    for_each_point_block_with_workspace(
        count, 16, "x",
        [&] {
            made.fetch_add(1);
            return Workspace{std::vector<double>(16, 1.0)};
        },
        [&](PointBlock range, Workspace& w) {
            if (w.buffer.size() != 16) unbuilt.fetch_add(1);
            for (std::size_t i = range.start; i < range.start + range.size; ++i)
                visits[i].fetch_add(1);
        });
    std::size_t wrong = 0;
    for (const auto& v : visits)
        if (v.load() != 1) ++wrong;
    check(wrong == 0, "for_each_point_block_with_workspace covers every point once (" +
                          std::to_string(wrong) + " points visited a different number of times)");
    check(made.load() >= 1 && made.load() <= 63,
          "one workspace per worker, never one per block (made " + std::to_string(made.load()) +
              " for 63 blocks)");
    check(unbuilt.load() == 0, "every block receives a built workspace");

    bool rethrown = false;
    std::atomic<int> ran{0};
    try {
        for_each_point_block_with_workspace(
            count, 16, "x", []() -> Workspace { throw std::runtime_error("workspace"); },
            [&](PointBlock, Workspace&) { ran.fetch_add(1); });
    } catch (const std::runtime_error& e) {
        rethrown = std::string(e.what()) == "workspace";
    }
    check(rethrown, "a throwing make_workspace() is rethrown after the loop");
    check(ran.load() == 0,
          "no block runs without a workspace (ran " + std::to_string(ran.load()) + ")");
}
} // namespace

int main() {
    return kumo_test::run([] {
        check(parallel_block_count(0, 128, "x") == 0, "0 points need 0 blocks");
        check(parallel_block_count(1, 128, "x") == 1, "1 point needs 1 block");
        check(parallel_block_count(256, 128, "x") == 2,
              "256 points in blocks of 128 need 2 blocks");
        check(parallel_block_count(257, 128, "x") == 3,
              "257 points in blocks of 128 need 3 blocks");

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

        check_for_each_block();
        check_point_blocks();
        check_point_blocks_with_workspace();

        if (kumo_test::failures == 0) std::cout << "parallel block helpers OK\n";
    });
}
