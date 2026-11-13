#pragma once
// Minimal check helpers for chem-kumo tests. Every failed check prints what
// was expected, so a failing test explains itself in the Meson log
// (`meson test --print-errorlogs`) instead of returning a bare exit code.
// Checks do not stop the test: all failures are reported, then
// kumo_test::exit_code() returns 1 if any check failed.

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <iostream>
#include <limits>
#include <span>
#include <sstream>
#include <string>
#include <string_view>

namespace kumo_test {

inline int failures = 0;

// Only the first failures are printed: a check inside a loop over grid points
// could otherwise flood the log. All failures are still counted.
inline constexpr int max_printed_failures = 100;

// Compact representation for messages; keeps small errors such as 3.2e-12
// readable, unlike std::to_string. Use a larger precision to show values that
// differ only in late digits.
inline std::string num(double value, int precision = 3) {
    std::ostringstream out;
    out.precision(precision);
    out << value;
    return out.str();
}

// Record a failure unless `ok` holds; `what` states the expected behaviour.
inline bool check(bool ok, const std::string& what) {
    if (!ok) {
        ++failures;
        if (failures <= max_printed_failures)
            std::cerr << "FAIL: " << what << '\n';
        else if (failures == max_printed_failures + 1)
            std::cerr << "(further failures are counted but not printed)\n";
    }
    return ok;
}

// Expect `call` to throw an exception of type Exception.
template <typename Exception, typename Call>
void check_throws(Call&& call, const std::string& what) {
    try {
        call();
    } catch (const Exception&) {
        return;
    } catch (const std::exception& e) {
        check(false, what + " (threw a different exception: " + e.what() + ")");
        return;
    }
    check(false, what + " (no exception thrown)");
}

// Expect `actual` close to `expected`: both finite and
//   |actual - expected| <= max(abs, rel * max(|actual|, |expected|)).
// With rel == abs == t this is the common t * max(1, |actual|, |expected|).
inline bool check_close(double actual, double expected, double rel, double abs,
                        std::string_view label) {
    const double difference = std::abs(actual - expected);
    const double tolerance = std::max(abs, rel * std::max(std::abs(actual), std::abs(expected)));
    if (std::isfinite(actual) && std::isfinite(expected) && difference <= tolerance) return true;
    return check(false, std::string(label) + ": expected " + num(expected, 17) + ", got " +
                            num(actual, 17) + " (|difference| " + num(difference) + ", tolerance " +
                            num(tolerance) + ")");
}

// Largest |actual[i] - expected[i]|, accumulated with std::max from 0 (as the
// tests' own loops did). A size mismatch is a failed check and returns
// infinity, so that any later tolerance check fails as well.
inline double max_abs_diff(std::span<const double> actual, std::span<const double> expected,
                           std::string_view label = "max_abs_diff") {
    if (!check(actual.size() == expected.size(), std::string(label) +
                                                     ": same number of values (expected " +
                                                     std::to_string(expected.size()) + ", got " +
                                                     std::to_string(actual.size()) + ")"))
        return std::numeric_limits<double>::infinity();
    double result = 0.0;
    for (std::size_t i = 0; i < actual.size(); ++i)
        result = std::max(result, std::abs(actual[i] - expected[i]));
    return result;
}

// (actual - expected).cwiseAbs().maxCoeff() for Eigen matrices and vectors;
// expressions are evaluated first. A shape mismatch is a failed check and
// returns infinity; empty matrices give 0.
template <typename A, typename B>
double max_abs_diff(const Eigen::MatrixBase<A>& actual, const Eigen::MatrixBase<B>& expected,
                    std::string_view label = "max_abs_diff") {
    const auto& a = actual.eval();
    const auto& b = expected.eval();
    if (!check(a.rows() == b.rows() && a.cols() == b.cols(),
               std::string(label) + ": same shape (expected " + std::to_string(b.rows()) + "x" +
                   std::to_string(b.cols()) + ", got " + std::to_string(a.rows()) + "x" +
                   std::to_string(a.cols()) + ")"))
        return std::numeric_limits<double>::infinity();
    if (a.size() == 0) return 0.0;
    return (a - b).cwiseAbs().maxCoeff();
}

// Check the number of command-line arguments (not counting the program
// name): on a mismatch print "usage: <prog> <usage>" and return false.
inline constexpr int unlimited = std::numeric_limits<int>::max();

inline bool require_args(int argc, char** argv, int min, int max, std::string_view usage) {
    const int given = argc - 1;
    if (given >= min && given <= max) return true;
    std::cerr << "usage: " << (argc > 0 && argv[0] != nullptr ? argv[0] : "test") << ' ' << usage
              << '\n';
    return false;
}

inline int exit_code() {
    if (failures == 0) return 0;
    std::cerr << failures << " check(s) failed\n";
    return 1;
}

// Run a test body: an unexpected exception is reported as a failure instead
// of terminating the program. Returns exit_code().
template <typename Body> int run(Body&& body) {
    try {
        body();
    } catch (const std::exception& e) {
        check(false, std::string("no unexpected exception (got: ") + e.what() + ")");
    } catch (...) {
        check(false, "no unexpected exception (got a non-standard exception)");
    }
    return exit_code();
}

} // namespace kumo_test
