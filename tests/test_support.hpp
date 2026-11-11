#pragma once
// Minimal check helpers for chem-kumo tests. Every failed check prints what
// was expected, so a failing test explains itself in the Meson log
// (`meson test --print-errorlogs`) instead of returning a bare exit code.
// Checks do not stop the test: all failures are reported, then
// kumo_test::exit_code() returns 1 if any check failed.

#include <iostream>
#include <sstream>
#include <string>

namespace kumo_test {

inline int failures = 0;

// Compact representation for messages; keeps small errors such as 3.2e-12
// readable, unlike std::to_string.
inline std::string num(double value) {
    std::ostringstream out;
    out.precision(3);
    out << value;
    return out.str();
}

// Record a failure unless `ok` holds; `what` states the expected behaviour.
inline bool check(bool ok, const std::string& what) {
    if (!ok) {
        std::cerr << "FAIL: " << what << '\n';
        ++failures;
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

inline int exit_code() {
    if (failures == 0) return 0;
    std::cerr << failures << " check(s) failed\n";
    return 1;
}

} // namespace kumo_test
