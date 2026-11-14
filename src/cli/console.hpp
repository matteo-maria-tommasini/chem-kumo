#pragma once
// Progress messages of the kumo command-line program (--verbose).
#include "chem/kumo/density.hpp"
#include <chrono>
#include <cstddef>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#ifdef _WIN32
#include <io.h>
#define KUMO_ISATTY(stream) _isatty(_fileno(stream))
#else
#include <unistd.h>
#define KUMO_ISATTY(stream) isatty(fileno(stream))
#endif

namespace kumo_cli {
// Progress messages for --verbose, written to standard error so that standard
// output keeps listing only the files written. On a terminal, evaluations
// also show a progress bar that is redrawn in place; otherwise only the
// start and the end of each step are reported.
class Console {
public:
    using clock = std::chrono::steady_clock;

    explicit Console(bool verbose)
        : verbose_(verbose), terminal_(verbose && KUMO_ISATTY(stderr)), start_(clock::now()) {}

    bool verbose() const noexcept { return verbose_; }

    void say(const std::string& line) const {
        if (verbose_) std::cerr << line << '\n' << std::flush;
    }

    /// Announce a step and start timing it.
    void begin(const std::string& step) {
        say(step);
        step_start_ = clock::now();
        last_percent_ = -1;
        bar_drawn_ = false;
    }

    /// Progress callback for the current step (empty unless on a terminal).
    chem::kumo::ProgressCallback progress() {
        if (!terminal_) return {};
        return [this](std::size_t done, std::size_t total) { draw(done, total); };
    }

    /// End the current step with a summary such as "Evaluated in".
    void end(const std::string& summary) {
        if (bar_drawn_) std::cerr << '\n';
        say("  " + summary + " " + seconds(since(step_start_)));
    }

    void finish() const { say("Finished in " + seconds(since(start_))); }

    static std::string seconds(double s) {
        std::ostringstream out;
        out << std::fixed << std::setprecision(s < 10 ? 2 : 1) << s << " s";
        return out.str();
    }

    /// 1543260 -> "1,543,260".
    static std::string grouped(std::size_t n) {
        auto digits = std::to_string(n);
        for (auto i = static_cast<std::ptrdiff_t>(digits.size()) - 3; i > 0; i -= 3)
            digits.insert(static_cast<std::size_t>(i), ",");
        return digits;
    }

private:
    static double since(clock::time_point t) {
        return std::chrono::duration<double>(clock::now() - t).count();
    }

    // "  [==========>         ]  52%  5.2 s elapsed, about 4.8 s left"
    void draw(std::size_t done, std::size_t total) {
        if (total == 0) return;
        const int percent = static_cast<int>(100 * done / total);
        if (percent == last_percent_) return;
        last_percent_ = percent;
        constexpr int width = 30;
        const int filled = static_cast<int>(width * done / total);
        std::string bar(static_cast<std::size_t>(filled), '=');
        if (filled < width) bar += '>';
        bar.resize(width, ' ');
        const double elapsed = since(step_start_);
        std::ostringstream line;
        line << "\r  [" << bar << "] " << std::setw(3) << percent << "%  " << seconds(elapsed)
             << " elapsed";
        if (done > 0 && done < total)
            line << ", about "
                 << seconds(elapsed * static_cast<double>(total - done) / static_cast<double>(done))
                 << " left";
        line << "    ";
        std::cerr << line.str() << std::flush;
        bar_drawn_ = true;
    }

    bool verbose_;
    bool terminal_;
    clock::time_point start_;
    clock::time_point step_start_ = clock::now();
    int last_percent_ = -1;
    bool bar_drawn_ = false;
};

// "a", "a and b", "a, b and c".
inline std::string join(const std::vector<std::string>& items) {
    std::string out;
    for (std::size_t i = 0; i < items.size(); ++i) {
        if (i > 0) out += i + 1 == items.size() ? " and " : ", ";
        out += items[i];
    }
    return out;
}
} // namespace kumo_cli
