#pragma once
// Optional MPFR proof-building primitives. This header requires MPFR.
#include <mpfr.h>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace chem::kumo::certified {

// Every input double denotes its exact real value. The returned interval
// encloses the real result, rounded outwards to IEEE binary64.
struct Interval { double lower, upper; };

class Number {
public:
    explicit Number(mpfr_prec_t precision = 256) { mpfr_init2(v_, precision); }
    ~Number() { mpfr_clear(v_); }
    Number(const Number&) = delete;
    Number& operator=(const Number&) = delete;
    mpfr_ptr get() noexcept { return v_; }
    mpfr_srcptr get() const noexcept { return v_; }
private:
    mpfr_t v_;
};

enum class Operation { Add, Multiply, Divide, SquareRoot, Exponential };

[[nodiscard]] inline Interval evaluate(Operation op, double x, double y = 0.0) {
    if (!std::isfinite(x) || !std::isfinite(y))
        throw std::invalid_argument("nonfinite MPFR input");
    if (op == Operation::Divide && y == 0.0)
        throw std::invalid_argument("division by zero");
    if (op == Operation::SquareRoot && x < 0.0)
        throw std::invalid_argument("negative square root");
    Number a, b, lo, hi;
    mpfr_set_d(a.get(), x, MPFR_RNDN);
    mpfr_set_d(b.get(), y, MPFR_RNDN);
    auto calculate = [&](mpfr_ptr out, mpfr_rnd_t rnd) {
        switch (op) {
            case Operation::Add: return mpfr_add(out,a.get(),b.get(),rnd);
            case Operation::Multiply: return mpfr_mul(out,a.get(),b.get(),rnd);
            case Operation::Divide: return mpfr_div(out,a.get(),b.get(),rnd);
            case Operation::SquareRoot: return mpfr_sqrt(out,a.get(),rnd);
            case Operation::Exponential: return mpfr_exp(out,a.get(),rnd);
        }
        return 0;
    };
    calculate(lo.get(),MPFR_RNDD);
    calculate(hi.get(),MPFR_RNDU);
    const double lower=mpfr_get_d(lo.get(),MPFR_RNDD);
    const double upper=mpfr_get_d(hi.get(),MPFR_RNDU);
    if (std::isnan(lower) || std::isnan(upper) || lower > upper)
        throw std::runtime_error("invalid MPFR interval");
    return {lower,upper};
}

} // namespace chem::kumo::certified
