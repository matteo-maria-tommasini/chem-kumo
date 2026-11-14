#include <chem/kumo.hpp>
#include "test_support.hpp"
#include <Eigen/LU>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
using namespace chem::kumo;

int main(int argc, char** argv) {
    if (!kumo_test::require_args(argc, argv, 4, 4, "fchk cube alpha|beta index")) return 2;
    try {
        const std::size_t index = std::stoul(argv[4]);
        const SpinChannel channel = std::string(argv[3]) == "alpha" ? SpinChannel::Alpha
                                    : std::string(argv[3]) == "beta"
                                        ? SpinChannel::Beta
                                        : throw std::invalid_argument("invalid channel");
        std::ifstream in(argv[2]);
        if (!in) throw std::runtime_error("cannot open cube");
        std::string line;
        std::getline(in, line);
        std::getline(in, line);
        Grid g;
        int natoms, nfields;
        if (!(in >> natoms >> g.origin[0] >> g.origin[1] >> g.origin[2] >> nfields) ||
            natoms >= 0 || nfields != 1)
            throw std::runtime_error("invalid orbital cube header");
        for (int axis = 0; axis < 3; ++axis) {
            int n;
            if (!(in >> n >> g.steps(0, axis) >> g.steps(1, axis) >> g.steps(2, axis)) || n <= 0)
                throw std::runtime_error("invalid cube grid");
            g.counts[axis] = n;
        }
        for (int a = 0; a < -natoms; ++a) {
            int z;
            double charge, x, y, zp;
            if (!(in >> z >> charge >> x >> y >> zp)) throw std::runtime_error("invalid atom");
        }
        int nmo, mo;
        if (!(in >> nmo >> mo) || nmo != 1 || mo != static_cast<int>(index))
            throw std::runtime_error("orbital index mismatch in cube");
        // Gaussian standard -3 grid stores 1/6 bohr rounded to six decimals.
        constexpr double h = 1.0 / 6.0;
        for (int a = 0; a < 3; ++a)
            if (std::abs(g.steps(a, a) - h) < 0.500001e-6) g.steps(a, a) = h;
        auto mol = Molecule::from_fchk(argv[1]);
        auto field = mol.orbital(g, index, channel);
        auto values = field.values();
        double max_plus = 0, max_minus = 0, l1_plus = 0, l1_minus = 0, ref_l1 = 0, sq_plus = 0,
               sq_minus = 0, sum_ref2 = 0, sum_calc2 = 0;
        for (double v : values) {
            double ref;
            if (!(in >> ref) || !std::isfinite(ref))
                throw std::runtime_error("truncated or invalid cube data");
            double ep = v - ref, em = v + ref;
            max_plus = std::max(max_plus, std::abs(ep));
            max_minus = std::max(max_minus, std::abs(em));
            l1_plus += std::abs(ep);
            l1_minus += std::abs(em);
            ref_l1 += std::abs(ref);
            sq_plus += ep * ep;
            sq_minus += em * em;
            sum_ref2 += ref * ref;
            sum_calc2 += v * v;
        }
        double extra;
        if (in >> extra) throw std::runtime_error("extra cube values");
        if (!(ref_l1 > 0)) throw std::runtime_error("zero reference amplitude");
        bool inverted = l1_minus < l1_plus;
        double max_abs = inverted ? max_minus : max_plus;
        double relative_l1 = (inverted ? l1_minus : l1_plus) / ref_l1;
        double rmse = std::sqrt((inverted ? sq_minus : sq_plus) / values.size());
        double volume = std::abs(g.steps.determinant());
        std::cout << std::setprecision(12) << "reference=" << argv[2] << " orbital=" << index
                  << " spin=" << argv[3] << " points=" << values.size()
                  << " phase=" << (inverted ? "inverted" : "same") << " max_abs_direct=" << max_plus
                  << " max_abs_phase_aligned=" << max_abs << " rmse=" << rmse
                  << " relative_L1=" << relative_l1 << " norm_kumo=" << sum_calc2 * volume
                  << " norm_gaussian=" << sum_ref2 * volume << '\n';
        // Tolerances: the Gaussian reference Cube files (tests/data/
        // mo_reference_cubes; see also tests/REFERENCE_CUBES.md) carry six
        // significant digits, a relative rounding of up to 5e-6. MO
        // amplitudes are largest near the nuclei, hence an absolute
        // max_abs <= 2e-4; the phase-aligned relative L1 error is bounded by
        // the same rounding, hence 2e-5.
        kumo_test::check(!(max_abs > 2e-4),
                         "phase-aligned max_abs <= 2e-4 against the Gaussian reference (got " +
                             kumo_test::num(max_abs) + ")");
        kumo_test::check(!(relative_l1 > 2e-5),
                         "phase-aligned relative_L1 <= 2e-5 against the Gaussian reference (got " +
                             kumo_test::num(relative_l1) + ")");
        return kumo_test::exit_code();
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << '\n';
        return 1;
    }
}
