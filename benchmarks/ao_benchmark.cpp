#include "chem/kumo/basis_evaluator.hpp"
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv) {
    if (argc < 2 || argc > 4) {
        std::cerr << "usage: kumo-ao-benchmark FILE.fchk [points=100000] [repeats=5]\n";
        return 2;
    }
    try {
        const int points = argc > 2 ? std::stoi(argv[2]) : 100000;
        const int repeats = argc > 3 ? std::stoi(argv[3]) : 5;
        if (points < 1 || repeats < 1) return 2;
        chem::kumo::FormattedCheckpoint fchk(argv[1]);
        chem::kumo::BasisEvaluator evaluator{chem::kumo::BasisSet(fchk)};
        Eigen::VectorXd values(evaluator.basis_size());
        double sink = 0;
        auto run = [&](bool reference) {
            auto begin = std::chrono::steady_clock::now();
            for (int rep = 0; rep < repeats; ++rep) {
                for (int i = 0; i < points; ++i) {
                    const double t = static_cast<double>(i);
                    const Eigen::Vector3d r(0.013*t - 7.0, 3.0*std::sin(0.019*t), 3.0*std::cos(0.023*t));
                    if (reference) evaluator.evaluate_into_reference(r, values);
                    else evaluator.evaluate_into(r, values);
                    sink += values[0];
                }
            }
            return std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count()/repeats;
        };
        const double old_time = run(true);
        const double new_time = run(false);
        std::cout << "basis=" << evaluator.basis_size() << " points=" << points << " repeats=" << repeats << '\n'
                  << "reference_ao_s=" << old_time << '\n' << "optimized_ao_s=" << new_time << '\n'
                  << "speedup=" << old_time/new_time << '\n' << "checksum=" << sink << '\n';
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
