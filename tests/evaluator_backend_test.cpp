#include "chem/kumo.hpp"
#include "test_support.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace chem::kumo::detail {
struct AoMoTestAccess {
    static auto run(const Molecule& m,
                    const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>>& points,
                    const std::vector<Spin>& components, std::size_t block, int order) {
        return m.density_derivatives_mo(points, components, block, order);
    }
};
} // namespace chem::kumo::detail

using namespace chem::kumo;

namespace {
// AO_MO Evaluator vs the direct MO kernel: both finite and within
// 2e-11 * max(1, |a|, |b|).
void close(double a, double b, const std::string& label) {
    kumo_test::check_close(a, b, 2e-11, 2e-11, label);
}

void run(const Molecule& m, const std::string& name) {
    Grid grid;
    grid.origin = Eigen::Vector3d(0.21, -0.16, 0.34);
    grid.steps << 0.23, 0.01, 0.0, -0.02, 0.19, 0.03, 0.01, 0.0, 0.17;
    grid.counts = Eigen::Vector3i(3, 2, 3);
    Eigen::Matrix<double, 3, Eigen::Dynamic> points(3, 18);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 2; ++j)
            for (int k = 0; k < 3; ++k)
                points.col((i * 2 + j) * 3 + k) =
                    grid.origin + grid.steps * Eigen::Vector3d(i, j, k);
    const std::array<DensityDescriptor, 4> ds{density::total, density::alpha, density::beta,
                                              density::spin};
    const std::array<Spin, 4> spins{Spin::Total, Spin::Alpha, Spin::Beta, Spin::Difference};
    for (std::size_t block : {std::size_t{1}, std::size_t{7}, std::size_t{32}}) {
        Evaluator mo(m, {.threads = 1, .block_size = block, .backend = Backend::AO_MO});
        Evaluator ao(m, {.threads = 1, .block_size = block, .backend = Backend::AO_P});
        Evaluator default_evaluator(m, {.threads = 1, .block_size = block});
        for (int order : {0, 1, 2}) {
            Request req;
            for (auto d : ds) {
                req.add(d.value());
                if (order >= 1) req.add(d.gradient()).add(d.tau());
                if (order >= 2) req.add(d.hessian()).add(d.laplacian());
            }
            const std::string where =
                name + " block " + std::to_string(block) + " order " + std::to_string(order);
            const auto result = mo.evaluate(grid, req);
            const auto default_result = default_evaluator.evaluate(grid, req);
            kumo_test::check(default_result.backend_used() == Backend::AO_MO,
                             where +
                                 ": the default backend selects AO_MO for an eligible molecule");
            kumo_test::check(result.backend_used() == Backend::AO_MO,
                             where + ": explicit AO_MO is reported as the backend used");
            const auto automatic =
                Evaluator(m, {.threads = 1, .block_size = block, .backend = Backend::AUTO})
                    .evaluate(grid, req);
            kumo_test::check(automatic.backend_used() == Backend::AO_MO,
                             where + ": AUTO selects AO_MO for an eligible molecule");
            const auto reference =
                detail::AoMoTestAccess::run(m, points, {spins.begin(), spins.end()}, block, order);
            for (std::size_t c = 0; c < 4; ++c) {
                const auto& ref = reference.at(spins[c]);
                const std::string label = where + " component " + std::to_string(c) +
                                          ": AO_MO Evaluator matches the direct MO kernel, ";
                for (std::size_t p = 0; p < 18; ++p) {
                    close(result.get(ds[c].value()).component()[p], ref[p].value,
                          label + "density");
                    if (order >= 1) {
                        for (int a = 0; a < 3; ++a)
                            close(result.get(ds[c].gradient()).component(a)[p], ref[p].gradient[a],
                                  label + "gradient");
                        close(result.get(ds[c].tau()).component()[p], ref[p].tau, label + "tau");
                    }
                    if (order >= 2) {
                        constexpr int rows[6] = {0, 0, 0, 1, 1, 2};
                        constexpr int cols[6] = {0, 1, 2, 1, 2, 2};
                        for (int a = 0; a < 6; ++a)
                            close(result.get(ds[c].hessian()).component(a)[p],
                                  ref[p].hessian(rows[a], cols[a]), label + "Hessian");
                        close(result.get(ds[c].laplacian()).component()[p], ref[p].laplacian,
                              label + "Laplacian");
                    }
                }
            }
            if (order == 2) {
                const auto ao_result = ao.evaluate(grid, req);
                kumo_test::check(ao_result.backend_used() == Backend::AO_P,
                                 where + ": explicit AO_P is reported as the backend used");
                kumo_test::check(ao_result.contains(density::total.hessian()),
                                 where + ": the AO_P result contains the total-density Hessian");
            }
        }
    }
    kumo_test::check_throws<std::invalid_argument>(
        [&] { Evaluator invalid(m, {.threads = 1, .backend = static_cast<Backend>(99)}); },
        name + ": Evaluator rejects an invalid backend value");
}
} // namespace

int main(int argc, char** argv) {
    if (!kumo_test::require_args(argc, argv, 1, kumo_test::unlimited, "FILE.fchk...")) return 2;
    return kumo_test::run([&] {
        for (int i = 1; i < argc; ++i) {
            auto m = Molecule::from_fchk(argv[i]);
            if (kumo_test::check(m.electronic_structure_info().ao_mo_eligible(),
                                 std::string(argv[i]) + ": the fixture is AO_MO eligible"))
                run(m, argv[i]);
        }
    });
}
