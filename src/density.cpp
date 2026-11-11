#include "backend_selection.hpp"
#include "chem/kumo/density.hpp"
#include "chem/kumo/fchk.hpp"
#include "derivative_layout.hpp"
#include "grid_geometry.hpp"
#include "mo_coefficients.hpp"
#include <algorithm>
#include <array>
#include <cstdint>
#include <exception>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <Eigen/LU>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>

namespace chem::kumo {
namespace {
Eigen::MatrixXd unpack(const std::vector<double>& packed, std::size_t n) {
    if (n > std::numeric_limits<std::size_t>::max() / (n + 1) || packed.size() != n * (n + 1) / 2)
        throw std::runtime_error("invalid packed AO density matrix size");
    Eigen::MatrixXd p(n, n);
    std::size_t q = 0;
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j <= i; ++j) {
            p(i, j) = packed[q];
            p(j, i) = packed[q];
            ++q;
        }
    if (!p.allFinite()) throw std::runtime_error("nonfinite AO density matrix");
    return p;
}

constexpr const char* alpha_count_label = "Number of alpha electrons";
constexpr const char* beta_count_label = "Number of beta electrons";

// Atoms from the FCHK geometry, in bohr.
std::vector<Atom> read_atoms(const FormattedCheckpoint& f) {
    const auto& zs = f.integers("Atomic numbers");
    const auto& xyz = f.reals("Current cartesian coordinates");
    if (zs.size() > std::numeric_limits<std::size_t>::max() / 3 || xyz.size() != 3 * zs.size() ||
        zs.empty())
        throw std::runtime_error("invalid FCHK geometry");
    const std::vector<double>* nuclear_charges = nullptr;
    if (f.contains("Nuclear charges")) {
        const auto& charges = f.reals("Nuclear charges");
        if (charges.size() != zs.size())
            throw std::runtime_error("invalid FCHK nuclear charge count");
        nuclear_charges = &charges;
    }
    std::vector<Atom> atoms;
    for (std::size_t i = 0; i < zs.size(); ++i) {
        if (zs[i] <= 0 || zs[i] > 118 || !std::isfinite(xyz[3 * i]) ||
            !std::isfinite(xyz[3 * i + 1]) || !std::isfinite(xyz[3 * i + 2]))
            throw std::runtime_error("invalid FCHK atom or coordinates");
        // FCHK atomic numbers identify the element. Use an explicit nuclear-charge
        // record when present; otherwise retain the historical all-electron Z
        // convention. An ECP-specific effective charge cannot be inferred here.
        const double charge =
            nuclear_charges == nullptr ? static_cast<double>(zs[i]) : (*nuclear_charges)[i];
        if (!std::isfinite(charge) || charge < 0.0)
            throw std::runtime_error("invalid FCHK nuclear charge");
        atoms.push_back({static_cast<int>(zs[i]), charge,
                         Eigen::Vector3d(xyz[3 * i], xyz[3 * i + 1], xyz[3 * i + 2])});
    }
    return atoms;
}

// Spin-difference AO density matrix. A closed-shell FCHK may omit it (zero
// matrix); an open-shell one, recognized from the electron counts, may not.
Eigen::MatrixXd read_spin_density(const FormattedCheckpoint& f, std::size_t n) {
    if (f.contains("Spin SCF Density")) return unpack(f.reals("Spin SCF Density"), n);
    if (f.contains(alpha_count_label) && f.contains(beta_count_label) &&
        f.integer(alpha_count_label) != f.integer(beta_count_label))
        throw std::runtime_error("open-shell FCHK lacks Spin SCF Density");
    return Eigen::MatrixXd::Zero(n, n);
}

// Numbers of alpha and beta electrons, each at most the basis size n.
std::pair<int, int> read_electron_counts(const FormattedCheckpoint& f, std::size_t n) {
    if (!f.contains(alpha_count_label) || !f.contains(beta_count_label))
        throw std::runtime_error("missing electron counts in FCHK");
    const auto na = f.integer(alpha_count_label);
    const auto nb = f.integer(beta_count_label);
    if (na < 0 || nb < 0 || na > static_cast<long long>(n) || nb > static_cast<long long>(n))
        throw std::runtime_error("invalid electron counts in FCHK");
    return {static_cast<int>(na), static_cast<int>(nb)};
}

// Orbital energies are optional; when both energies and MOs are present,
// there must be one energy per MO.
std::vector<double> read_orbital_energies(const FormattedCheckpoint& f, const std::string& label,
                                          const Eigen::MatrixXd& mos, const char* channel) {
    if (!f.contains(label)) return {};
    const auto& energies = f.reals(label);
    if (mos.cols() && !energies.empty() && energies.size() != static_cast<std::size_t>(mos.cols()))
        throw std::runtime_error(std::string(channel) + " orbital energy count mismatch");
    return energies;
}

// Decide once, at load time, whether the occupied MOs reproduce the FCHK
// density matrices (AO-MO eligible). Missing or too few MOs leave the status
// Unavailable; an ineligible representation never prevents AO-P evaluation.
ElectronicStructureInfo assess_orbital_representation(const Eigen::MatrixXd& total,
                                                      const Eigen::MatrixXd& spin,
                                                      const Eigen::MatrixXd& alpha_mos,
                                                      const Eigen::MatrixXd& beta_mos, int nalpha,
                                                      int nbeta, bool unrestricted, std::size_t n) {
    ElectronicStructureInfo info;
    info.occupied_alpha = static_cast<std::size_t>(nalpha);
    info.occupied_beta = static_cast<std::size_t>(nbeta);
    const bool has_alpha =
        alpha_mos.cols() >= nalpha && alpha_mos.rows() == static_cast<Eigen::Index>(n);
    const bool has_beta = !unrestricted || (beta_mos.cols() >= nbeta &&
                                            beta_mos.rows() == static_cast<Eigen::Index>(n));
    if (!has_alpha || !has_beta) return info;
    const auto decision = detail::validate_orbital_backend(
        total, spin, alpha_mos, beta_mos, info.occupied_alpha, info.occupied_beta, unrestricted);
    info.total_absolute_error = decision.total_absolute_error;
    info.total_relative_error = decision.total_relative_error;
    info.spin_absolute_error = decision.spin_absolute_error;
    info.spin_relative_error = decision.spin_relative_error;
    info.orbital_status = decision.backend == detail::ElectronicBackend::AoMo
                              ? OrbitalRepresentationStatus::Valid
                              : OrbitalRepresentationStatus::Inconsistent;
    return info;
}

// AO density matrix of one spin component. Total and spin difference are the
// canonical FCHK matrices and are returned by reference; alpha and beta are
// built in `storage` as (P_total +/- P_spin) / 2, never from MO occupations.
const Eigen::MatrixXd& select_density_matrix(const Eigen::MatrixXd& total,
                                             const Eigen::MatrixXd& spin, Spin component,
                                             Eigen::MatrixXd& storage) {
    switch (component) {
    case Spin::Total: return total;
    case Spin::Difference: return spin;
    case Spin::Alpha: storage = 0.5 * (total + spin); return storage;
    case Spin::Beta: storage = 0.5 * (total - spin); return storage;
    }
    throw std::invalid_argument("invalid density spin component");
}

// Write a Gaussian Cube file. With an orbital index, follow the Gaussian MO
// convention: negative atom count, then a line with NMO = 1 and the MO label.
void write_field_cube(const std::filesystem::path& path, const std::string& comment,
                      const Grid& grid, const std::vector<Atom>& atoms,
                      const std::vector<double>& values, std::size_t orbital_index = 0) {
    const bool orbital = orbital_index != 0;
    std::ofstream out(path);
    if (!out)
        throw std::runtime_error(orbital ? "cannot open orbital cube output"
                                         : "cannot open cube output");
    out << "chem-kumo " CHEM_KUMO_VERSION "\n"
        << comment << "\n"
        << std::fixed << std::setprecision(6);
    const auto natoms = static_cast<long long>(atoms.size());
    out << std::setw(5) << (orbital ? -natoms : natoms);
    for (int k = 0; k < 3; ++k)
        out << std::setw(13) << grid.origin[k];
    out << '\n';
    for (int a = 0; a < 3; ++a) {
        out << std::setw(5) << grid.counts[a];
        for (int k = 0; k < 3; ++k)
            out << std::setw(13) << grid.steps(k, a);
        out << '\n';
    }
    for (const auto& a : atoms) {
        out << std::setw(5) << a.atomic_number << std::setw(13) << a.nuclear_charge;
        for (int k = 0; k < 3; ++k)
            out << std::setw(13) << a.position[k];
        out << '\n';
    }
    if (orbital) out << std::setw(5) << 1 << std::setw(5) << orbital_index << '\n';
    // Six values per line, and a new line at the end of every k column.
    out << std::scientific << std::setprecision(5);
    for (std::size_t i = 0; i < values.size(); ++i) {
        out << std::setw(13) << values[i];
        if ((i + 1) % 6 == 0 || (i + 1) % static_cast<std::size_t>(grid.counts[2]) == 0)
            out << '\n';
    }
    if (!out) throw std::runtime_error(orbital ? "orbital cube write failed" : "cube write failed");
}
} // namespace

Grid Grid::around(const std::vector<Atom>& atoms, double margin, double spacing) {
    if (atoms.empty() || !std::isfinite(margin) || margin < 0 || !std::isfinite(spacing) ||
        spacing <= 0)
        throw std::invalid_argument("invalid grid parameters");
    if (!atoms.front().position.allFinite() || atoms.front().atomic_number <= 0)
        throw std::invalid_argument("invalid atom in grid");
    Eigen::Vector3d lo = atoms.front().position, hi = lo;
    for (const auto& a : atoms) {
        if (!a.position.allFinite() || a.atomic_number <= 0)
            throw std::invalid_argument("invalid atom in grid");
        lo = lo.cwiseMin(a.position);
        hi = hi.cwiseMax(a.position);
    }
    Grid g;
    g.origin = lo.array() - margin;
    for (int k = 0; k < 3; ++k) {
        const double count = std::ceil((hi[k] - lo[k] + 2 * margin) / spacing) + 1;
        if (!std::isfinite(count) || count > std::numeric_limits<int>::max())
            throw std::invalid_argument("grid too large");
        g.counts[k] = std::max(2, static_cast<int>(count));
    }
    g.steps = spacing * Eigen::Matrix3d::Identity();
    return g;
}

Field::Field(Grid g, std::vector<Atom> atoms, std::vector<double> v)
    : grid_(std::move(g)), atoms_(std::move(atoms)), values_(std::move(v)) {
    if (values_.size() != detail::validated_point_count(grid_))
        throw std::invalid_argument("field size mismatch");
}

double Field::operator()(int i, int j, int k) const {
    if (i < 0 || j < 0 || k < 0 || i >= grid_.counts[0] || j >= grid_.counts[1] ||
        k >= grid_.counts[2])
        throw std::out_of_range("field index");
    return values_.at((static_cast<std::size_t>(i) * grid_.counts[1] + j) * grid_.counts[2] + k);
}

void Field::write_cube(const std::filesystem::path& path, const std::string& comment) const {
    write_field_cube(path, comment, grid_, atoms_, values_);
}

void Field::write_orbital_cube(const std::filesystem::path& path, std::size_t orbital_index,
                               const std::string& comment) const {
    if (orbital_index == 0 || orbital_index > 9999)
        throw std::invalid_argument("invalid orbital index for Cube");
    write_field_cube(path, comment, grid_, atoms_, values_, orbital_index);
}

Molecule::Molecule(std::vector<Atom> a, BasisEvaluator b, Eigen::MatrixXd t, Eigen::MatrixXd s)
    : atoms_(std::move(a)), evaluator_(std::move(b)), total_(std::move(t)), spin_(std::move(s)) {}

// Records are read and validated in a fixed order (geometry, basis, density
// matrices, electron counts, MOs, energies), so a file with several defects
// always reports the first one in this order.
Molecule Molecule::from_fchk(const std::filesystem::path& path) {
    FormattedCheckpoint f(path);
    auto atoms = read_atoms(f);
    BasisEvaluator evaluator{BasisSet(f)};
    const auto n = evaluator.basis_size();
    auto total = unpack(f.reals("Total SCF Density"), n);
    auto spin = read_spin_density(f, n);
    Molecule result(std::move(atoms), std::move(evaluator), std::move(total), std::move(spin));
    std::tie(result.nalpha_, result.nbeta_) = read_electron_counts(f, n);

    // MO coefficients are optional here: densities only need the AO matrices.
    if (f.contains("Alpha MO coefficients"))
        result.alpha_mos_ = detail::read_mo_coefficients(f, "Alpha MO coefficients", n);
    if (f.contains("Beta MO coefficients"))
        result.beta_mos_ = detail::read_mo_coefficients(f, "Beta MO coefficients", n);
    result.unrestricted_ = result.beta_mos_.cols() != 0;
    if (result.nalpha_ != result.nbeta_ && !result.unrestricted_)
        throw std::runtime_error("open-shell FCHK lacks Beta MO coefficients");
    result.alpha_energies_ =
        read_orbital_energies(f, "Alpha Orbital Energies", result.alpha_mos_, "alpha");
    result.beta_energies_ =
        read_orbital_energies(f, "Beta Orbital Energies", result.beta_mos_, "beta");

    result.electronic_info_ = assess_orbital_representation(
        result.total_, result.spin_, result.alpha_mos_, result.beta_mos_, result.nalpha_,
        result.nbeta_, result.unrestricted_, n);
    return result;
}

Eigen::MatrixXd Molecule::density_matrix(Spin component) const {
    Eigen::MatrixXd storage;
    return select_density_matrix(total_, spin_, component, storage);
}

DensityDerivatives Molecule::density_derivatives(const Eigen::Vector3d& point,
                                                 Spin component) const {
    if (!point.allFinite()) throw std::invalid_argument("density evaluation point must be finite");
    // P is symmetric. The alpha/beta split preserves the FCHK density
    // convention and permits the same contractions for all four fields.
    Eigen::MatrixXd storage;
    const Eigen::MatrixXd& p = select_density_matrix(total_, spin_, component, storage);
    const auto ao = evaluator_.evaluate_derivatives(point);
    const Eigen::VectorXd y = p * ao.values;
    const Eigen::Matrix<double, Eigen::Dynamic, 3> pg = p * ao.gradient;
    DensityDerivatives result;
    result.value = ao.values.dot(y);
    result.gradient = 2.0 * ao.gradient.transpose() * y;
    for (int a = 0; a < 3; ++a) {
        for (int b = a; b < 3; ++b) {
            const double value = 2.0 * (ao.hessian.col(detail::ao_hessian_column[a][b]).dot(y) +
                                        ao.gradient.col(a).dot(pg.col(b)));
            result.hessian(a, b) = value;
            result.hessian(b, a) = value;
        }
        result.tau += 0.5 * ao.gradient.col(a).dot(pg.col(a));
    }
    result.laplacian = result.hessian.trace();
    return result;
}

std::vector<DensityDerivatives> Molecule::density_derivatives_blocked(
    const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>>& points, Spin component,
    std::size_t block_size) const {
    if (block_size == 0) throw std::invalid_argument("block size must be positive");
    if (!points.allFinite()) throw std::invalid_argument("density points must be finite");
    Eigen::MatrixXd storage;
    const Eigen::MatrixXd& p = select_density_matrix(total_, spin_, component, storage);
    const auto n = static_cast<Eigen::Index>(basis_size());
    const auto count = static_cast<std::size_t>(points.cols());
    std::vector<DensityDerivatives> result(count);
    // Independent blocks write to disjoint output slots. All Eigen panels
    // and AO derivative buffers are private to each iteration.
    const std::size_t nblocks = count / block_size + (count % block_size != 0);
    if (nblocks > static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max()))
        throw std::length_error("too many derivative blocks for parallel evaluation");
    std::exception_ptr worker_error;
    // Each OpenMP worker owns a reusable panel workspace. In serial builds
    // this region executes once, with exactly the same buffer lifetime.
#ifdef _OPENMP
#pragma omp parallel if (nblocks > 1)
#endif
    {
        const Eigen::Index capacity = static_cast<Eigen::Index>(std::min(block_size, count));
        Eigen::MatrixXd input(n, 4 * capacity);
        Eigen::MatrixXd product(n, 4 * capacity);
        Eigen::MatrixXd hxx(n, capacity), hyy(n, capacity), hzz(n, capacity);
        Eigen::MatrixXd hxy(n, capacity), hxz(n, capacity), hyz(n, capacity);
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
        for (std::int64_t block = 0; block < static_cast<std::int64_t>(nblocks); ++block) {
            try {
                const std::size_t start = static_cast<std::size_t>(block) * block_size;
                const std::size_t size = std::min(block_size, count - start);
                const auto b = static_cast<Eigen::Index>(size);
                // Pack AO values and first derivatives into one column-major panel.
                // A single GEMM applies P to all four panels, avoiding four separate
                // matrix products and improving reuse of the (potentially large) P.
                auto phi = input.middleCols(0, b);
                auto dx = input.middleCols(b, b);
                auto dy = input.middleCols(2 * b, b);
                auto dz = input.middleCols(3 * b, b);
                for (Eigen::Index j = 0; j < b; ++j) {
                    const auto ao = evaluator_.evaluate_derivatives(
                        points.col(static_cast<Eigen::Index>(start) + j));
                    phi.col(j) = ao.values;
                    dx.col(j) = ao.gradient.col(0);
                    dy.col(j) = ao.gradient.col(1);
                    dz.col(j) = ao.gradient.col(2);
                    hxx.col(j) = ao.hessian.col(0);
                    hyy.col(j) = ao.hessian.col(1);
                    hzz.col(j) = ao.hessian.col(2);
                    hxy.col(j) = ao.hessian.col(3);
                    hxz.col(j) = ao.hessian.col(4);
                    hyz.col(j) = ao.hessian.col(5);
                }
                product.leftCols(4 * b).noalias() = p * input.leftCols(4 * b);
                const auto y = product.middleCols(0, b);
                const Eigen::MatrixXd* h[6] = {&hxx, &hyy, &hzz, &hxy, &hxz, &hyz};
                for (Eigen::Index j = 0; j < b; ++j) {
                    auto& out = result[start + static_cast<std::size_t>(j)];
                    out.value = phi.col(j).dot(y.col(j));
                    for (int a = 0; a < 3; ++a) {
                        const auto ga = input.col((a + 1) * b + j);
                        out.gradient[a] = 2.0 * ga.dot(y.col(j));
                        out.tau += 0.5 * ga.dot(product.col((a + 1) * b + j));
                        for (int c = a; c < 3; ++c) {
                            const double v =
                                2.0 * (h[detail::ao_hessian_column[a][c]]->col(j).dot(y.col(j)) +
                                       ga.dot(product.col((c + 1) * b + j)));
                            out.hessian(a, c) = out.hessian(c, a) = v;
                        }
                    }
                    out.laplacian = out.hessian.trace();
                }
            } catch (...) {
#ifdef _OPENMP
#pragma omp critical(kumo_derivatives_block_error)
#endif
                {
                    if (!worker_error) worker_error = std::current_exception();
                }
            }
        }
    } // worker-local workspace lifetime
    if (worker_error) std::rethrow_exception(worker_error);
    return result;
}

// Shared AO evaluation for all requested SCF spin components.
// Total and spin are the canonical FCHK matrices; alpha/beta are linear
// combinations of their contractions, never reconstructed from MO occupations.
std::map<Spin, std::vector<DensityDerivatives>> Molecule::density_derivatives_shared(
    const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>>& points,
    const std::vector<Spin>& components, std::size_t block_size, int ao_order) const {
    if (ao_order != 1 && ao_order != 2)
        throw std::invalid_argument("shared AO derivative order must be 1 or 2");
    if (block_size == 0) throw std::invalid_argument("block size must be positive");
    if (!points.allFinite()) throw std::invalid_argument("density points must be finite");
    const std::size_t count = static_cast<std::size_t>(points.cols());
    std::map<Spin, std::vector<DensityDerivatives>> result;
    bool total_needed = false, spin_needed = false;
    for (Spin component : components) {
        switch (component) {
        case Spin::Total: total_needed = true; break;
        case Spin::Difference: spin_needed = true; break;
        case Spin::Alpha:
        case Spin::Beta: total_needed = spin_needed = true; break;
        default: throw std::invalid_argument("invalid density spin component");
        }
        result.try_emplace(component, count);
    }
    if (result.empty() || count == 0) return result;
    const auto n = static_cast<Eigen::Index>(basis_size());
    const std::size_t nblocks = count / block_size + (count % block_size != 0);
    if (nblocks > static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max()))
        throw std::length_error("too many derivative blocks for parallel evaluation");
    std::exception_ptr worker_error;
#ifdef _OPENMP
#pragma omp parallel if (nblocks > 1)
#endif
    {
        const Eigen::Index capacity = static_cast<Eigen::Index>(std::min(block_size, count));
        Eigen::MatrixXd input(n, 4 * capacity);
        Eigen::MatrixXd product_total, product_spin;
        if (total_needed) product_total.resize(n, 4 * capacity);
        if (spin_needed) product_spin.resize(n, 4 * capacity);
        const Eigen::Index hcapacity = ao_order == 2 ? capacity : 0;
        Eigen::MatrixXd hxx(n, hcapacity), hyy(n, hcapacity), hzz(n, hcapacity);
        Eigen::MatrixXd hxy(n, hcapacity), hxz(n, hcapacity), hyz(n, hcapacity);
        const Eigen::MatrixXd* h[6] = {&hxx, &hyy, &hzz, &hxy, &hxz, &hyz};
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
        for (std::int64_t block = 0; block < static_cast<std::int64_t>(nblocks); ++block) {
            try {
                const std::size_t start = static_cast<std::size_t>(block) * block_size;
                const Eigen::Index b =
                    static_cast<Eigen::Index>(std::min(block_size, count - start));
                auto phi = input.middleCols(0, b);
                for (Eigen::Index j = 0; j < b; ++j) {
                    const auto ao = evaluator_.evaluate_derivatives(
                        points.col(static_cast<Eigen::Index>(start) + j), ao_order);
                    phi.col(j) = ao.values;
                    for (int a = 0; a < 3; ++a)
                        input.col((a + 1) * b + j) = ao.gradient.col(a);
                    if (ao_order == 2) {
                        hxx.col(j) = ao.hessian.col(0);
                        hyy.col(j) = ao.hessian.col(1);
                        hzz.col(j) = ao.hessian.col(2);
                        hxy.col(j) = ao.hessian.col(3);
                        hxz.col(j) = ao.hessian.col(4);
                        hyz.col(j) = ao.hessian.col(5);
                    }
                }
                if (total_needed)
                    product_total.leftCols(4 * b).noalias() = total_ * input.leftCols(4 * b);
                if (spin_needed)
                    product_spin.leftCols(4 * b).noalias() = spin_ * input.leftCols(4 * b);
                // Construct the canonical contractions once, then combine their
                // linear primitive properties for alpha/beta as needed.
                for (Eigen::Index j = 0; j < b; ++j) {
                    DensityDerivatives t, s;
                    const auto contract = [&](const Eigen::MatrixXd& product,
                                              DensityDerivatives& out) {
                        const auto y = product.col(j);
                        out.value = phi.col(j).dot(y);
                        for (int a = 0; a < 3; ++a) {
                            const auto ga = input.col((a + 1) * b + j);
                            out.gradient[a] = 2.0 * ga.dot(y);
                            out.tau += 0.5 * ga.dot(product.col((a + 1) * b + j));
                            if (ao_order == 2)
                                for (int c = a; c < 3; ++c) {
                                    const double v =
                                        2.0 * (h[detail::ao_hessian_column[a][c]]->col(j).dot(y) +
                                               ga.dot(product.col((c + 1) * b + j)));
                                    out.hessian(a, c) = out.hessian(c, a) = v;
                                }
                        }
                        out.laplacian = out.hessian.trace();
                    };
                    if (total_needed) contract(product_total, t);
                    if (spin_needed) contract(product_spin, s);
                    const auto combine = [&](double sign) {
                        DensityDerivatives d;
                        d.value = 0.5 * (t.value + sign * s.value);
                        d.gradient = 0.5 * (t.gradient + sign * s.gradient);
                        d.hessian = 0.5 * (t.hessian + sign * s.hessian);
                        d.laplacian = d.hessian.trace();
                        d.tau = 0.5 * (t.tau + sign * s.tau);
                        return d;
                    };
                    const auto slot = start + static_cast<std::size_t>(j);
                    if (auto it = result.find(Spin::Total); it != result.end())
                        it->second[slot] = t;
                    if (auto it = result.find(Spin::Difference); it != result.end())
                        it->second[slot] = s;
                    if (auto it = result.find(Spin::Alpha); it != result.end())
                        it->second[slot] = combine(1.0);
                    if (auto it = result.find(Spin::Beta); it != result.end())
                        it->second[slot] = combine(-1.0);
                }
            } catch (...) {
#ifdef _OPENMP
#pragma omp critical(kumo_shared_derivatives_error)
#endif
                {
                    if (!worker_error) worker_error = std::current_exception();
                }
            }
        }
    }
    if (worker_error) std::rethrow_exception(worker_error);
    return result;
}

// Validated occupied-MO contraction for electron-density derivatives.
// Evaluator selects this backend when the orbital representation is eligible;
// AO-P remains available as the general density-matrix reference.
std::map<Spin, std::vector<DensityDerivatives>> Molecule::density_derivatives_mo(
    const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>>& points,
    const std::vector<Spin>& components, std::size_t block_size, int ao_order) const {
    if (!electronic_info_.ao_mo_eligible())
        throw std::invalid_argument("AO-MO backend is not eligible for this molecule");
    if (ao_order < 0 || ao_order > 2)
        throw std::invalid_argument("AO-MO derivative order must be 0, 1 or 2");
    if (block_size == 0 || block_size > 65536)
        throw std::invalid_argument("invalid AO-MO block size");
    if (!points.allFinite()) throw std::invalid_argument("nonfinite AO-MO points");
    const std::size_t count = static_cast<std::size_t>(points.cols());
    std::map<Spin, std::vector<DensityDerivatives>> result;
    for (Spin c : components) {
        switch (c) {
        case Spin::Total:
        case Spin::Alpha:
        case Spin::Beta:
        case Spin::Difference: break;
        default: throw std::invalid_argument("invalid AO-MO spin component");
        }
        result.try_emplace(c, count);
    }
    if (result.empty() || count == 0) return result;
    const bool need_a = result.contains(Spin::Alpha) || result.contains(Spin::Total) ||
                        result.contains(Spin::Difference);
    const bool need_b = result.contains(Spin::Beta) || result.contains(Spin::Total) ||
                        result.contains(Spin::Difference);
    const auto n = static_cast<Eigen::Index>(basis_size());
    const auto cap = static_cast<Eigen::Index>(std::min(count, block_size));
    const auto nblocks = count / block_size + (count % block_size != 0);
    if (nblocks > static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max()))
        throw std::length_error("too many AO-MO blocks");
    const auto& cb = unrestricted_ ? beta_mos_ : alpha_mos_;
    const Eigen::MatrixXd ca_occ = alpha_mos_.leftCols(nalpha_);
    const Eigen::MatrixXd cb_occ = cb.leftCols(nbeta_);
    std::exception_ptr worker_error;
#ifdef _OPENMP
#pragma omp parallel if (nblocks > 1)
#endif
    {
        Eigen::MatrixXd ao(n, (ao_order == 0 ? 1 : 4) * cap);
        Eigen::MatrixXd hxx(n, ao_order == 2 ? cap : 0), hyy(n, ao_order == 2 ? cap : 0);
        Eigen::MatrixXd hzz(n, ao_order == 2 ? cap : 0), hxy(n, ao_order == 2 ? cap : 0);
        Eigen::MatrixXd hxz(n, ao_order == 2 ? cap : 0), hyz(n, ao_order == 2 ? cap : 0);
        Eigen::MatrixXd* ah[6] = {&hxx, &hyy, &hzz, &hxy, &hxz, &hyz};
        // Each spin channel is projected once per block, using all occupied MOs.
        const auto project = [&](const Eigen::MatrixXd& c, Eigen::Index b,
                                 std::array<Eigen::MatrixXd, 10>& projected) {
            if (c.cols() == 0) return;
            for (int a = 0; a < (ao_order == 0 ? 1 : 4); ++a)
                projected[a].noalias() = c.transpose() * ao.middleCols(a * b, b);
            if (ao_order == 2)
                for (int a = 0; a < 6; ++a)
                    projected[a + 4].noalias() = c.transpose() * ah[a]->leftCols(b);
        };
        const auto contract = [&](const std::array<Eigen::MatrixXd, 10>& v, Eigen::Index j) {
            DensityDerivatives d;
            if (v[0].size() == 0) return d;
            d.value = v[0].col(j).squaredNorm();
            if (ao_order == 0) return d;
            for (int a = 0; a < 3; ++a) {
                d.gradient[a] = 2.0 * v[0].col(j).dot(v[a + 1].col(j));
                d.tau += 0.5 * v[a + 1].col(j).squaredNorm();
                if (ao_order == 2)
                    for (int c = a; c < 3; ++c) {
                        const double value =
                            2.0 * (v[0].col(j).dot(v[4 + detail::ao_hessian_column[a][c]].col(j)) +
                                   v[a + 1].col(j).dot(v[c + 1].col(j)));
                        d.hessian(a, c) = d.hessian(c, a) = value;
                    }
            }
            d.laplacian = d.hessian.trace();
            return d;
        };
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
        for (std::int64_t block = 0; block < static_cast<std::int64_t>(nblocks); ++block) {
            try {
                const std::size_t start = static_cast<std::size_t>(block) * block_size;
                const auto b = static_cast<Eigen::Index>(std::min(block_size, count - start));
                for (Eigen::Index j = 0; j < b; ++j) {
                    const auto point = points.col(static_cast<Eigen::Index>(start) + j);
                    if (ao_order == 0) {
                        ao.col(j) = evaluator_.evaluate(point);
                    } else {
                        const auto v = evaluator_.evaluate_derivatives(point, ao_order);
                        ao.col(j) = v.values;
                        for (int a = 0; a < 3; ++a)
                            ao.col((a + 1) * b + j) = v.gradient.col(a);
                        if (ao_order == 2)
                            for (int a = 0; a < 6; ++a)
                                ah[a]->col(j) = v.hessian.col(a);
                    }
                }
                std::array<Eigen::MatrixXd, 10> va, vb;
                if (need_a) project(ca_occ, b, va);
                if (need_b) project(cb_occ, b, vb);
                for (Eigen::Index j = 0; j < b; ++j) {
                    const auto a = need_a ? contract(va, j) : DensityDerivatives{};
                    const auto beta = need_b ? contract(vb, j) : DensityDerivatives{};
                    const auto combine = [&](double sign) {
                        DensityDerivatives d;
                        d.value = a.value + sign * beta.value;
                        d.gradient = a.gradient + sign * beta.gradient;
                        d.hessian = a.hessian + sign * beta.hessian;
                        d.laplacian = d.hessian.trace();
                        d.tau = a.tau + sign * beta.tau;
                        return d;
                    };
                    const auto slot = start + static_cast<std::size_t>(j);
                    if (auto it = result.find(Spin::Alpha); it != result.end())
                        it->second[slot] = a;
                    if (auto it = result.find(Spin::Beta); it != result.end())
                        it->second[slot] = beta;
                    if (auto it = result.find(Spin::Total); it != result.end())
                        it->second[slot] = combine(1.0);
                    if (auto it = result.find(Spin::Difference); it != result.end())
                        it->second[slot] = combine(-1.0);
                }
            } catch (...) {
#ifdef _OPENMP
#pragma omp critical(kumo_mo_derivatives_error)
#endif
                {
                    if (!worker_error) worker_error = std::current_exception();
                }
            }
        }
    }
    if (worker_error) std::rethrow_exception(worker_error);
    return result;
}

Field Molecule::density_pointwise(const Grid& grid, Spin spin) const {
    const auto count = detail::validated_point_count(grid);
    std::vector<double> values(count);
    Eigen::MatrixXd storage;
    const Eigen::MatrixXd& p = select_density_matrix(total_, spin_, spin, storage);
    for (std::size_t q = 0; q < count; ++q) {
        const Eigen::VectorXd ao = evaluator_.evaluate(detail::grid_point(grid, q));
        values[q] = ao.dot(p * ao);
    }
    return Field(grid, atoms_, std::move(values));
}

Field Molecule::density(const Grid& grid, Spin spin) const {
    auto fields = densities_blocked(grid, {spin});
    return std::move(fields.at(spin));
}

SpinDensities Molecule::spin_densities(const Grid& grid) const {
    auto fields = densities_blocked(grid, {Spin::Total, Spin::Alpha, Spin::Beta, Spin::Difference});
    return {std::move(fields.at(Spin::Total)), std::move(fields.at(Spin::Alpha)),
            std::move(fields.at(Spin::Beta)), std::move(fields.at(Spin::Difference))};
}

std::map<Spin, Field> Molecule::densities(const Grid& grid,
                                          std::initializer_list<Spin> components) const {
    return densities_blocked(grid, components);
}

std::map<Spin, Field> Molecule::densities_blocked(const Grid& grid,
                                                  std::initializer_list<Spin> components,
                                                  std::size_t block_size) const {
    const auto count = detail::validated_point_count(grid);
    if (block_size == 0 || block_size > 65536)
        throw std::invalid_argument("block_size must be in [1,65536]");
    std::map<Spin, std::vector<double>> buffers;
    for (Spin component : components)
        buffers.try_emplace(component, count);
    if (buffers.empty()) return {};
    const bool need_a = buffers.contains(Spin::Alpha) || buffers.contains(Spin::Total) ||
                        buffers.contains(Spin::Difference);
    const bool need_b = buffers.contains(Spin::Beta) || buffers.contains(Spin::Total) ||
                        buffers.contains(Spin::Difference);
    const Eigen::MatrixXd pa = need_a ? Eigen::MatrixXd(0.5 * (total_ + spin_)) : Eigen::MatrixXd();
    const Eigen::MatrixXd pb = need_b ? Eigen::MatrixXd(0.5 * (total_ - spin_)) : Eigen::MatrixXd();
    const Eigen::Index m = static_cast<Eigen::Index>(basis_size());
    // Each block owns its AO panels and writes to disjoint output indices.
    // Keep the same arithmetic and point ordering in serial and OpenMP builds.
    const std::size_t nblocks = count / block_size + (count % block_size != 0);
    // OpenMP's canonical loop uses a signed counter. Avoid narrowing overflow.
    if (nblocks > static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max()))
        throw std::length_error("too many grid blocks for parallel evaluation");
    std::exception_ptr worker_error;
    // Exceptions must not escape an OpenMP region. Preserve the first error
    // and rethrow it on the caller thread after all workers have finished.
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (nblocks > 1)
#endif
    for (std::int64_t block = 0; block < static_cast<std::int64_t>(nblocks); ++block) {
        try {
            const std::size_t start = static_cast<std::size_t>(block) * block_size;
            const auto b = static_cast<Eigen::Index>(std::min(block_size, count - start));
            Eigen::Matrix<double, 3, Eigen::Dynamic> points(3, b);
            for (Eigen::Index t = 0; t < b; ++t)
                points.col(t) = detail::grid_point(grid, start + static_cast<std::size_t>(t));
            Eigen::MatrixXd phi(m, b);
            evaluator_.evaluate_batch(points, phi);
            Eigen::VectorXd a = Eigen::VectorXd::Zero(b), bb = Eigen::VectorXd::Zero(b);
            if (need_a) {
                Eigen::MatrixXd y = pa * phi;
                a = (phi.array() * y.array()).colwise().sum().transpose().matrix();
            }
            if (need_b) {
                Eigen::MatrixXd y = pb * phi;
                bb = (phi.array() * y.array()).colwise().sum().transpose().matrix();
            }
            for (Eigen::Index t = 0; t < b; ++t) {
                const auto q = start + static_cast<std::size_t>(t);
                if (auto it = buffers.find(Spin::Alpha); it != buffers.end()) it->second[q] = a[t];
                if (auto it = buffers.find(Spin::Beta); it != buffers.end()) it->second[q] = bb[t];
                if (auto it = buffers.find(Spin::Total); it != buffers.end())
                    it->second[q] = a[t] + bb[t];
                if (auto it = buffers.find(Spin::Difference); it != buffers.end())
                    it->second[q] = a[t] - bb[t];
            }
        } catch (...) {
#ifdef _OPENMP
#pragma omp critical(kumo_density_block_error)
#endif
            {
                if (!worker_error) worker_error = std::current_exception();
            }
        }
    }
    if (worker_error) std::rethrow_exception(worker_error);
    std::map<Spin, Field> result;
    for (auto& [spin, values] : buffers)
        result.emplace(spin, Field(grid, atoms_, std::move(values)));
    return result;
}
} // namespace chem::kumo
