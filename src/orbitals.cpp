#include <Eigen/LU>
#include "chem/kumo/density.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <cstdint>
#include <exception>

namespace chem::kumo {

// The MO amplitude is a linear contraction of FCHK coefficients with AO
// values. Keep the signed result: squaring would destroy nodal information.

namespace {
void check_channel(Spin s) {
    if (s != Spin::Alpha && s != Spin::Beta)
        throw std::invalid_argument("MO spin channel must be alpha or beta");
}
std::size_t grid_size(const Grid& grid) {
    std::size_t size = 1;
    for (int k=0;k<3;++k) {
        if (grid.counts[k] <= 0 || static_cast<std::size_t>(grid.counts[k]) > std::numeric_limits<std::size_t>::max()/size)
            throw std::invalid_argument("invalid orbital grid size");
        size *= static_cast<std::size_t>(grid.counts[k]);
    }
    if (!grid.origin.allFinite() || !grid.steps.allFinite() ||
        std::abs(grid.steps.determinant()) < 1e-16)
        throw std::invalid_argument("invalid orbital grid geometry");
    return size;
}
}
std::size_t Molecule::orbital_count(Spin channel) const {
    check_channel(channel);
    const auto& matrix = (channel == Spin::Beta && unrestricted_) ? beta_mos_ : alpha_mos_;
    return static_cast<std::size_t>(matrix.cols());
}
int Molecule::occupied_count(Spin channel) const {
    check_channel(channel);
    return channel == Spin::Alpha ? nalpha_ : nbeta_;
}
const std::vector<double>& Molecule::orbital_energies(Spin channel) const {
    check_channel(channel);
    return (channel == Spin::Beta && unrestricted_) ? beta_energies_ : alpha_energies_;
}
OrbitalDerivatives Molecule::orbital_derivatives(
    const Eigen::Vector3d& point, std::size_t index, Spin channel) const {
    check_channel(channel);
    if (!point.allFinite())
        throw std::invalid_argument("MO evaluation point must be finite");
    const auto& coeff = (channel == Spin::Beta && unrestricted_) ? beta_mos_ : alpha_mos_;
    if (index == 0 || index > static_cast<std::size_t>(coeff.cols()))
        throw std::out_of_range("MO index outside available coefficient range");
    const auto ao = evaluator_.evaluate_derivatives(point);
    const auto c = coeff.col(static_cast<Eigen::Index>(index - 1));
    OrbitalDerivatives result;
    result.value = c.dot(ao.values);
    result.gradient = ao.gradient.transpose() * c;
    const Eigen::Matrix<double, 6, 1> h = ao.hessian.transpose() * c;
    result.hessian(0,0) = h[0];
    result.hessian(1,1) = h[1];
    result.hessian(2,2) = h[2];
    result.hessian(0,1) = result.hessian(1,0) = h[3];
    result.hessian(0,2) = result.hessian(2,0) = h[4];
    result.hessian(1,2) = result.hessian(2,1) = h[5];
    return result;
}
std::vector<OrbitalDerivatives> Molecule::orbital_derivatives_blocked(
    const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>>& points,
    std::size_t index, Spin channel, std::size_t block_size) const {
    check_channel(channel);
    const auto& coeff = (channel == Spin::Beta && unrestricted_) ? beta_mos_ : alpha_mos_;
    if (index == 0 || index > static_cast<std::size_t>(coeff.cols()))
        throw std::out_of_range("MO index outside available coefficient range");
    if (block_size == 0) throw std::invalid_argument("orbital derivative block size must be positive");
    if (!points.allFinite()) throw std::invalid_argument("MO evaluation points must be finite");
    const std::size_t count = static_cast<std::size_t>(points.cols());
    std::vector<OrbitalDerivatives> result(count);
    const std::size_t nblocks = count / block_size + (count % block_size != 0);
    if (nblocks > static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max()))
        throw std::length_error("too many orbital derivative blocks for parallel evaluation");
    std::exception_ptr worker_error;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(nblocks > 1)
#endif
    for (std::int64_t block=0; block<static_cast<std::int64_t>(nblocks); ++block) {
        try {
            const std::size_t start = static_cast<std::size_t>(block) * block_size;
            const std::size_t size = std::min(block_size, count-start);
            for (std::size_t j=0; j<size; ++j)
                result[start+j] = orbital_derivatives(points.col(static_cast<Eigen::Index>(start+j)), index, channel);
        } catch (...) {
#ifdef _OPENMP
#pragma omp critical(kumo_orbital_derivatives_block_error)
#endif
            {
                if (!worker_error) worker_error = std::current_exception();
            }
        }
    }
    if (worker_error) std::rethrow_exception(worker_error);
    return result;
}
Field Molecule::orbital(const Grid& grid, std::size_t index, Spin channel) const {
    auto fields = orbitals(grid, {index}, channel);
    return std::move(fields.front());
}
std::vector<Field> Molecule::orbitals(const Grid& grid,
    const std::vector<std::size_t>& indices, Spin channel) const {
    check_channel(channel);
    const auto& coeff = (channel == Spin::Beta && unrestricted_) ? beta_mos_ : alpha_mos_;
    if (coeff.cols() == 0) throw std::runtime_error("FCHK has no MO coefficients");
    for (auto i:indices)
        if (i == 0 || i > static_cast<std::size_t>(coeff.cols()))
            throw std::out_of_range("MO index outside available coefficient range");
    if (indices.empty()) return {};
    const auto count = grid_size(grid);
    std::vector<std::vector<double>> buffers(indices.size(), std::vector<double>(count));
    Eigen::MatrixXd selected(coeff.rows(), static_cast<Eigen::Index>(indices.size()));
    for (std::size_t i=0;i<indices.size();++i)
        selected.col(static_cast<Eigen::Index>(i)) = coeff.col(static_cast<Eigen::Index>(indices[i]-1));
    const int ny=grid.counts[1], nz=grid.counts[2];
    constexpr std::size_t block_size=128;
    const std::size_t nblocks = count / block_size + (count % block_size != 0);
    if (nblocks > static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max()))
        throw std::length_error("too many orbital blocks for parallel evaluation");
    std::exception_ptr worker_error;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(nblocks > 1)
#endif
    for (std::int64_t block=0; block<static_cast<std::int64_t>(nblocks); ++block) {
        try {
        const std::size_t start = static_cast<std::size_t>(block) * block_size;
        const auto size=std::min(block_size,count-start);
        const auto b=static_cast<Eigen::Index>(size);
        Eigen::Matrix<double,3,Eigen::Dynamic> points(3,b);
        for (Eigen::Index t=0;t<b;++t) {
            const auto q=start+static_cast<std::size_t>(t);
            const int k=static_cast<int>(q%nz);
            const int j=static_cast<int>((q/nz)%ny);
            const int i=static_cast<int>(q/(static_cast<std::size_t>(ny)*nz));
            points.col(t)=grid.origin+grid.steps*Eigen::Vector3d(i,j,k);
        }
        Eigen::MatrixXd phi(static_cast<Eigen::Index>(basis_size()),b);
        evaluator_.evaluate_batch(points,phi);
        const Eigen::MatrixXd psi=selected.transpose()*phi;
        for (std::size_t mo=0;mo<indices.size();++mo)
            for (Eigen::Index t=0;t<b;++t)
                buffers[mo][start+static_cast<std::size_t>(t)]=psi(static_cast<Eigen::Index>(mo),t);
        } catch (...) {
#ifdef _OPENMP
#pragma omp critical(kumo_orbital_block_error)
#endif
            {
                if (!worker_error) worker_error = std::current_exception();
            }
        }
    }
    if (worker_error) std::rethrow_exception(worker_error);
    std::vector<Field> result;
    result.reserve(indices.size());
    for (auto& buffer:buffers) result.emplace_back(grid,atoms_,std::move(buffer));
    return result;
}
}
