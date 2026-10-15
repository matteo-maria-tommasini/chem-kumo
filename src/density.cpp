#include "chem/kumo/density.hpp"
#include "chem/kumo/fchk.hpp"
#include <algorithm>
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

namespace chem::kumo {
namespace {
std::size_t checked_size(const Grid& g) {
    std::size_t n = 1;
    for (int k = 0; k < 3; ++k) {
        if (g.counts[k] <= 0 || static_cast<std::size_t>(g.counts[k]) > std::numeric_limits<std::size_t>::max()/n)
            throw std::invalid_argument("invalid or oversized grid");
        n *= static_cast<std::size_t>(g.counts[k]);
    }
    if (!g.origin.allFinite() || !g.steps.allFinite() || std::abs(g.steps.determinant()) < 1e-16)
        throw std::invalid_argument("invalid grid geometry");
    return n;
}
Eigen::MatrixXd unpack(const std::vector<double>& packed, std::size_t n) {
    if (n > std::numeric_limits<std::size_t>::max() / (n + 1) || packed.size() != n*(n+1)/2) throw std::runtime_error("invalid packed AO density matrix size");
    Eigen::MatrixXd p(n,n);
    std::size_t q=0;
    for (std::size_t i=0;i<n;++i)
        for (std::size_t j=0;j<=i;++j) { p(i,j)=packed[q]; p(j,i)=packed[q]; ++q; }
    if (!p.allFinite()) throw std::runtime_error("nonfinite AO density matrix");
    return p;
}
}
Grid Grid::around(const std::vector<Atom>& atoms, double margin, double spacing) {
    if (atoms.empty() || !std::isfinite(margin) || margin < 0 || !std::isfinite(spacing) || spacing <= 0)
        throw std::invalid_argument("invalid grid parameters");
    if (!atoms.front().position.allFinite() || atoms.front().atomic_number <= 0) throw std::invalid_argument("invalid atom in grid");
    Eigen::Vector3d lo=atoms.front().position, hi=lo;
    for (const auto& a: atoms) {
        if (!a.position.allFinite() || a.atomic_number <= 0) throw std::invalid_argument("invalid atom in grid");
        lo=lo.cwiseMin(a.position); hi=hi.cwiseMax(a.position);
    }
    Grid g;
    g.origin=lo.array()-margin;
    for (int k=0;k<3;++k) {
        const double count=std::ceil((hi[k]-lo[k]+2*margin)/spacing)+1;
        if (!std::isfinite(count) || count > std::numeric_limits<int>::max()) throw std::invalid_argument("grid too large");
        g.counts[k]=std::max(2,static_cast<int>(count));
    }
    g.steps=spacing*Eigen::Matrix3d::Identity();
    return g;
}
Field::Field(Grid g, std::vector<Atom> atoms, std::vector<double> v):grid_(std::move(g)),atoms_(std::move(atoms)),values_(std::move(v)) {
    if (values_.size()!=checked_size(grid_)) throw std::invalid_argument("field size mismatch");
}
double Field::operator()(int i,int j,int k) const {
    if (i<0 || j<0 || k<0 || i>=grid_.counts[0] || j>=grid_.counts[1] || k>=grid_.counts[2]) throw std::out_of_range("field index");
    return values_.at((static_cast<std::size_t>(i)*grid_.counts[1]+j)*grid_.counts[2]+k);
}
void Field::write_cube(const std::filesystem::path& path, const std::string& comment) const {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot open cube output");
    out << "chem-kumo " CHEM_KUMO_VERSION "\n" << comment << "\n" << std::fixed << std::setprecision(6);
    out << std::setw(5) << atoms_.size();
    for (int k=0;k<3;++k) out << std::setw(13) << grid_.origin[k];
    out << '\n';
    for (int a=0;a<3;++a) {
        out << std::setw(5) << grid_.counts[a];
        for (int k=0;k<3;++k) out << std::setw(13) << grid_.steps(k,a);
        out << '\n';
    }
    for (const auto& a:atoms_) {
        out << std::setw(5) << a.atomic_number << std::setw(13) << static_cast<double>(a.atomic_number);
        for (int k=0;k<3;++k) out << std::setw(13) << a.position[k];
        out << '\n';
    }
    out << std::scientific << std::setprecision(5);
    for (std::size_t i=0;i<values_.size();++i) {
        out << std::setw(13) << values_[i];
        if ((i+1)%6==0 || (i+1)%static_cast<std::size_t>(grid_.counts[2])==0) out << '\n';
    }
    if (!out) throw std::runtime_error("cube write failed");
}
void Field::write_orbital_cube(const std::filesystem::path& path, std::size_t orbital_index, const std::string& comment) const {
    if (orbital_index == 0 || orbital_index > 9999) throw std::invalid_argument("invalid orbital index for Cube");
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot open orbital cube output");
    out << "chem-kumo " CHEM_KUMO_VERSION "\n" << comment << "\n" << std::fixed << std::setprecision(6);
    out << std::setw(5) << -static_cast<long long>(atoms_.size());
    for (int k=0;k<3;++k) out << std::setw(13) << grid_.origin[k];
    out << '\n';
    for (int a=0;a<3;++a) {
        out << std::setw(5) << grid_.counts[a];
        for (int k=0;k<3;++k) out << std::setw(13) << grid_.steps(k,a);
        out << '\n';
    }
    for (const auto& a:atoms_) {
        out << std::setw(5) << a.atomic_number << std::setw(13) << static_cast<double>(a.atomic_number);
        for (int k=0;k<3;++k) out << std::setw(13) << a.position[k];
        out << '\n';
    }
    // Gaussian MO Cube convention: negative atom count, then NMO and MO labels.
    out << std::setw(5) << 1 << std::setw(5) << orbital_index << '\n';
    out << std::scientific << std::setprecision(5);
    for (std::size_t i=0;i<values_.size();++i) {
        out << std::setw(13) << values_[i];
        if ((i+1)%6==0 || (i+1)%static_cast<std::size_t>(grid_.counts[2])==0) out << '\n';
    }
    if (!out) throw std::runtime_error("orbital cube write failed");
}
Molecule::Molecule(std::vector<Atom> a, BasisEvaluator b, Eigen::MatrixXd t, Eigen::MatrixXd s)
    :atoms_(std::move(a)),evaluator_(std::move(b)),total_(std::move(t)),spin_(std::move(s)) {}
Molecule Molecule::from_fchk(const std::filesystem::path& path) {
    FormattedCheckpoint f(path);
    const auto& zs=f.integers("Atomic numbers");
    const auto& xyz=f.reals("Current cartesian coordinates");
    if (zs.size() > std::numeric_limits<std::size_t>::max()/3 || xyz.size()!=3*zs.size() || zs.empty()) throw std::runtime_error("invalid FCHK geometry");
    std::vector<Atom> atoms;
    for (std::size_t i=0;i<zs.size();++i) {
        if (zs[i] <= 0 || zs[i] > 118 || !std::isfinite(xyz[3*i]) || !std::isfinite(xyz[3*i+1]) || !std::isfinite(xyz[3*i+2]))
            throw std::runtime_error("invalid FCHK atom or coordinates");
        atoms.push_back({static_cast<int>(zs[i]), Eigen::Vector3d(xyz[3*i],xyz[3*i+1],xyz[3*i+2])});
    }
    BasisEvaluator evaluator{BasisSet(f)};
    const auto n=evaluator.basis_size();
    auto total=unpack(f.reals("Total SCF Density"),n);
    Eigen::MatrixXd spin=Eigen::MatrixXd::Zero(n,n);
    if (f.contains("Spin SCF Density")) spin=unpack(f.reals("Spin SCF Density"),n);
    else if (f.contains("Number of alpha electrons") && f.contains("Number of beta electrons") &&
             f.integer("Number of alpha electrons")!=f.integer("Number of beta electrons"))
        throw std::runtime_error("open-shell FCHK lacks Spin SCF Density");
    Molecule result(std::move(atoms),std::move(evaluator),std::move(total),std::move(spin));
    if (!f.contains("Number of alpha electrons") || !f.contains("Number of beta electrons"))
        throw std::runtime_error("missing electron counts in FCHK");
    const auto na = f.integer("Number of alpha electrons");
    const auto nb = f.integer("Number of beta electrons");
    if (na < 0 || nb < 0 || na > static_cast<long long>(n) || nb > static_cast<long long>(n))
        throw std::runtime_error("invalid electron counts in FCHK");
    result.nalpha_ = static_cast<int>(na);
    result.nbeta_ = static_cast<int>(nb);
    const auto load_mos = [&](const std::string& label, Eigen::MatrixXd& dest) {
        if (!f.contains(label)) return;
        const auto& coefficients = f.reals(label);
        if (n == 0 || coefficients.size() % n != 0 || coefficients.size()/n > n)
            throw std::runtime_error("invalid MO coefficient dimensions: " + label);
        const std::size_t count = coefficients.size()/n;
        dest.resize(static_cast<Eigen::Index>(n), static_cast<Eigen::Index>(count));
        // Gaussian FCHK: consecutive coefficients for each MO, AO index fastest.
        for (std::size_t i=0; i<count; ++i)
            for (std::size_t mu=0; mu<n; ++mu)
                dest(static_cast<Eigen::Index>(mu),static_cast<Eigen::Index>(i)) = coefficients[i*n+mu];
        if (!dest.allFinite()) throw std::runtime_error("nonfinite MO coefficients");
    };
    load_mos("Alpha MO coefficients", result.alpha_mos_);
    load_mos("Beta MO coefficients", result.beta_mos_);
    result.unrestricted_ = result.beta_mos_.cols() != 0;
    if (result.nalpha_ != result.nbeta_ && !result.unrestricted_)
        throw std::runtime_error("open-shell FCHK lacks Beta MO coefficients");
    if (f.contains("Alpha Orbital Energies")) result.alpha_energies_ = f.reals("Alpha Orbital Energies");
    if (f.contains("Beta Orbital Energies")) result.beta_energies_ = f.reals("Beta Orbital Energies");
    if (result.alpha_mos_.cols() && !result.alpha_energies_.empty() &&
        result.alpha_energies_.size() != static_cast<std::size_t>(result.alpha_mos_.cols()))
        throw std::runtime_error("alpha orbital energy count mismatch");
    if (result.beta_mos_.cols() && !result.beta_energies_.empty() &&
        result.beta_energies_.size() != static_cast<std::size_t>(result.beta_mos_.cols()))
        throw std::runtime_error("beta orbital energy count mismatch");
    return result;
}
DensityDerivatives Molecule::density_derivatives(
    const Eigen::Vector3d& point, Spin component) const {
    if (!point.allFinite())
        throw std::invalid_argument("density evaluation point must be finite");
    // P is symmetric. The alpha/beta split preserves the FCHK density
    // convention and permits the same contractions for all four fields.
    const Eigen::MatrixXd* p = &total_;
    Eigen::MatrixXd selected;
    switch (component) {
        case Spin::Total: break;
        case Spin::Difference: p = &spin_; break;
        case Spin::Alpha:
            selected = 0.5 * (total_ + spin_);
            p = &selected;
            break;
        case Spin::Beta:
            selected = 0.5 * (total_ - spin_);
            p = &selected;
            break;
        default: throw std::invalid_argument("invalid density spin component");
    }
    const auto ao = evaluator_.evaluate_derivatives(point);
    const Eigen::VectorXd y = (*p) * ao.values;
    const Eigen::Matrix<double, Eigen::Dynamic, 3> pg = (*p) * ao.gradient;
    DensityDerivatives result;
    result.value = ao.values.dot(y);
    result.gradient = 2.0 * ao.gradient.transpose() * y;
    // AO Hessian storage order: xx, yy, zz, xy, xz, yz.
    constexpr int indices[3][3] = {{0, 3, 4}, {3, 1, 5}, {4, 5, 2}};
    for (int a = 0; a < 3; ++a) {
        for (int b = a; b < 3; ++b) {
            const double value = 2.0 * (
                ao.hessian.col(indices[a][b]).dot(y) +
                ao.gradient.col(a).dot(pg.col(b)));
            result.hessian(a,b) = value;
            result.hessian(b,a) = value;
        }
        result.tau += 0.5 * ao.gradient.col(a).dot(pg.col(a));
    }
    result.laplacian = result.hessian.trace();
    return result;
}
std::vector<DensityDerivatives> Molecule::density_derivatives_blocked(
    const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>>& points,
    Spin component, std::size_t block_size) const {
    if (block_size == 0) throw std::invalid_argument("block size must be positive");
    if (!points.allFinite()) throw std::invalid_argument("density points must be finite");
    Eigen::MatrixXd selected;
    const Eigen::MatrixXd* p = &total_;
    switch (component) {
        case Spin::Total: break;
        case Spin::Difference: p = &spin_; break;
        case Spin::Alpha: selected = 0.5 * (total_ + spin_); p = &selected; break;
        case Spin::Beta: selected = 0.5 * (total_ - spin_); p = &selected; break;
        default: throw std::invalid_argument("invalid density spin component");
    }
    const auto n = static_cast<Eigen::Index>(basis_size());
    const auto count = static_cast<std::size_t>(points.cols());
    std::vector<DensityDerivatives> result(count);
    constexpr int index[3][3] = {{0,3,4},{3,1,5},{4,5,2}};
    for (std::size_t start = 0; start < count;) {
        const std::size_t size = std::min(block_size, count-start);
        const auto b = static_cast<Eigen::Index>(size);
        // Pack AO values and first derivatives into one column-major panel.
        // A single GEMM applies P to all four panels, avoiding four separate
        // matrix products and improving reuse of the (potentially large) P.
        Eigen::MatrixXd input(n, 4*b);
        auto phi = input.middleCols(0, b);
        auto dx  = input.middleCols(b, b);
        auto dy  = input.middleCols(2*b, b);
        auto dz  = input.middleCols(3*b, b);
        Eigen::MatrixXd hxx(n,b), hyy(n,b), hzz(n,b), hxy(n,b), hxz(n,b), hyz(n,b);
        for (Eigen::Index j=0; j<b; ++j) {
            const auto ao = evaluator_.evaluate_derivatives(points.col(static_cast<Eigen::Index>(start)+j));
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
        Eigen::MatrixXd product(n, 4*b);
        product.noalias() = (*p) * input;
        const auto y  = product.middleCols(0, b);
        const Eigen::MatrixXd* h[6] = {&hxx,&hyy,&hzz,&hxy,&hxz,&hyz};
        for (Eigen::Index j=0; j<b; ++j) {
            auto& out = result[start+static_cast<std::size_t>(j)];
            out.value = phi.col(j).dot(y.col(j));
            for (int a=0; a<3; ++a) {
                const auto ga = input.col((a+1)*b+j);
                out.gradient[a] = 2.0*ga.dot(y.col(j));
                out.tau += 0.5*ga.dot(product.col((a+1)*b+j));
                for (int c=a; c<3; ++c) {
                    const double v = 2.0*(h[index[a][c]]->col(j).dot(y.col(j)) +
                                            ga.dot(product.col((c+1)*b+j)));
                    out.hessian(a,c) = out.hessian(c,a) = v;
                }
            }
            out.laplacian = out.hessian.trace();
        }
        start += size;
    }
    return result;
}
Field Molecule::density_pointwise(const Grid& grid, Spin spin) const {
    const auto count=checked_size(grid);
    std::vector<double> values(count);
    Eigen::MatrixXd p;
    switch(spin) {
        case Spin::Total: p=total_; break;
        case Spin::Alpha: p=0.5*(total_+spin_); break;
        case Spin::Beta: p=0.5*(total_-spin_); break;
        case Spin::Difference: p=spin_; break;
    }
    const int ny=grid.counts[1], nz=grid.counts[2];
    for (std::size_t q=0;q<count;++q) {
        const auto k=static_cast<int>(q%nz);
        const auto j=static_cast<int>((q/nz)%ny);
        const auto i=static_cast<int>(q/(static_cast<std::size_t>(ny)*nz));
        Eigen::Vector3d r=grid.origin+grid.steps*Eigen::Vector3d(i,j,k);
        Eigen::VectorXd ao=evaluator_.evaluate(r);
        values[q]=ao.dot(p*ao);
    }
    return Field(grid,atoms_,std::move(values));
}
}

namespace chem::kumo {
Field Molecule::density(const Grid& grid, Spin spin) const {
    auto fields = densities_blocked(grid, {spin});
    return std::move(fields.at(spin));
}

SpinDensities Molecule::spin_densities(const Grid& grid) const {
    auto fields = densities_blocked(grid,
        {Spin::Total, Spin::Alpha, Spin::Beta, Spin::Difference});
    return {std::move(fields.at(Spin::Total)),
            std::move(fields.at(Spin::Alpha)),
            std::move(fields.at(Spin::Beta)),
            std::move(fields.at(Spin::Difference))};
}
}

namespace chem::kumo {
std::map<Spin, Field> Molecule::densities(const Grid& grid,
                                        std::initializer_list<Spin> components) const {
    return densities_blocked(grid, components);
}
}

namespace chem::kumo {
std::map<Spin, Field> Molecule::densities_blocked(
    const Grid& grid, std::initializer_list<Spin> components, std::size_t block_size) const {
    const auto count = checked_size(grid);
    if (block_size == 0 || block_size > 65536)
        throw std::invalid_argument("block_size must be in [1,65536]");
    std::map<Spin, std::vector<double>> buffers;
    for (Spin component : components) buffers.try_emplace(component, count);
    if (buffers.empty()) return {};
    const bool need_a = buffers.contains(Spin::Alpha) || buffers.contains(Spin::Total) || buffers.contains(Spin::Difference);
    const bool need_b = buffers.contains(Spin::Beta) || buffers.contains(Spin::Total) || buffers.contains(Spin::Difference);
    const Eigen::MatrixXd pa = need_a ? Eigen::MatrixXd(0.5*(total_+spin_)) : Eigen::MatrixXd();
    const Eigen::MatrixXd pb = need_b ? Eigen::MatrixXd(0.5*(total_-spin_)) : Eigen::MatrixXd();
    const Eigen::Index m = static_cast<Eigen::Index>(basis_size());
    const int ny = grid.counts[1], nz = grid.counts[2];
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
#pragma omp parallel for schedule(static) if(nblocks > 1)
#endif
    for (std::int64_t block = 0; block < static_cast<std::int64_t>(nblocks); ++block) {
        try {
        const std::size_t start = static_cast<std::size_t>(block) * block_size;
        const auto b = static_cast<Eigen::Index>(std::min(block_size, count-start));
        Eigen::Matrix<double,3,Eigen::Dynamic> points(3,b);
        for (Eigen::Index t=0; t<b; ++t) {
            const std::size_t q=start+static_cast<std::size_t>(t);
            const int k=static_cast<int>(q%nz);
            const int j=static_cast<int>((q/nz)%ny);
            const int i=static_cast<int>(q/(static_cast<std::size_t>(ny)*nz));
            points.col(t)=grid.origin+grid.steps*Eigen::Vector3d(i,j,k);
        }
        Eigen::MatrixXd phi(m,b);
        evaluator_.evaluate_batch(points,phi);
        Eigen::VectorXd a=Eigen::VectorXd::Zero(b), bb=Eigen::VectorXd::Zero(b);
        if (need_a) {
            Eigen::MatrixXd y=pa*phi;
            a=(phi.array()*y.array()).colwise().sum().transpose().matrix();
        }
        if (need_b) {
            Eigen::MatrixXd y=pb*phi;
            bb=(phi.array()*y.array()).colwise().sum().transpose().matrix();
        }
        for (Eigen::Index t=0; t<b; ++t) {
            const auto q=start+static_cast<std::size_t>(t);
            if (auto it=buffers.find(Spin::Alpha);it!=buffers.end()) it->second[q]=a[t];
            if (auto it=buffers.find(Spin::Beta);it!=buffers.end()) it->second[q]=bb[t];
            if (auto it=buffers.find(Spin::Total);it!=buffers.end()) it->second[q]=a[t]+bb[t];
            if (auto it=buffers.find(Spin::Difference);it!=buffers.end()) it->second[q]=a[t]-bb[t];
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
    std::map<Spin,Field> result;
    for (auto& [spin, values]:buffers)
        result.emplace(spin,Field(grid,atoms_,std::move(values)));
    return result;
}
}
