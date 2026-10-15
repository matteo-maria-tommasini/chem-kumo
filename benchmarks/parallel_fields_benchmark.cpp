#include "chem/kumo.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace {
using Clock = std::chrono::steady_clock;
using chem::kumo::Molecule;
using chem::kumo::Grid;
using chem::kumo::Spin;

Grid make_grid(const std::vector<chem::kumo::Atom>& atoms, std::size_t requested, double margin) {
    if (atoms.empty()) throw std::invalid_argument("empty molecule");
    const double side = std::ceil(std::cbrt(static_cast<double>(requested)));
    if (side < 2 || side > 1000) throw std::invalid_argument("invalid grid side");
    const int n = static_cast<int>(side);
    Grid grid;
    grid.counts = Eigen::Vector3i(n,n,n);
    Eigen::Vector3d lo = atoms.front().position, hi = lo;
    for (const auto& atom : atoms) {
        lo = lo.cwiseMin(atom.position);
        hi = hi.cwiseMax(atom.position);
    }
    lo.array() -= margin;
    hi.array() += margin;
    grid.origin = lo;
    grid.steps.setZero();
    for (int k=0;k<3;++k) grid.steps(k,k) = (hi[k]-lo[k])/(n-1);
    return grid;
}

Eigen::Matrix<double,3,Eigen::Dynamic> grid_points(const Grid& grid) {
    const std::size_t nx=static_cast<std::size_t>(grid.counts[0]);
    const std::size_t ny=static_cast<std::size_t>(grid.counts[1]);
    const std::size_t nz=static_cast<std::size_t>(grid.counts[2]);
    const std::size_t count=nx*ny*nz;
    Eigen::Matrix<double,3,Eigen::Dynamic> points(3,static_cast<Eigen::Index>(count));
    for (std::size_t q=0;q<count;++q) {
        const auto k=static_cast<int>(q%nz);
        const auto j=static_cast<int>((q/nz)%ny);
        const auto i=static_cast<int>(q/(ny*nz));
        points.col(static_cast<Eigen::Index>(q))=grid.origin+grid.steps*Eigen::Vector3d(i,j,k);
    }
    return points;
}

// All reported components enter the comparison; never compare only a checksum.
void append(std::vector<double>& dest, const chem::kumo::DensityDerivatives& x) {
    dest.push_back(x.value);
    for (int k=0;k<3;++k) dest.push_back(x.gradient[k]);
    for (int i=0;i<3;++i) for (int j=0;j<3;++j) dest.push_back(x.hessian(i,j));
    dest.push_back(x.laplacian);
    dest.push_back(x.tau);
}
void append(std::vector<double>& dest, const chem::kumo::OrbitalDerivatives& x) {
    dest.push_back(x.value);
    for (int k=0;k<3;++k) dest.push_back(x.gradient[k]);
    for (int i=0;i<3;++i) for (int j=0;j<3;++j) dest.push_back(x.hessian(i,j));
}

std::vector<double> evaluate(const Molecule& mol, const Grid& grid,
    const Eigen::Matrix<double,3,Eigen::Dynamic>& points,
    const std::string& mode, std::size_t block, std::size_t mo_index,
    const std::vector<std::size_t>& mo_indices) {
    std::vector<double> out;
    if (mode == "density-derivatives") {
        const auto values = mol.density_derivatives_blocked(points, Spin::Total, block);
        out.reserve(values.size()*15);
        for (const auto& x: values) append(out,x);
    } else if (mode == "mo-derivatives") {
        const auto values = mol.orbital_derivatives_blocked(points, mo_index, Spin::Alpha, block);
        out.reserve(values.size()*13);
        for (const auto& x: values) append(out,x);
    } else if (mode == "mo-single" || mode == "mo-multiple") {
        const std::vector<std::size_t> indices = mode == "mo-single"
            ? std::vector<std::size_t>{mo_index} : mo_indices;
        const auto fields = mol.orbitals(grid,indices,Spin::Alpha);
        out.reserve(points.cols()*static_cast<Eigen::Index>(indices.size()));
        for (const auto& field:fields) {
            const auto values=field.values();
            out.insert(out.end(),values.begin(),values.end());
        }
    } else {
        throw std::invalid_argument("unknown mode: " + mode);
    }
    return out;
}

double max_abs_difference(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size()!=b.size()) throw std::runtime_error("output size mismatch");
    double error=0;
    for (std::size_t i=0;i<a.size();++i) {
        if (!std::isfinite(a[i]) || !std::isfinite(b[i]))
            throw std::runtime_error("non-finite field value");
        error=std::max(error,std::abs(a[i]-b[i]));
    }
    return error;
}

void set_threads(int n) {
#ifdef _OPENMP
    omp_set_dynamic(0);
    omp_set_num_threads(n);
#else
    if (n!=1) throw std::runtime_error("OpenMP unavailable");
#endif
}

template<class Work> double best_time(int repeats, Work&& work) {
    work(); // warm-up excluded from timing
    double best=std::numeric_limits<double>::infinity();
    for (int r=0;r<repeats;++r) {
        const auto start=Clock::now();
        work();
        best=std::min(best,std::chrono::duration<double>(Clock::now()-start).count());
    }
    return best;
}
}

int main(int argc, char** argv) {
    if (argc < 3 || argc > 7) {
        std::cerr << "Usage: kumo-parallel-fields-benchmark FILE.fchk[.gz] "
                     "MODE [points=10000] [repeats=3] [block=128] [margin_bohr=2.5]\n"
                     "MODE: density-derivatives | mo-single | mo-multiple | mo-derivatives\n";
        return 2;
    }
    try {
        const std::string mode=argv[2];
        const auto requested=argc>3 ? std::stoull(argv[3]) : 10000ULL;
        const int repeats=argc>4 ? std::stoi(argv[4]) : 3;
        const auto block=argc>5 ? std::stoull(argv[5]) : 128ULL;
        const double margin=argc>6 ? std::stod(argv[6]) : 2.5;
        if (requested<8 || requested>1000000 || repeats<1 || repeats>20 ||
            block<1 || block>65536 || !std::isfinite(margin) || margin<=0)
            throw std::invalid_argument("invalid benchmark arguments");
        const auto mol=Molecule::from_fchk(argv[1]);
        const auto grid=make_grid(mol.atoms(),requested,margin);
        const auto points=grid_points(grid);
        const auto count=static_cast<std::size_t>(points.cols());
        const auto occupied=static_cast<std::size_t>(std::max(1,mol.occupied_count(Spin::Alpha)));
        const auto nmo=mol.orbital_count(Spin::Alpha);
        if (nmo==0) throw std::runtime_error("no alpha orbitals");
        const auto mo_index=std::min(occupied,nmo);
        std::vector<std::size_t> mo_indices;
        for (std::size_t i=1;i<=std::min<std::size_t>(8,nmo);++i) mo_indices.push_back(i);
        std::cout << std::setprecision(12) << "file=" << argv[1] << " mode=" << mode
                  << " basis=" << mol.basis_size() << " atoms=" << mol.atoms().size()
                  << " points=" << count << " repeats=" << repeats
                  << " block_size=" << block << " margin_bohr=" << margin
                  << " mo_index=" << mo_index << " mo_count=" << mo_indices.size() << '\n';
#ifdef _OPENMP
        std::cout << "openmp=enabled\n";
#else
        std::cout << "openmp=disabled\n";
#endif
        set_threads(1);
        const auto reference=evaluate(mol,grid,points,mode,block,mo_index,mo_indices);
        std::vector<double> current;
        auto work=[&] { current=evaluate(mol,grid,points,mode,block,mo_index,mo_indices); };
        const double baseline=best_time(repeats,work);
        const double first_error=max_abs_difference(reference,current);
        if (first_error>1e-10) throw std::runtime_error("serial mismatch");
        std::cout << "threads=1 time_s=" << baseline << " speedup=1 efficiency=1 max_abs=" << first_error << '\n';
#ifdef _OPENMP
        for (int threads : {2,4,8,16}) {
            set_threads(threads);
            const double seconds=best_time(repeats,work);
            const double error=max_abs_difference(reference,current);
            if (error>1e-10) throw std::runtime_error("parallel mismatch");
            const double speedup=baseline/seconds;
            std::cout << "threads=" << threads << " time_s=" << seconds
                      << " speedup=" << speedup << " efficiency=" << speedup/threads
                      << " max_abs=" << error << '\n';
        }
#endif
        std::cout << "samples=" << reference.size() << " checksum=" << reference.at(reference.size()/2) << '\n';
    } catch (const std::exception& e) {
        std::cerr << "benchmark error: " << e.what() << '\n';
        return 1;
    }
}
