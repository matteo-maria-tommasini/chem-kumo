#include "chem/kumo/local_sparse_budget.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

int main() {
    using chem::kumo::detail::LocalBudgetCharge;
    using chem::kumo::detail::LocalSparseBudget;
    constexpr Eigen::Index n=1410;
    const std::vector<double> budgets{1e-8,1e-10};
    LocalSparseBudget local(n,budgets);
    std::vector<Eigen::MatrixXd> dense(2,Eigen::MatrixXd::Zero(n,n));
    std::size_t checks=0, accepted=0, rejected=0;
    for (std::size_t trial=0;trial<300;++trial) {
        std::vector<LocalBudgetCharge> charges;
        for (std::size_t j=0;j<24;++j) {
            const Eigen::Index a=static_cast<Eigen::Index>((trial*17+j*23)%n);
            const Eigen::Index b=static_cast<Eigen::Index>((trial*29+j*31)%n);
            const auto lo=std::min(a,b), hi=std::max(a,b);
            const std::size_t k=j%2;
            const double value=(trial%9==0 ? 1e-9 : 1e-12)*(k==0 ? 1.0 : 0.01);
            charges.push_back({lo,hi,k,value});
            if (j%5==0) charges.push_back({lo,hi,k,value});
        }
        auto proposed=dense;
        for (const auto& c:charges) proposed[c.channel](c.row,c.col)+=c.amount;
        bool dense_ok=true;
        for (std::size_t k=0;k<budgets.size();++k)
            for (Eigen::Index i=0;i<n;++i)
                for (Eigen::Index j=i;j<n;++j)
                    if (proposed[k](i,j)>budgets[k]) dense_ok=false;
        const bool sparse_ok=local.try_charge(charges);
        if (sparse_ok!=dense_ok) throw std::runtime_error("dense/local admission mismatch");
        if (sparse_ok) { dense=std::move(proposed); ++accepted; }
        else ++rejected;
        ++checks;
    }
    for (std::size_t k=0;k<budgets.size();++k)
        if ((dense[k]-local.spent()[k]).cwiseAbs().maxCoeff()>1e-22)
            throw std::runtime_error("ledger mismatch");
    bool invalid_rejected=false;
    try { (void)local.try_charge({{0,0,0,std::numeric_limits<double>::quiet_NaN()}}); }
    catch (const std::invalid_argument&) { invalid_rejected=true; }
    if (!invalid_rejected) throw std::runtime_error("invalid value accepted");
    std::cout<<"checks="<<checks<<" accepted="<<accepted<<" rejected="<<rejected
             <<" channels=2 nao="<<n<<" status=OK\n";
}
