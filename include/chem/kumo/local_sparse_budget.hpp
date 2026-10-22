#pragma once

#include <Eigen/Core>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace chem::kumo::detail {

struct LocalBudgetCharge {
    Eigen::Index row;
    Eigen::Index col;
    std::size_t channel;
    double amount;
};

// Each charge is a nonnegative contribution to an upper-triangle AO element.
// Repeated coordinates are allowed: they are accumulated before admission.
class LocalSparseBudget {
public:
    LocalSparseBudget(Eigen::Index nao, std::vector<double> budgets)
        : spent_(budgets.size(), Eigen::MatrixXd::Zero(nao, nao)),
          budgets_(std::move(budgets)), nao_(nao) {
        if (nao < 0 || budgets_.empty())
            throw std::invalid_argument("invalid budget dimensions");
        for (double b : budgets_)
            if (!(b >= 0.0) || !std::isfinite(b))
                throw std::invalid_argument("invalid channel budget");
    }

    // Transactional: rejected quartets leave all channel ledgers unchanged.
    bool try_charge(const std::vector<LocalBudgetCharge>& charges) {
        for (const auto& c : charges) validate(c);
        // Coalesce duplicates; no full AO matrix traversal or temporary dense matrix.
        std::vector<LocalBudgetCharge> merged;
        merged.reserve(charges.size());
        for (const auto& c : charges) {
            bool found = false;
            for (auto& m : merged) {
                if (m.channel == c.channel && m.row == c.row && m.col == c.col) {
                    const long double sum = static_cast<long double>(m.amount) + c.amount;
                    if (sum > std::numeric_limits<double>::max()) return false;
                    m.amount = static_cast<double>(sum);
                    if (m.amount < sum)
                        m.amount = std::nextafter(m.amount, std::numeric_limits<double>::infinity());
                    found = true;
                    break;
                }
            }
            if (!found) merged.push_back(c);
        }
        for (const auto& c : merged) {
            const long double proposed = static_cast<long double>(spent_[c.channel](c.row,c.col)) + c.amount;
            if (proposed > budgets_[c.channel]) return false;
        }
        for (const auto& c : merged) {
            const long double sum = static_cast<long double>(spent_[c.channel](c.row,c.col)) + c.amount;
            double rounded = static_cast<double>(sum);
            if (rounded < sum)
                rounded = std::nextafter(rounded, std::numeric_limits<double>::infinity());
            spent_[c.channel](c.row,c.col) = rounded;
        }
        return true;
    }

    const std::vector<Eigen::MatrixXd>& spent() const noexcept { return spent_; }

private:
    void validate(const LocalBudgetCharge& c) const {
        if (c.channel >= budgets_.size() || c.row < 0 || c.col < c.row || c.col >= nao_ ||
            !(c.amount >= 0.0) || !std::isfinite(c.amount))
            throw std::invalid_argument("invalid local budget charge");
    }
    std::vector<Eigen::MatrixXd> spent_;
    std::vector<double> budgets_;
    Eigen::Index nao_;
};

} // namespace chem::kumo::detail
