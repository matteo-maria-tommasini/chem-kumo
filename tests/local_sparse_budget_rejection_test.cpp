#include "chem/kumo/local_sparse_budget.hpp"
#include <Eigen/Core>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

using chem::kumo::detail::LocalBudgetCharge;
using chem::kumo::detail::LocalSparseBudget;

static void require(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}

static void unchanged(const LocalSparseBudget& ledger,
                      const std::vector<Eigen::MatrixXd>& before) {
    const auto& after = ledger.spent();
    require(after.size() == before.size(), "channel count changed");
    for (std::size_t k = 0; k < before.size(); ++k)
        require((after[k].array() == before[k].array()).all(),
                "rejected transaction changed the ledger");
}

int main() {
    constexpr Eigen::Index nao = 1410;
    LocalSparseBudget ledger(nao, {1.0, 0.5});
    std::size_t checks = 0, accepted = 0, rejected = 0;

    auto accept = [&](const std::vector<LocalBudgetCharge>& charges) {
        require(ledger.try_charge(charges), "expected acceptance");
        ++checks;
        ++accepted;
    };
    auto reject = [&](const std::vector<LocalBudgetCharge>& charges) {
        const auto before = ledger.spent();
        require(!ledger.try_charge(charges), "expected rejection");
        unchanged(ledger, before);
        ++checks;
        ++rejected;
    };

    // Exact equality to the budget must be accepted.
    accept({{0, 0, 0, 1.0}});
    require(ledger.spent()[0](0, 0) == 1.0, "boundary charge mismatch");

    // One over-budget charge must be rejected without changing state.
    reject({{0, 0, 0, 0.25}});
    reject({{1, 1, 1, 0.75}});

    // A valid charge preceding an invalid charge cannot be partially committed.
    reject({{2, 3, 0, 0.25}, {1, 1, 1, 0.75}});
    require(ledger.spent()[0](2, 3) == 0.0, "partial commit");

    // Individually admissible duplicates exceed the budget when combined.
    reject({{4, 5, 1, 0.3}, {4, 5, 1, 0.3}});

    // Same AO coordinate in two channels: failure in one must roll back both.
    reject({{6, 7, 0, 0.25}, {6, 7, 1, 0.6}});

    // A rejected transaction must not poison subsequent valid updates.
    accept({{4, 5, 1, 0.25}, {4, 5, 1, 0.25}, {6, 7, 0, 0.125}});
    require(ledger.spent()[1](4, 5) == 0.5, "duplicate sum mismatch");
    require(ledger.spent()[0](6, 7) == 0.125, "recovery mismatch");
    reject({{4, 5, 1, 0.125}});

    // Empty transaction and zero charge are valid no-ops.
    const auto snapshot = ledger.spent();
    accept({});
    accept({{8, 9, 0, 0.0}});
    unchanged(ledger, snapshot);

    // Invalid input must throw and preserve all ledgers.
    const std::vector<std::vector<LocalBudgetCharge>> invalid = {
        {{-1, 0, 0, 0.1}},
        {{2, 1, 0, 0.1}},
        {{0, nao, 0, 0.1}},
        {{0, 0, 2, 0.1}},
        {{0, 0, 0, -0.1}},
        {{0, 0, 0, std::numeric_limits<double>::quiet_NaN()}},
        {{0, 0, 0, std::numeric_limits<double>::infinity()}},
        {{10, 11, 0, 0.1}, {0, 0, 0, std::numeric_limits<double>::quiet_NaN()}}
    };
    for (const auto& charges : invalid) {
        const auto before = ledger.spent();
        bool threw = false;
        try { (void)ledger.try_charge(charges); }
        catch (const std::invalid_argument&) { threw = true; }
        require(threw, "invalid transaction was not rejected");
        unchanged(ledger, before);
        ++checks;
    }

    // Overflow while merging duplicate amounts must reject atomically.
    reject({{12, 13, 0, std::numeric_limits<double>::max()},
            {12, 13, 0, std::numeric_limits<double>::max()}});

    std::cout << "checks=" << checks
              << " accepted=" << accepted
              << " rejected=" << rejected
              << " invalid=" << invalid.size()
              << " channels=2 nao=" << nao
              << " atomicity=OK status=OK\n";
}
