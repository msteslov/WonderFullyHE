#include "m2424/bootstrap_candidates.hpp"

#include <cmath>
#include <cstdio>
#include <stdexcept>

int main() {
    const auto& candidates = m2424::bootstrapCandidates();
    bool allValid = candidates.size() >= 12;
    for (const auto& candidate : candidates) {
        allValid = allValid && m2424::isBootstrapCandidateValid(candidate);
    }

    const auto& selected = m2424::bootstrapCandidateById("balanced_8192_s55");
    const auto incomplete = m2424::forecastBootstrapFeasibility(selected, {});
    const auto feasible = m2424::forecastBootstrapFeasibility(selected, {
        1e-11, 1e-11, 5e-11, 1e-10, 5e-11
    });
    const auto infeasible = m2424::forecastBootstrapFeasibility(selected, {
        1e-11, 1e-11, 3e-10, 1e-10, 5e-11
    });
    const double propagated = m2424::propagatedBootstrapError({
        1e-11, 2e-11, 3e-11, 4e-11, 5e-11, 2.0, 3.0, 6e-11, 7e-11
    });

    bool unknownRejected = false;
    try {
        (void)m2424::bootstrapCandidateById("missing");
    } catch (const std::invalid_argument&) {
        unknownRejected = true;
    }

    const bool ok = allValid && !incomplete.complete && !incomplete.passes
        && feasible.complete && feasible.passes && feasible.predictedTotalError <= m2424::kTargetAbsoluteError
        && infeasible.complete && !infeasible.passes && unknownRejected
        && std::abs(propagated - 6.1e-10) < 1e-22;
    std::printf("[test_bootstrap_candidates] count=%zu %s\n", candidates.size(), ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
