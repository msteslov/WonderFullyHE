#include "m2424/bootstrap_precision_model.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

bool near(double lhs, double rhs, double tolerance) {
    return std::abs(lhs - rhs) <= tolerance;
}

} // namespace

int main() {
    bool ok = true;

    const auto budget = m2424::make_bootstrap_error_budget(1e-9, 3);
    ok = ok && budget.cycles == 3;
    ok = ok && near(budget.per_cycle_budget, 1e-9 / 6.0, 1e-24);
    ok = ok && budget.per_cycle_budget > 1e-10 && budget.per_cycle_budget < 3e-10;
    ok = ok && budget.dft_roundtrip_budget > 0.0;
    ok = ok && budget.evalmod_budget > 0.0;
    ok = ok && budget.modraise_budget > 0.0;
    ok = ok && budget.dft_roundtrip_budget < budget.per_cycle_budget;
    ok = ok && near(budget.dft_roundtrip_budget + budget.evalmod_budget + budget.modraise_budget,
                    budget.per_cycle_budget,
                    1e-24);

    bool threw = false;
    try {
        (void)m2424::make_bootstrap_error_budget(1e-9, 0);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    ok = ok && threw;

    const std::vector<m2424::DftPrecisionMeasurement> current_measurements{
        {"precision_boot_deep_ckks", 4, 50.0, 50.0, 2, 15, 12, 3, 1.97e-11, 3.07e-9, 1.79e-8},
        {"precision_boot_deep_ckks", 4, 50.0, 55.0, 2, 15, 12, 3, 1.10e-11, 1.28e-9, 2.52e-9},
        {"precision_boot_deep_ckks", 4, 50.0, 60.0, 2, 15, 12, 3, 9.58e-12, 1.22e-9, 1.65e-9}
    };
    const auto fit = m2424::fit_dft_precision_floor(current_measurements, 2e-10);
    ok = ok && fit.quantization_coefficient > 0.0;
    ok = ok && fit.noise_floor > 0.0;
    ok = ok && fit.predicted_best_error > 2e-10;
    ok = ok && fit.floor_blocks_target;
    ok = ok && fit.blocker.find("key-switch/rotation-noise-floor") != std::string::npos;

    const auto cost = m2424::estimate_bootstrap_dft_cost(4, m2424::BootstrapDftType::HomomorphicEncode, 60.0);
    ok = ok && cost.slots == 4;
    ok = ok && cost.layer_count > 0;
    ok = ok && cost.diagonal_term_count > 0;
    ok = ok && cost.plaintext_multiplication_count == cost.diagonal_term_count;
    ok = ok && cost.rescale_count == cost.layer_count;

    m2424::BootstrapPrecisionPlanningRequest request;
    request.target_total_error = 1e-9;
    request.cycles = 3;
    request.slots = 4;
    request.candidate_profiles.push_back(m2424::CkksProfile{32768, {60, 50, 60}, std::exp2(50.0), 4});
    request.candidate_transform_scales = {50.0, 55.0, 60.0};
    const auto planning = m2424::plan_bootstrap_precision(request, current_measurements);
    ok = ok && !planning.feasible;
    ok = ok && planning.dft_fit.floor_blocks_target;
    ok = ok && !planning.blocker.empty();

    const m2424::ComplexVector sample{
        {1e-5, -2e-6},
        {-3e-6, 4e-6},
        {2e-6, 1e-6},
        {-1e-6, -3e-6}
    };
    const auto small_stc = m2424::FactorizedLinearTransform(m2424::make_small_slots4_stc_plan(60.0));
    const auto small_cts = m2424::FactorizedLinearTransform(m2424::make_small_slots4_cts_plan(60.0));
    const auto small_roundtrip = small_cts.apply_plain(small_stc.apply_plain(sample));
    double small_plain_error = 0.0;
    for (std::size_t i = 0; i < sample.size(); ++i) {
        small_plain_error = std::max(small_plain_error, std::abs(small_roundtrip[i] - sample[i]));
    }
    ok = ok && small_plain_error <= 1e-12;

    if (!ok) {
        std::cerr << "bootstrap precision model test failed\n";
        return 1;
    }
    return 0;
}
