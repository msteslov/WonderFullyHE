#include "m2424/bootstrap_plan.hpp"

#include <cmath>
#include <cstdio>
#include <vector>

namespace {

m2424::CipherInfo make_info(double scale_log2,
                            std::size_t chain_index,
                            std::size_t coeff_modulus_size,
                            double coeff_modulus_log2) {
    m2424::CipherInfo info;
    info.scale = std::exp2(scale_log2);
    info.chain_index = chain_index;
    info.coeff_modulus_size = coeff_modulus_size;
    info.coeff_modulus_log2 = coeff_modulus_log2;
    return info;
}

m2424::BootstrapScaleDesign design(double period_log2,
                                   double plain_scale_log2,
                                   double target_scale_log2,
                                   const std::vector<int>& active_bits,
                                   const m2424::CipherInfo& info,
                                   double max_abs_before_normalization,
                                   std::size_t min_chain_remaining) {
    return m2424::make_bootstrap_scale_design(
        m2424::BootstrapPeriodMode::ManualPowerOfTwo,
        period_log2,
        period_log2,
        m2424::BootstrapScalingStrategy::DecomposedPlainMultiplyRescale,
        plain_scale_log2,
        target_scale_log2,
        m2424::EvalModDegree::P3,
        active_bits,
        info,
        max_abs_before_normalization,
        min_chain_remaining,
        2.0);
}

} // namespace

int main() {
    const auto period_blocked = design(0.0,
                                       40.0,
                                       60.0,
                                       {60, 40, 40, 60},
                                       make_info(40.0, 3, 4, 200.0),
                                       1.0,
                                       1);

    const auto scale_blocked = design(220.0,
                                      50.0,
                                      60.0,
                                      {60, 40},
                                      make_info(40.0, 1, 2, 100.0),
                                      1.0e60,
                                      1);

    const auto capacity_blocked = design(0.0,
                                         40.0,
                                         40.0,
                                         {1, 1, 1, 1, 1, 1},
                                         make_info(40.0, 5, 6, 6.0),
                                         1.0e-5,
                                         1);

    const auto ready = design(20.0,
                              40.0,
                              60.0,
                              {60, 40, 40, 40, 40, 40, 60},
                              make_info(40.0, 6, 7, 320.0),
                              1.0e-1,
                              1);

    const bool ok =
        period_blocked.status == m2424::BootstrapScaleDesignStatus::PeriodModelBlocked
        && scale_blocked.status == m2424::BootstrapScaleDesignStatus::ScaleStrategyBlocked
        && capacity_blocked.status == m2424::BootstrapScaleDesignStatus::EvalModCapacityBlocked
        && ready.status == m2424::BootstrapScaleDesignStatus::ReadyForEvalModP3
        && ready.evalmod_capacity_ok
        && ready.scale_strategy_ok
        && ready.magnitude_ok;

    std::printf("[test_bootstrap_scale_design] status_cases=%s period=%s scale=%s capacity=%s ready=%s\n",
                ok ? "PASS" : "FAIL",
                m2424::to_string(period_blocked.status),
                m2424::to_string(scale_blocked.status),
                m2424::to_string(capacity_blocked.status),
                m2424::to_string(ready.status));
    return ok ? 0 : 1;
}
