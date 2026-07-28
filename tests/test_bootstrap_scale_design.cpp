#include "m2424/bootstrap_plan.hpp"

#include <cmath>
#include <cstdio>
#include <vector>

namespace {

m2424::CipherInfo make_info(double scale_log2,
                            std::size_t chainIndex,
                            std::size_t coeffModulusSize,
                            double coeffModulusLog2) {
    m2424::CipherInfo info;
    info.scale = std::exp2(scale_log2);
    info.chainIndex = chainIndex;
    info.coeffModulusSize = coeffModulusSize;
    info.coeffModulusLog2 = coeffModulusLog2;
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

m2424::CipherInfo refresh_info(std::size_t chainIndex) {
    return make_info(45.0, chainIndex, chainIndex + 1, 45.0 * static_cast<double>(chainIndex + 1));
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
    const auto impossible_window = m2424::bootstrap_period_feasibility_window(
        300.0,
        40.0,
        60.0,
        std::exp2(242.0),
        m2424::EvalModPolynomial::approximation_bound,
        2.0);
    const auto possible_window = m2424::bootstrap_period_feasibility_window(
        360.0,
        40.0,
        60.0,
        std::exp2(242.0),
        m2424::EvalModPolynomial::approximation_bound,
        2.0);
    m2424::CkksProfile profile;
    profile.coeffModulusBits = {60, 40, 40, 60};
    const auto active_bits = m2424::active_coeff_modulus_bits(
        profile,
        make_info(40.0, 2, 3, 140.0));
    const auto refresh_scale_gate = m2424::plan_bootstrap_refresh_scale_gate({
        profile,
        make_info(40.0, 1, 2, 100.0),
        make_info(40.0, 3, 4, 200.0),
        make_info(40.0, 3, 4, 200.0),
        1.0e-5
    });
    const auto refresh_scale_search = m2424::search_bootstrap_refresh_scale_gate({
        profile,
        make_info(40.0, 1, 2, 100.0),
        make_info(40.0, 3, 4, 200.0),
        make_info(40.0, 3, 4, 200.0),
        1.0e-5,
        {
            m2424::BootstrapPeriodMode::TotalCoeffModulus,
            m2424::BootstrapPeriodMode::ManualPowerOfTwo
        },
        {20.0, 40.0, 160.0},
        {0.0},
        {},
        0.0,
        {40.0, 80.0},
        {40.0, 60.0}
    });
    const auto dense_layout = m2424::plan_bootstrap_layout({
        16,
        16384,
        m2424::SecurityLevel::TC128,
        m2424::BootstrapTransformBackend::DenseDiagonal,
        252.0,
        40.0,
        160.0,
        60.0,
        2.0,
        1,
        3,
        1,
        1,
        60,
        40,
        60
    });
    const auto fft_layout = m2424::plan_bootstrap_layout({
        16,
        32768,
        m2424::SecurityLevel::TC128,
        m2424::BootstrapTransformBackend::FftLike,
        132.0,
        40.0,
        160.0,
        60.0,
        2.0
    });
    const auto lattigo_like_mod1_levels = m2424::estimated_bootstrap_mod1_levels({
        m2424::BootstrapMod1Type::CosDiscrete,
        31,
        3
    });
    const auto openfhe_like_mod1_levels = m2424::estimated_bootstrap_mod1_levels({
        m2424::BootstrapMod1Type::CosDiscrete,
        91,
        6
    });
    const auto cos_layout = m2424::plan_bootstrap_layout({
        16,
        32768,
        m2424::SecurityLevel::TC128,
        m2424::BootstrapTransformBackend::FftLike,
        132.0,
        40.0,
        160.0,
        60.0,
        2.0,
        0,
        0,
        0,
        1,
        60,
        40,
        60,
        {m2424::BootstrapMod1Type::CosDiscrete, 31, 3}
    });

    const auto parameter_ready = m2424::plan_bootstrap_parameters({
        16,
        32768,
        m2424::SecurityLevel::TC128,
        m2424::BootstrapCircuitOrder::SlotsToCoeffsFirst,
        m2424::BootstrapTransformBackend::FftLike,
        0.0,
        40.0,
        40.0,
        40.0,
        2.0,
        1,
        60,
        40,
        60,
        {m2424::BootstrapMod1Type::CosDiscrete, 31, 3}
    });
    const auto parameter_unsupported = m2424::plan_bootstrap_parameters({
        32,
        32768,
        m2424::SecurityLevel::TC128,
        m2424::BootstrapCircuitOrder::SlotsToCoeffsFirst,
        m2424::BootstrapTransformBackend::FftLike,
        0.0,
        40.0,
        40.0,
        40.0,
        2.0,
        1,
        60,
        40,
        60,
        {m2424::BootstrapMod1Type::CosDiscrete, 31, 3}
    });
    const auto parameter_security_blocked = m2424::plan_bootstrap_parameters({
        16,
        16384,
        m2424::SecurityLevel::TC128,
        m2424::BootstrapCircuitOrder::SlotsToCoeffsFirst,
        m2424::BootstrapTransformBackend::DenseDiagonal,
        252.0,
        40.0,
        160.0,
        60.0,
        2.0,
        1,
        60,
        40,
        60,
        {m2424::BootstrapMod1Type::CosDiscrete, 31, 3}
    });

    m2424::CkksOperationBudget small_budget;
    small_budget.ciphertext_muls = 2;
    const auto fits_refresh_plan = m2424::plan_bootstrap_refresh({
        refresh_info(3),
        small_budget,
        1e-9,
        4096,
        128,
        m2424::ParameterOptimizeFor::Speed,
        1
    });

    m2424::CkksOperationBudget large_budget;
    large_budget.ciphertext_muls = 2;
    const auto required_refresh_plan = m2424::plan_bootstrap_refresh({
        refresh_info(1),
        large_budget,
        1e-9,
        4096,
        128,
        m2424::ParameterOptimizeFor::Speed,
        0
    });

    m2424::CkksOperationBudget impossible_budget;
    impossible_budget.evalmod_p3 = 10;
    const auto blocked_refresh_plan = m2424::plan_bootstrap_refresh({
        refresh_info(1),
        impossible_budget,
        1e-15,
        4096,
        128,
        m2424::ParameterOptimizeFor::Speed,
        0
    });

    const bool ok =
        period_blocked.status == m2424::BootstrapScaleDesignStatus::PeriodModelBlocked
        && scale_blocked.status == m2424::BootstrapScaleDesignStatus::ScaleStrategyBlocked
        && !scale_blocked.scale_plan.blocker.empty()
        && capacity_blocked.status == m2424::BootstrapScaleDesignStatus::EvalModCapacityBlocked
        && ready.status == m2424::BootstrapScaleDesignStatus::ReadyForEvalModP3
        && ready.evalmod_capacity_ok
        && ready.scale_strategy_ok
        && ready.magnitude_ok
        && !impossible_window.possible
        && impossible_window.margin_log2 < 0.0
        && possible_window.possible
        && possible_window.margin_log2 >= 0.0
        && active_bits == std::vector<int>({60, 40, 40})
        && refresh_scale_gate.status == m2424::BootstrapScaleDesignStatus::ScaleStrategyBlocked
        && refresh_scale_gate.period_log2 == 160.0
        && refresh_scale_search.candidates == 16
        && refresh_scale_search.ready == (refresh_scale_search.ready_candidates > 0)
        && refresh_scale_search.best_design.period_log2 > 0.0
        && dense_layout.total_levels > dense_layout.normalization_levels
        && dense_layout.security_budget_bits == 438
        && fft_layout.coeff_to_slot_levels > dense_layout.coeff_to_slot_levels
        && fft_layout.profile.slots == fft_layout.polyModulusDegree / 2
        && !fft_layout.blocker.empty()
        && lattigo_like_mod1_levels == 8
        && openfhe_like_mod1_levels == 13
        && cos_layout.evalmod_levels == lattigo_like_mod1_levels
        && cos_layout.mod1_model.type == m2424::BootstrapMod1Type::CosDiscrete
        && parameter_ready.status == m2424::BootstrapParameterPlanningStatus::Ready
        && parameter_ready.blocker == "none"
        && parameter_ready.expected_output_chain_index == 1
        && !parameter_ready.rotationSteps.empty()
        && parameter_unsupported.status == m2424::BootstrapParameterPlanningStatus::BlockedByUnsupportedSlots
        && parameter_security_blocked.status == m2424::BootstrapParameterPlanningStatus::BlockedBySecurityBudget
        && fits_refresh_plan.status == m2424::BootstrapRefreshPlanningStatus::ComputeFitsWithoutRefresh
        && !fits_refresh_plan.needs_refresh
        && required_refresh_plan.status == m2424::BootstrapRefreshPlanningStatus::RefreshRequired
        && required_refresh_plan.needs_refresh
        && blocked_refresh_plan.status == m2424::BootstrapRefreshPlanningStatus::RefreshPlanBlocked
        && blocked_refresh_plan.needs_refresh;

    std::printf("[test_bootstrap_scale_design] status_cases=%s period=%s scale=%s capacity=%s ready=%s params=%s\n",
                ok ? "PASS" : "FAIL",
                m2424::to_string(period_blocked.status),
                m2424::to_string(scale_blocked.status),
                m2424::to_string(capacity_blocked.status),
                m2424::to_string(ready.status),
                m2424::to_string(parameter_ready.status));
    return ok ? 0 : 1;
}
