#include "m2424/bootstrap_scaling.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace m2424 {

const char* to_string(BootstrapPeriodMode mode) noexcept {
    switch (mode) {
    case BootstrapPeriodMode::NoBootstrapPeriod:
        return "NoBootstrapPeriod";
    case BootstrapPeriodMode::SourceCoeffModulus:
        return "SourceCoeffModulus";
    case BootstrapPeriodMode::TotalCoeffModulus:
        return "TotalCoeffModulus";
    case BootstrapPeriodMode::LastPrime:
        return "LastPrime";
    case BootstrapPeriodMode::DroppedPrimeProduct:
        return "DroppedPrimeProduct";
    case BootstrapPeriodMode::ManualPowerOfTwo:
        return "ManualPowerOfTwo";
    }
    return "unknown";
}

double bootstrapPeriodLog2(BootstrapPeriodMode mode,
                             double manual_period_log2,
                             const std::vector<int>& coeffModulusBits,
                             const CipherInfo& before_mod_raise,
                             const CipherInfo& after_mod_raise) {
    switch (mode) {
    case BootstrapPeriodMode::NoBootstrapPeriod:
        return 0.0;
    case BootstrapPeriodMode::SourceCoeffModulus:
        if (!std::isfinite(before_mod_raise.scale) || before_mod_raise.scale <= 0.0) {
            throw std::invalid_argument("source ciphertext scale must be positive and finite");
        }
        return before_mod_raise.coeffModulusLog2 - std::log2(before_mod_raise.scale);
    case BootstrapPeriodMode::TotalCoeffModulus:
        return after_mod_raise.coeffModulusLog2 - std::log2(after_mod_raise.scale);
    case BootstrapPeriodMode::LastPrime: {
        if (coeffModulusBits.empty()) {
            throw std::invalid_argument("coeffModulusBits must not be empty");
        }
        if (before_mod_raise.coeffModulusSize >= coeffModulusBits.size()) {
            return static_cast<double>(coeffModulusBits.back());
        }
        return static_cast<double>(coeffModulusBits[before_mod_raise.coeffModulusSize]);
    }
    case BootstrapPeriodMode::DroppedPrimeProduct:
        return std::max(0.0, after_mod_raise.coeffModulusLog2 - before_mod_raise.coeffModulusLog2);
    case BootstrapPeriodMode::ManualPowerOfTwo:
        return manual_period_log2;
    }
    return 0.0;
}

BootstrapScalingFactors make_bootstrap_scaling_factors(double amplitude_factor,
                                                       double bootstrapPeriodLog2,
                                                       double plain_scale_log2) {
    if (!std::isfinite(amplitude_factor) || amplitude_factor <= 0.0) {
        throw std::invalid_argument("amplitude factor must be positive and finite");
    }
    if (!std::isfinite(bootstrapPeriodLog2)) {
        throw std::invalid_argument("bootstrap period log2 must be finite");
    }
    if (!std::isfinite(plain_scale_log2)) {
        throw std::invalid_argument("plain scale log2 must be finite");
    }

    BootstrapScalingFactors factors;
    factors.bootstrapPeriodLog2 = bootstrapPeriodLog2;
    factors.normalization_factor_log2 = bootstrapPeriodLog2 > 0.0
        ? -bootstrapPeriodLog2
        : std::log2(amplitude_factor);
    factors.plain_scale_log2 = plain_scale_log2;
    factors.factor_times_plain_scale_log2 =
        factors.normalization_factor_log2 + factors.plain_scale_log2;
    factors.representable = factors.factor_times_plain_scale_log2 >= 0.0;
    factors.factor = std::exp2(factors.normalization_factor_log2);
    return factors;
}

BootstrapScalarApplication apply_bootstrap_scalar_decomposed(SealAdapter& adapter,
                                                             const Cipher& input,
                                                             double factor_log2,
                                                             double plain_scale_log2) {
    if (!std::isfinite(factor_log2)) {
        throw std::invalid_argument("bootstrap scalar factor log2 must be finite");
    }
    if (!std::isfinite(plain_scale_log2) || plain_scale_log2 <= 0.0) {
        throw std::invalid_argument("plain scale log2 must be positive and finite");
    }

    const auto start_info = adapter.info(input);
    BootstrapScalarApplication application;
    constexpr double scale_capacity_margin_log2 = 2.0;
    if (factor_log2 + plain_scale_log2 >= 0.0) {
        application.result = multiplyPlainAndRescale(adapter, 
            input,
            adapter.encodeScalarAtScaleFor(std::exp2(factor_log2), std::exp2(plain_scale_log2), input));
        application.chunks = 1;
        application.levels_consumed = start_info.chainIndex - adapter.info(application.result).chainIndex;
        return application;
    }

    if (factor_log2 >= 0.0) {
        throw std::runtime_error("positive bootstrap scalar decomposition is not supported");
    }

    auto current = input;
    double remaining_abs_log2 = -factor_log2;
    constexpr double decomposition_epsilon_log2 = 1e-6;
    while (remaining_abs_log2 > decomposition_epsilon_log2) {
        const auto info = adapter.info(current);
        if (info.chainIndex == 0) {
            throw std::runtime_error("not enough levels for bootstrap scalar decomposition");
        }
        const double current_scale_log2 = std::log2(info.scale);
        const double capacity_log2 =
            info.coeffModulusLog2 - current_scale_log2 - scale_capacity_margin_log2;
        const double chunk_plain_scale_log2 = std::min(plain_scale_log2, capacity_log2);
        if (chunk_plain_scale_log2 <= 0.0) {
            throw std::runtime_error("not enough scale capacity for bootstrap scalar decomposition");
        }
        const double chunk_abs_log2 = std::min(remaining_abs_log2, chunk_plain_scale_log2);
        const double chunk_log2 = -chunk_abs_log2;
        current = multiplyPlainAndRescale(adapter, 
            current,
            adapter.encodeScalarAtScaleFor(
                std::exp2(chunk_log2), std::exp2(chunk_plain_scale_log2), current));
        remaining_abs_log2 = std::max(0.0, remaining_abs_log2 - chunk_abs_log2);
        ++application.chunks;
    }

    application.result = current;
    application.levels_consumed = start_info.chainIndex - adapter.info(application.result).chainIndex;
    return application;
}

BootstrapScaleSquash squash_bootstrap_scale(SealAdapter& adapter,
                                            const Cipher& input,
                                            double max_scale_log2,
                                            std::size_t min_chain_remaining) {
    if (!std::isfinite(max_scale_log2) || max_scale_log2 <= 0.0) {
        throw std::invalid_argument("max scale log2 must be positive and finite");
    }

    auto current = input;
    const auto start_info = adapter.info(current);
    BootstrapScaleSquash squash;
    while (std::log2(adapter.info(current).scale) > max_scale_log2) {
        const auto info = adapter.info(current);
        if (info.chainIndex <= min_chain_remaining) {
            throw std::runtime_error("not enough levels for bootstrap scale squash");
        }
        current = adapter.rescaleToNext(current);
        ++squash.levels_consumed;
    }
    squash.result = current;
    squash.levels_consumed = start_info.chainIndex - adapter.info(current).chainIndex;
    return squash;
}

BootstrapScaleStrategyPlan plan_bootstrap_scale_strategy(const std::vector<int>& active_coeff_modulus_bits,
                                                         const CipherInfo& start_info,
                                                         double factor_log2,
                                                         double max_plain_scale_log2,
                                                         double target_scale_log2,
                                                         std::size_t min_chain_remaining) {
    if (!std::isfinite(factor_log2)) {
        throw std::invalid_argument("bootstrap scalar factor log2 must be finite");
    }
    if (!std::isfinite(max_plain_scale_log2) || max_plain_scale_log2 <= 0.0) {
        throw std::invalid_argument("max plain scale log2 must be positive and finite");
    }
    if (!std::isfinite(target_scale_log2) || target_scale_log2 <= 0.0) {
        throw std::invalid_argument("target scale log2 must be positive and finite");
    }
    if (!std::isfinite(start_info.scale) || start_info.scale <= 0.0) {
        throw std::invalid_argument("start ciphertext scale must be positive and finite");
    }

    BootstrapScaleStrategyPlan plan;
    plan.factor_abs_log2 = factor_log2 < 0.0 ? -factor_log2 : factor_log2;
    plan.start_scale_log2 = std::log2(start_info.scale);
    plan.target_scale_log2 = target_scale_log2;
    plan.max_plain_scale_log2 = max_plain_scale_log2;
    plan.start_chain_index = start_info.chainIndex;
    plan.min_chain_remaining = min_chain_remaining;
    plan.scalar_representable = plan.factor_abs_log2 <=
        max_plain_scale_log2 * static_cast<double>(start_info.chainIndex);

    if (start_info.chainIndex < min_chain_remaining) {
        plan.blocker = "insufficient_start_chain";
        return plan;
    }
    plan.max_consumable_levels = start_info.chainIndex - min_chain_remaining;

    if (plan.factor_abs_log2 == 0.0) {
        plan.scalar_levels_needed = 0;
    } else {
        plan.scalar_levels_needed = static_cast<std::size_t>(
            std::ceil(plan.factor_abs_log2 / max_plain_scale_log2));
    }

    plan.required_drop_log2 = std::max(0.0, plan.start_scale_log2 + plan.factor_abs_log2 - target_scale_log2);
    const std::size_t available_drop_levels =
        std::min(plan.max_consumable_levels, active_coeff_modulus_bits.size());
    for (std::size_t i = 0; i < available_drop_levels; ++i) {
        plan.available_drop_log2 += static_cast<double>(active_coeff_modulus_bits[active_coeff_modulus_bits.size() - 1 - i]);
    }

    double cumulative_drop = 0.0;
    std::size_t levels_for_target = 0;
    while (levels_for_target < available_drop_levels && cumulative_drop < plan.required_drop_log2) {
        cumulative_drop += static_cast<double>(
            active_coeff_modulus_bits[active_coeff_modulus_bits.size() - 1 - levels_for_target]);
        ++levels_for_target;
    }
    plan.total_levels_needed = std::max(plan.scalar_levels_needed, levels_for_target);
    if (plan.total_levels_needed >= plan.scalar_levels_needed) {
        plan.scale_squash_levels_needed = plan.total_levels_needed - plan.scalar_levels_needed;
    }

    double scalar_drop = 0.0;
    const std::size_t scalar_drop_levels = std::min(plan.scalar_levels_needed, active_coeff_modulus_bits.size());
    for (std::size_t i = 0; i < scalar_drop_levels; ++i) {
        scalar_drop += static_cast<double>(active_coeff_modulus_bits[active_coeff_modulus_bits.size() - 1 - i]);
    }
    plan.scale_after_scalar_log2 = plan.start_scale_log2 + plan.factor_abs_log2 - scalar_drop;

    double total_drop = 0.0;
    const std::size_t total_drop_levels = std::min(plan.total_levels_needed, active_coeff_modulus_bits.size());
    for (std::size_t i = 0; i < total_drop_levels; ++i) {
        total_drop += static_cast<double>(active_coeff_modulus_bits[active_coeff_modulus_bits.size() - 1 - i]);
    }
    plan.scale_after_squash_log2 = plan.start_scale_log2 + plan.factor_abs_log2 - total_drop;

    if (!plan.scalar_representable) {
        plan.blocker = "scalar_not_representable";
    } else if (plan.scalar_levels_needed > plan.max_consumable_levels) {
        plan.missing_scalar_levels = plan.scalar_levels_needed - plan.max_consumable_levels;
        plan.blocker = "not_enough_levels_for_scalar";
    } else if (plan.required_drop_log2 > plan.available_drop_log2) {
        plan.missing_drop_log2 = plan.required_drop_log2 - plan.available_drop_log2;
        plan.blocker = "not_enough_levels_for_scale";
    } else if (plan.total_levels_needed > plan.max_consumable_levels) {
        plan.missing_total_levels = plan.total_levels_needed - plan.max_consumable_levels;
        plan.blocker = "not_enough_levels_for_scalar_and_scale";
    } else {
        plan.feasible = true;
        plan.blocker = "none";
    }
    return plan;
}

BootstrapEvalModCapacityPlan plan_evalmod_first_multiply_capacity(
    const std::vector<int>& active_coeff_modulus_bits,
    const BootstrapScaleStrategyPlan& scale_plan,
    double margin_log2) {
    if (!std::isfinite(margin_log2) || margin_log2 < 0.0) {
        throw std::invalid_argument("EvalMod capacity margin must be finite and non-negative");
    }

    BootstrapEvalModCapacityPlan plan;
    plan.margin_log2 = margin_log2;
    plan.first_product_scale_log2 = 2.0 * scale_plan.scale_after_squash_log2;

    if (!scale_plan.feasible) {
        plan.blocker = "scale_strategy_not_feasible";
        return plan;
    }
    if (scale_plan.total_levels_needed > active_coeff_modulus_bits.size()) {
        plan.blocker = "not_enough_active_moduli";
        return plan;
    }

    const std::size_t remaining_moduli = active_coeff_modulus_bits.size() - scale_plan.total_levels_needed;
    for (std::size_t i = 0; i < remaining_moduli; ++i) {
        plan.remaining_coeff_modulus_log2 += static_cast<double>(active_coeff_modulus_bits[i]);
    }
    plan.first_multiply_ready =
        plan.first_product_scale_log2 + margin_log2 <= plan.remaining_coeff_modulus_log2;
    plan.blocker = plan.first_multiply_ready ? "none" : "first_evalmod_multiply_scale";
    return plan;
}

} // namespace m2424
