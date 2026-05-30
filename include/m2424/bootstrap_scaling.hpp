#pragma once

#include "m2424/seal_adapter.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace m2424 {

enum class BootstrapPeriodMode {
    NoBootstrapPeriod,
    SourceCoeffModulus,
    TotalCoeffModulus,
    LastPrime,
    DroppedPrimeProduct,
    ManualPowerOfTwo
};

struct BootstrapScalingFactors {
    double bootstrap_period_log2{};
    double normalization_factor_log2{};
    double plain_scale_log2{};
    double factor_times_plain_scale_log2{};
    bool representable{};
    double factor{};
};

struct BootstrapScalarApplication {
    Cipher result;
    std::size_t chunks{};
    std::size_t levels_consumed{};
};

struct BootstrapScaleSquash {
    Cipher result;
    std::size_t levels_consumed{};
};

struct BootstrapScaleStrategyPlan {
    double factor_abs_log2{};
    double start_scale_log2{};
    double target_scale_log2{};
    double max_plain_scale_log2{};
    std::size_t start_chain_index{};
    std::size_t min_chain_remaining{};
    std::size_t max_consumable_levels{};
    std::size_t scalar_levels_needed{};
    std::size_t total_levels_needed{};
    std::size_t scale_squash_levels_needed{};
    double required_drop_log2{};
    double available_drop_log2{};
    double missing_drop_log2{};
    double scale_after_scalar_log2{};
    double scale_after_squash_log2{};
    std::size_t missing_scalar_levels{};
    std::size_t missing_total_levels{};
    bool scalar_representable{};
    bool feasible{};
    std::string blocker;
};

struct BootstrapEvalModCapacityPlan {
    double remaining_coeff_modulus_log2{};
    double first_product_scale_log2{};
    double margin_log2{};
    bool first_multiply_ready{};
    std::string blocker;
};

const char* to_string(BootstrapPeriodMode mode) noexcept;

double bootstrap_period_log2(BootstrapPeriodMode mode,
                             double manual_period_log2,
                             const std::vector<int>& coeff_modulus_bits,
                             const CipherInfo& before_mod_raise,
                             const CipherInfo& after_mod_raise);

BootstrapScalingFactors make_bootstrap_scaling_factors(double amplitude_factor,
                                                       double bootstrap_period_log2,
                                                       double plain_scale_log2);

BootstrapScalarApplication apply_bootstrap_scalar_decomposed(SealAdapter& adapter,
                                                             const Cipher& input,
                                                             double factor_log2,
                                                             double plain_scale_log2);

BootstrapScaleSquash squash_bootstrap_scale(SealAdapter& adapter,
                                            const Cipher& input,
                                            double max_scale_log2,
                                            std::size_t min_chain_remaining);

BootstrapScaleStrategyPlan plan_bootstrap_scale_strategy(const std::vector<int>& active_coeff_modulus_bits,
                                                         const CipherInfo& start_info,
                                                         double factor_log2,
                                                         double max_plain_scale_log2,
                                                         double target_scale_log2,
                                                         std::size_t min_chain_remaining);

BootstrapEvalModCapacityPlan plan_evalmod_first_multiply_capacity(
    const std::vector<int>& active_coeff_modulus_bits,
    const BootstrapScaleStrategyPlan& scale_plan,
    double margin_log2);

} // namespace m2424
