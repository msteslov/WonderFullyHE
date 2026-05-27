#pragma once

#include "m2424/seal_adapter.hpp"

#include <vector>

namespace m2424 {

enum class BootstrapPeriodMode {
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

const char* to_string(BootstrapPeriodMode mode) noexcept;

double bootstrap_period_log2(BootstrapPeriodMode mode,
                             double manual_period_log2,
                             const std::vector<int>& coeff_modulus_bits,
                             const CipherInfo& before_mod_raise,
                             const CipherInfo& after_mod_raise);

BootstrapScalingFactors make_bootstrap_scaling_factors(double amplitude_factor,
                                                       double bootstrap_period_log2,
                                                       double plain_scale_log2);

} // namespace m2424
