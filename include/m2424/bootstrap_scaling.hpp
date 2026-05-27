#pragma once

#include "m2424/seal_adapter.hpp"

#include <cstddef>
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

} // namespace m2424
