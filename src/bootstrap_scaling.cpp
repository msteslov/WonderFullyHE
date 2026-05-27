#include "m2424/bootstrap_scaling.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace m2424 {

const char* to_string(BootstrapPeriodMode mode) noexcept {
    switch (mode) {
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

double bootstrap_period_log2(BootstrapPeriodMode mode,
                             double manual_period_log2,
                             const std::vector<int>& coeff_modulus_bits,
                             const CipherInfo& before_mod_raise,
                             const CipherInfo& after_mod_raise) {
    switch (mode) {
    case BootstrapPeriodMode::TotalCoeffModulus:
        return after_mod_raise.coeff_modulus_log2 - std::log2(after_mod_raise.scale);
    case BootstrapPeriodMode::LastPrime: {
        if (coeff_modulus_bits.empty()) {
            throw std::invalid_argument("coeff_modulus_bits must not be empty");
        }
        if (before_mod_raise.coeff_modulus_size >= coeff_modulus_bits.size()) {
            return static_cast<double>(coeff_modulus_bits.back());
        }
        return static_cast<double>(coeff_modulus_bits[before_mod_raise.coeff_modulus_size]);
    }
    case BootstrapPeriodMode::DroppedPrimeProduct:
        return std::max(0.0, after_mod_raise.coeff_modulus_log2 - before_mod_raise.coeff_modulus_log2);
    case BootstrapPeriodMode::ManualPowerOfTwo:
        return manual_period_log2;
    }
    return 0.0;
}

BootstrapScalingFactors make_bootstrap_scaling_factors(double amplitude_factor,
                                                       double bootstrap_period_log2,
                                                       double plain_scale_log2) {
    if (!std::isfinite(amplitude_factor) || amplitude_factor <= 0.0) {
        throw std::invalid_argument("amplitude factor must be positive and finite");
    }
    if (!std::isfinite(bootstrap_period_log2)) {
        throw std::invalid_argument("bootstrap period log2 must be finite");
    }
    if (!std::isfinite(plain_scale_log2)) {
        throw std::invalid_argument("plain scale log2 must be finite");
    }

    BootstrapScalingFactors factors;
    factors.bootstrap_period_log2 = bootstrap_period_log2;
    factors.normalization_factor_log2 = std::log2(amplitude_factor) - bootstrap_period_log2;
    factors.plain_scale_log2 = plain_scale_log2;
    factors.factor_times_plain_scale_log2 =
        factors.normalization_factor_log2 + factors.plain_scale_log2;
    factors.representable = factors.factor_times_plain_scale_log2 >= 0.0;
    factors.factor = std::exp2(factors.normalization_factor_log2);
    return factors;
}

} // namespace m2424
