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

double bootstrap_period_log2(BootstrapPeriodMode mode,
                             double manual_period_log2,
                             const std::vector<int>& coeff_modulus_bits,
                             const CipherInfo& before_mod_raise,
                             const CipherInfo& after_mod_raise) {
    switch (mode) {
    case BootstrapPeriodMode::NoBootstrapPeriod:
        return 0.0;
    case BootstrapPeriodMode::SourceCoeffModulus:
        if (!std::isfinite(before_mod_raise.scale) || before_mod_raise.scale <= 0.0) {
            throw std::invalid_argument("source ciphertext scale must be positive and finite");
        }
        return before_mod_raise.coeff_modulus_log2 - std::log2(before_mod_raise.scale);
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
    factors.normalization_factor_log2 = bootstrap_period_log2 > 0.0
        ? -bootstrap_period_log2
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
        application.result = adapter.mul_plain_rescale(
            input,
            adapter.encode_scalar_at_scale_like(std::exp2(factor_log2), std::exp2(plain_scale_log2), input));
        application.chunks = 1;
        application.levels_consumed = start_info.chain_index - adapter.info(application.result).chain_index;
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
        if (info.chain_index == 0) {
            throw std::runtime_error("not enough levels for bootstrap scalar decomposition");
        }
        const double current_scale_log2 = std::log2(info.scale);
        const double capacity_log2 =
            info.coeff_modulus_log2 - current_scale_log2 - scale_capacity_margin_log2;
        const double chunk_plain_scale_log2 = std::min(plain_scale_log2, capacity_log2);
        if (chunk_plain_scale_log2 <= 0.0) {
            throw std::runtime_error("not enough scale capacity for bootstrap scalar decomposition");
        }
        const double chunk_abs_log2 = std::min(remaining_abs_log2, chunk_plain_scale_log2);
        const double chunk_log2 = -chunk_abs_log2;
        current = adapter.mul_plain_rescale(
            current,
            adapter.encode_scalar_at_scale_like(
                std::exp2(chunk_log2), std::exp2(chunk_plain_scale_log2), current));
        remaining_abs_log2 = std::max(0.0, remaining_abs_log2 - chunk_abs_log2);
        ++application.chunks;
    }

    application.result = current;
    application.levels_consumed = start_info.chain_index - adapter.info(application.result).chain_index;
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
        if (info.chain_index <= min_chain_remaining) {
            throw std::runtime_error("not enough levels for bootstrap scale squash");
        }
        current = adapter.rescale_to_next(current);
        ++squash.levels_consumed;
    }
    squash.result = current;
    squash.levels_consumed = start_info.chain_index - adapter.info(current).chain_index;
    return squash;
}

} // namespace m2424
