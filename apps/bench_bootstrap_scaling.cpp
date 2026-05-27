#include "m2424/diagonal_transform.hpp"
#include "m2424/profiles.hpp"
#include "m2424/seal_adapter.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>
#include <limits>
#include <string>
#include <vector>

namespace {

enum class BootstrapPeriodMode {
    TotalCoeffModulus,
    LastPrime,
    DroppedPrimeProduct,
    ManualPowerOfTwo
};

struct PeriodCase {
    BootstrapPeriodMode mode{};
    double manual_period_log2{};
};

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

std::vector<double> make_input(std::size_t slots, double amplitude) {
    std::vector<double> input;
    input.reserve(slots);
    for (std::size_t i = 0; i < slots; ++i) {
        input.push_back(amplitude * std::sin(static_cast<double>(i) / 4.0));
    }
    return input;
}

m2424::ComplexVector head(m2424::ComplexVector values, std::size_t n) {
    values.resize(std::min(values.size(), n));
    return values;
}

double max_abs_value(const m2424::ComplexVector& values) {
    double result = 0.0;
    for (const auto& value : values) {
        result = std::max(result, std::abs(value));
    }
    return result;
}

double max_complex_error(const m2424::ComplexVector& expected, const m2424::ComplexVector& actual) {
    const std::size_t n = std::min(expected.size(), actual.size());
    double result = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        result = std::max(result, std::abs(expected[i] - actual[i]));
    }
    return result;
}

m2424::ComplexVector scaled(const m2424::ComplexVector& values, double factor) {
    m2424::ComplexVector result;
    result.reserve(values.size());
    for (const auto& value : values) {
        result.push_back(value * factor);
    }
    return result;
}

std::vector<int> coeff_to_slot_rotation_steps(std::size_t slots) {
    auto coeff_to_slot = m2424::DiagonalLinearTransform::from_matrix(
        m2424::canonical_embedding_matrix(slots));
    return coeff_to_slot.rotation_steps();
}

double normalization_factor_for(const m2424::ComplexVector& values) {
    constexpr double evalmod_target = 0.0009765625 * 0.5;
    const double max_abs = max_abs_value(values);
    if (max_abs == 0.0) {
        return 1.0;
    }
    return max_abs > evalmod_target ? evalmod_target / max_abs : 1.0;
}

double bootstrap_period_log2(const PeriodCase& period_case,
                             const m2424::CkksProfile& profile,
                             const m2424::CipherInfo& before_mod_raise,
                             const m2424::CipherInfo& after_mod_raise) {
    switch (period_case.mode) {
    case BootstrapPeriodMode::TotalCoeffModulus:
        return after_mod_raise.coeff_modulus_log2 - std::log2(after_mod_raise.scale);
    case BootstrapPeriodMode::LastPrime: {
        if (before_mod_raise.coeff_modulus_size >= profile.coeff_modulus_bits.size()) {
            return static_cast<double>(profile.coeff_modulus_bits.back());
        }
        return static_cast<double>(profile.coeff_modulus_bits[before_mod_raise.coeff_modulus_size]);
    }
    case BootstrapPeriodMode::DroppedPrimeProduct:
        return std::max(0.0, after_mod_raise.coeff_modulus_log2 - before_mod_raise.coeff_modulus_log2);
    case BootstrapPeriodMode::ManualPowerOfTwo:
        return period_case.manual_period_log2;
    }
    return 0.0;
}

} // namespace

int main() {
    constexpr std::size_t slots = 16;
    constexpr double tolerance = 2e-5;
    const auto profile = m2424::profiles::boot_ckks();
    const double plain_scale_log2 = std::log2(profile.scale);

    const std::vector<double> amplitudes{1e-5, 1e-4, 5e-4, 1e-3, 1e-2, 1e-1};
    const std::vector<std::size_t> level_drops{1, 2};
    std::vector<PeriodCase> period_cases{
        {BootstrapPeriodMode::TotalCoeffModulus, 0.0},
        {BootstrapPeriodMode::LastPrime, 0.0},
        {BootstrapPeriodMode::DroppedPrimeProduct, 0.0}
    };
    for (double manual : {40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 120.0, 140.0}) {
        period_cases.push_back({BootstrapPeriodMode::ManualPowerOfTwo, manual});
    }

    auto rotation_steps = coeff_to_slot_rotation_steps(slots);
    auto adapter = m2424::SealAdapter::create(profile);
    adapter.keygen(rotation_steps, true);
    auto coeff_to_slot = m2424::DiagonalLinearTransform::from_matrix(
        m2424::canonical_embedding_matrix(slots));

    std::printf("case,period_mode,manual_period_log2,amplitude,chain_before,chain_after,cipher_scale_log2,bootstrap_period_log2,normalization_factor_log2,plain_scale_log2,factor_times_plain_scale_log2,representable,expected_max_abs_after_normalization,actual_max_abs_after_normalization,normalization_error,status\n");

    std::size_t case_id = 0;
    for (double amplitude : amplitudes) {
        for (std::size_t level_drop : level_drops) {
            const auto input = make_input(slots, amplitude);
            auto current = adapter.encrypt(adapter.encode(input));
            for (std::size_t i = 0; i < level_drop; ++i) {
                current = adapter.mul_plain_rescale(current, adapter.encode_scalar_like(1.0, current));
            }
            const auto before_mod_raise = adapter.info(current);
            current = adapter.mod_raise_to_first(current);
            const auto after_mod_raise = adapter.info(current);
            current = coeff_to_slot.apply(adapter, current);

            m2424::ComplexVector before_normalization;
            try {
                before_normalization = head(adapter.decode_complex(adapter.decrypt(current)), slots);
            } catch (...) {
                before_normalization.assign(slots, {0.0, 0.0});
            }
            const double amplitude_factor = normalization_factor_for(before_normalization);

            for (const auto& period_case : period_cases) {
                ++case_id;
                const auto before_normalization_info = adapter.info(current);
                const double cipher_scale_log2 = std::log2(before_normalization_info.scale);
                const double period_log2 = bootstrap_period_log2(
                    period_case, profile, before_mod_raise, after_mod_raise);
                const double normalization_factor_log2 = std::log2(amplitude_factor) - period_log2;
                const double factor_times_plain_scale_log2 =
                    normalization_factor_log2 + plain_scale_log2;
                const bool representable = factor_times_plain_scale_log2 >= 0.0;

                double expected_max = 0.0;
                double actual_max = 0.0;
                double normalization_error = 0.0;
                const char* status = "FAIL";

                if (!representable) {
                    status = "SCALAR_NOT_REPRESENTABLE";
                } else {
                    try {
                        const double normalization_factor = std::exp2(normalization_factor_log2);
                        const auto expected = scaled(before_normalization, normalization_factor);
                        expected_max = max_abs_value(expected);
                        const auto plain = adapter.encode_scalar_at_scale_like(
                            normalization_factor, profile.scale, current);
                        auto normalized = adapter.mul_plain_rescale(current, plain);
                        auto actual = head(adapter.decode_complex(adapter.decrypt(normalized)), slots);
                        actual_max = max_abs_value(actual);
                        normalization_error = max_complex_error(expected, actual);
                        status = normalization_error <= tolerance ? "PASS" : "FAIL";
                    } catch (const std::exception&) {
                        status = "FAIL";
                    }
                }

                std::printf("%zu,%s,%.0f,%.6e,%zu,%zu,%.6e,%.6e,%.6e,%.6e,%.6e,%s,%.6e,%.6e,%.6e,%s\n",
                            case_id,
                            to_string(period_case.mode),
                            period_case.mode == BootstrapPeriodMode::ManualPowerOfTwo
                                ? period_case.manual_period_log2
                                : 0.0,
                            amplitude,
                            before_mod_raise.chain_index,
                            before_normalization_info.chain_index,
                            cipher_scale_log2,
                            period_log2,
                            normalization_factor_log2,
                            plain_scale_log2,
                            factor_times_plain_scale_log2,
                            representable ? "true" : "false",
                            expected_max,
                            actual_max,
                            normalization_error,
                            status);
            }
        }
    }

    return 0;
}
