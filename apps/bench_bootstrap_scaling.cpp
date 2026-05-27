#include "m2424/diagonal_transform.hpp"
#include "m2424/bootstrap_scaling.hpp"
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

struct PeriodCase {
    m2424::BootstrapPeriodMode mode{};
    double manual_period_log2{};
};

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

m2424::Cipher apply_scalar_decomposed(m2424::SealAdapter& adapter,
                                      const m2424::Cipher& input,
                                      double factor_log2,
                                      double plain_scale_log2) {
    constexpr double scale_capacity_margin_log2 = 2.0;
    if (factor_log2 + plain_scale_log2 >= 0.0) {
        return adapter.mul_plain_rescale(
            input,
            adapter.encode_scalar_at_scale_like(std::exp2(factor_log2), std::exp2(plain_scale_log2), input));
    }

    if (factor_log2 >= 0.0) {
        throw std::runtime_error("positive bootstrap scalar decomposition is not supported");
    }

    auto current = input;
    double remaining_abs_log2 = -factor_log2;
    while (remaining_abs_log2 > 1e-9) {
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
        remaining_abs_log2 -= chunk_abs_log2;
    }
    return current;
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

} // namespace

int main() {
    constexpr std::size_t slots = 16;
    constexpr double tolerance = 2e-5;
    const auto profile = m2424::profiles::boot_ckks();

    const std::vector<double> amplitudes{1e-5, 1e-4, 5e-4, 1e-3, 1e-2, 1e-1};
    const std::vector<std::size_t> level_drops{1, 2};
    const std::vector<double> plain_scale_log2_values{
        40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 120.0,
        140.0, 160.0, 180.0, 200.0, 240.0, 280.0, 320.0, 400.0, 600.0
    };
    std::vector<PeriodCase> period_cases{
        {m2424::BootstrapPeriodMode::NoBootstrapPeriod, 0.0},
        {m2424::BootstrapPeriodMode::SourceCoeffModulus, 0.0},
        {m2424::BootstrapPeriodMode::TotalCoeffModulus, 0.0},
        {m2424::BootstrapPeriodMode::LastPrime, 0.0},
        {m2424::BootstrapPeriodMode::DroppedPrimeProduct, 0.0}
    };
    for (double manual : {40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 120.0, 140.0, 220.0, 258.0, 260.0, 300.0}) {
        period_cases.push_back({m2424::BootstrapPeriodMode::ManualPowerOfTwo, manual});
    }

    auto rotation_steps = coeff_to_slot_rotation_steps(slots);
    auto adapter = m2424::SealAdapter::create(profile);
    adapter.keygen(rotation_steps, true);
    auto coeff_to_slot = m2424::DiagonalLinearTransform::from_matrix(
        m2424::canonical_embedding_matrix(slots));

    std::printf("case,period_mode,manual_period_log2,plain_scale_log2,amplitude,level_drop,chain_before,chain_after,cipher_scale_log2,bootstrap_period_log2,amplitude_factor_log2,period_contribution_log2,normalization_factor_log2,factor_times_plain_scale_log2,required_plain_scale_log2,plain_scale_margin_log2,representable,coeff_to_slot_max_abs_before_normalization,expected_max_abs_after_normalization,actual_max_abs_after_normalization,normalization_error,exception,status\n");

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
            const double amplitude_factor_log2 = std::log2(amplitude_factor);
            const double coeff_to_slot_max_abs = max_abs_value(before_normalization);

            for (const auto& period_case : period_cases) {
                for (double plain_scale_log2 : plain_scale_log2_values) {
                    ++case_id;
                    const auto before_normalization_info = adapter.info(current);
                    const double cipher_scale_log2 = std::log2(before_normalization_info.scale);
                    const double period_log2 = m2424::bootstrap_period_log2(
                        period_case.mode,
                        period_case.manual_period_log2,
                        profile.coeff_modulus_bits,
                        before_mod_raise,
                        after_mod_raise);
                    const auto scaling = m2424::make_bootstrap_scaling_factors(
                        amplitude_factor, period_log2, plain_scale_log2);
                    const double normalization_factor_log2 = scaling.normalization_factor_log2;
                    const double normalization_factor = scaling.factor;
                    const double factor_times_plain_scale_log2 = scaling.factor_times_plain_scale_log2;
                    const double required_plain_scale_log2 = -normalization_factor_log2;
                    const double plain_scale_margin_log2 = factor_times_plain_scale_log2;

                    double expected_max = 0.0;
                    double actual_max = 0.0;
                    double normalization_error = 0.0;
                    std::string exception;
                    bool scalar_representable = false;
                    const char* status = "FAIL";

                    try {
                        const auto expected = scaled(before_normalization, normalization_factor);
                        expected_max = max_abs_value(expected);
                        auto normalized = apply_scalar_decomposed(
                            adapter, current, normalization_factor_log2, plain_scale_log2);
                        scalar_representable = true;
                        auto actual = head(adapter.decode_complex(adapter.decrypt(normalized)), slots);
                        actual_max = max_abs_value(actual);
                        normalization_error = max_complex_error(expected, actual);
                        status = normalization_error <= tolerance ? "PASS" : "FAIL";
                    } catch (const std::exception& e) {
                        exception = e.what();
                        std::replace(exception.begin(), exception.end(), ',', ';');
                        status = std::string(exception) == "bootstrap scalar chunk is not representable"
                            ? "SCALAR_NOT_REPRESENTABLE"
                            : "FAIL";
                    }

                    std::printf("%zu,%s,%.0f,%.0f,%.6e,%zu,%zu,%zu,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%s,%.6e,%.6e,%.6e,%.6e,%s,%s\n",
                                case_id,
                                m2424::to_string(period_case.mode),
                                period_case.mode == m2424::BootstrapPeriodMode::ManualPowerOfTwo
                                    ? period_case.manual_period_log2
                                    : 0.0,
                                plain_scale_log2,
                                amplitude,
                                level_drop,
                                before_mod_raise.chain_index,
                                before_normalization_info.chain_index,
                                cipher_scale_log2,
                                period_log2,
                                amplitude_factor_log2,
                                period_log2,
                                normalization_factor_log2,
                                factor_times_plain_scale_log2,
                                required_plain_scale_log2,
                                plain_scale_margin_log2,
                                scalar_representable ? "true" : "false",
                                coeff_to_slot_max_abs,
                                expected_max,
                                actual_max,
                                normalization_error,
                                exception.c_str(),
                                status);
                }
            }
        }
    }

    return 0;
}
