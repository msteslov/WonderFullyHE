#include "m2424/bootstrap_plan.hpp"
#include "m2424/diagonal_transform.hpp"
#include "m2424/eval_mod.hpp"
#include "m2424/profiles.hpp"
#include "m2424/seal_adapter.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>
#include <string>
#include <vector>

namespace {

std::vector<double> make_input(std::size_t slots, double amplitude) {
    std::vector<double> input;
    input.reserve(slots);
    for (std::size_t i = 0; i < slots; ++i) {
        input.push_back(amplitude * std::sin(static_cast<double>(i) / 4.0));
    }
    return input;
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

m2424::ComplexMatrix scaled_matrix(m2424::ComplexMatrix matrix, double factor) {
    for (auto& row : matrix) {
        for (auto& value : row) {
            value *= factor;
        }
    }
    return matrix;
}

m2424::ComplexVector head(m2424::ComplexVector values, std::size_t n) {
    values.resize(std::min(values.size(), n));
    return values;
}

bool run_case(const char* profile_name,
              const m2424::CkksProfile& profile,
              double prescale_log2,
              double transform_plain_scale_log2,
              double normalization_plain_scale_log2) {
    constexpr std::size_t slots = 16;
    constexpr double tolerance = 2e-5;
    constexpr double amplitude = 1e-5;
    constexpr std::size_t level_drop = 2;
    constexpr double target_scale_log2 = 60.0;
    constexpr std::size_t min_chain_remaining = 3;
    constexpr double evalmod_capacity_margin_log2 = 2.0;

    auto coeff_to_slot = m2424::DiagonalLinearTransform::from_matrix(
        scaled_matrix(m2424::canonical_embedding_matrix(slots), std::exp2(-prescale_log2)));

    auto adapter = m2424::SealAdapter::create(profile);
    adapter.keygen(coeff_to_slot.rotation_steps(), true);

    auto current = adapter.encrypt(adapter.encode(make_input(slots, amplitude)));
    for (std::size_t i = 0; i < level_drop; ++i) {
        current = adapter.mul_plain_rescale(current, adapter.encode_scalar_like(1.0, current));
    }
    current = adapter.mod_raise_to_first(current);
    current = coeff_to_slot.apply_at_plain_scale(adapter, current, std::exp2(transform_plain_scale_log2));

    const auto start_info = adapter.info(current);
    const auto before_normalization = head(adapter.decode_complex(adapter.decrypt(current)), slots);
    const double max_abs_before_normalization = max_abs_value(before_normalization);
    const auto active_bits = m2424::active_coeff_modulus_bits(profile, start_info);
    const auto window = m2424::bootstrap_period_feasibility_window(
        start_info.coeff_modulus_log2,
        std::log2(start_info.scale),
        target_scale_log2,
        max_abs_before_normalization,
        m2424::EvalModPolynomial::approximation_bound,
        evalmod_capacity_margin_log2);

    const double candidate_period_log2 = std::ceil(window.min_period_for_magnitude_log2);
    const auto design = m2424::make_bootstrap_scale_design(
        m2424::BootstrapPeriodMode::ManualPowerOfTwo,
        candidate_period_log2,
        candidate_period_log2,
        m2424::BootstrapScalingStrategy::DecomposedPlainMultiplyRescale,
        normalization_plain_scale_log2,
        target_scale_log2,
        m2424::EvalModDegree::P3,
        active_bits,
        start_info,
        max_abs_before_normalization,
        min_chain_remaining,
        evalmod_capacity_margin_log2);

    const double expected_max_abs_after_normalization =
        max_abs_before_normalization * std::exp2(-candidate_period_log2);
    double normalization_error = 0.0;
    double scale_log2_after_normalization = 0.0;
    std::size_t chain_after_normalization = 0;
    double scale_log2_after_squash = 0.0;
    std::size_t chain_after_squash = 0;
    std::size_t scale_squash_levels_consumed = 0;
    bool scalar_pass = false;
    bool evalmod_interval_ready = false;
    bool scale_squash_ready = false;
    bool p3_ready = false;
    std::string physical_exception;

    try {
        const auto expected = scaled(before_normalization, std::exp2(-candidate_period_log2));
        auto normalized = m2424::apply_bootstrap_scalar_decomposed(
            adapter, current, -candidate_period_log2, normalization_plain_scale_log2);
        const auto normalized_info = adapter.info(normalized.result);
        scale_log2_after_normalization = std::log2(normalized_info.scale);
        chain_after_normalization = normalized_info.chain_index;
        const auto actual = head(adapter.decode_complex(adapter.decrypt(normalized.result)), slots);
        normalization_error = max_complex_error(expected, actual);
        scalar_pass = normalization_error <= tolerance;
        evalmod_interval_ready = expected_max_abs_after_normalization <= m2424::EvalModPolynomial::approximation_bound;
        auto squashed = m2424::squash_bootstrap_scale(
            adapter, normalized.result, target_scale_log2, min_chain_remaining);
        const auto squashed_info = adapter.info(squashed.result);
        scale_log2_after_squash = std::log2(squashed_info.scale);
        chain_after_squash = squashed_info.chain_index;
        scale_squash_levels_consumed = squashed.levels_consumed;
        scale_squash_ready = scale_log2_after_squash <= target_scale_log2
            && chain_after_squash >= min_chain_remaining;
        p3_ready = scalar_pass && evalmod_interval_ready && scale_squash_ready;
    } catch (const std::exception& e) {
        physical_exception = e.what();
        std::replace(physical_exception.begin(), physical_exception.end(), ',', ';');
    }

    std::printf("%s,%.0f,%.0f,%.0f,%zu,%.6e,%.6e,%.6e,%.6e,%.6e,%s,%.0f,%.6e,%zu,%zu,%zu,%zu,%.6e,%.6e,%.6e,%.6e,%s,%s,%.6e,%zu,%.6e,%zu,%zu,%s,%s,%s,%s,%s,%s\n",
                profile_name,
                prescale_log2,
                transform_plain_scale_log2,
                normalization_plain_scale_log2,
                start_info.chain_index,
                std::log2(start_info.scale),
                start_info.coeff_modulus_log2,
                max_abs_before_normalization,
                window.min_period_for_magnitude_log2,
                window.max_period_for_evalmod_capacity_log2,
                window.possible ? "true" : "false",
                candidate_period_log2,
                expected_max_abs_after_normalization,
                design.required_levels,
                design.chain_remaining_after_strategy,
                design.scale_plan.scalar_levels_needed,
                design.scale_plan.scale_squash_levels_needed,
                design.scale_plan.required_drop_log2,
                design.scale_plan.available_drop_log2,
                design.scale_plan.scale_after_squash_log2,
                design.evalmod_capacity.remaining_coeff_modulus_log2,
                m2424::to_string(design.status),
                design.blocker.c_str(),
                normalization_error,
                chain_after_normalization,
                scale_log2_after_normalization,
                chain_after_squash,
                scale_squash_levels_consumed,
                scalar_pass ? "true" : "false",
                evalmod_interval_ready ? "true" : "false",
                scale_squash_ready ? "true" : "false",
                p3_ready ? "true" : "false",
                physical_exception.c_str(),
                physical_exception.empty() ? "PASS" : "FAIL");
    return p3_ready;
}

} // namespace

int main() {
    std::printf("profile,prescale_log2,transform_plain_scale_log2,normalization_plain_scale_log2,chain_index,start_scale_log2,start_coeff_modulus_log2,max_abs_before_normalization,min_period_for_magnitude_log2,max_period_for_evalmod_capacity_log2,period_window_possible,candidate_period_log2,expected_max_abs_after_normalization,required_levels,chain_remaining_after_strategy,scalar_levels_needed,scale_squash_levels_needed,required_drop_log2,available_drop_log2,scale_after_squash_log2,remaining_coeff_modulus_log2_after_strategy,scale_design_status,blocker,normalization_error,chain_after_normalization,scale_log2_after_normalization,chain_after_squash,scale_squash_levels_consumed,scalar_pass,evalmod_interval_ready,scale_squash_ready,p3_ready,exception,status\n");
    std::size_t cases = 0;
    std::size_t p3_ready_cases = 0;
    for (double prescale_log2 : {0.0, 40.0, 80.0, 120.0}) {
        for (double transform_plain_scale_log2 : {40.0, 80.0, 120.0, 160.0}) {
            ++cases;
            try {
                if (run_case("boot_ckks",
                             m2424::profiles::boot_ckks(),
                             prescale_log2,
                             transform_plain_scale_log2,
                             160.0)) {
                    ++p3_ready_cases;
                }
            } catch (const std::exception& e) {
                auto message = std::string(e.what());
                std::replace(message.begin(), message.end(), ',', ';');
                std::printf("boot_ckks,%.0f,%.0f,160,0,0,0,0,0,0,false,0,0,0,0,0,0,0,0,0,0,exception,%s,0,0,0,0,0,false,false,false,false,%s,FAIL\n",
                            prescale_log2,
                            transform_plain_scale_log2,
                            message.c_str(),
                            message.c_str());
            }
        }
    }
    std::printf("summary,boot_ckks,%zu,%zu\n", cases, p3_ready_cases);
    return 0;
}
