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

void run_case(const char* profile_name,
              const m2424::CkksProfile& profile,
              double prescale_log2,
              double plain_scale_log2) {
    constexpr std::size_t slots = 16;
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
    current = coeff_to_slot.apply(adapter, current);

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
        plain_scale_log2,
        target_scale_log2,
        m2424::EvalModDegree::P3,
        active_bits,
        start_info,
        max_abs_before_normalization,
        min_chain_remaining,
        evalmod_capacity_margin_log2);

    const double expected_max_abs_after_normalization =
        max_abs_before_normalization * std::exp2(-candidate_period_log2);
    std::printf("%s,%.0f,%.0f,%zu,%.6e,%.6e,%.6e,%.6e,%.6e,%s,%.0f,%.6e,%zu,%zu,%zu,%zu,%.6e,%.6e,%.6e,%.6e,%s,%s,,PASS\n",
                profile_name,
                prescale_log2,
                plain_scale_log2,
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
                design.blocker.c_str());
}

} // namespace

int main() {
    std::printf("profile,prescale_log2,plain_scale_log2,chain_index,start_scale_log2,start_coeff_modulus_log2,max_abs_before_normalization,min_period_for_magnitude_log2,max_period_for_evalmod_capacity_log2,period_window_possible,candidate_period_log2,expected_max_abs_after_normalization,required_levels,chain_remaining_after_strategy,scalar_levels_needed,scale_squash_levels_needed,required_drop_log2,available_drop_log2,scale_after_squash_log2,remaining_coeff_modulus_log2_after_strategy,scale_design_status,blocker,exception,status\n");
    for (double prescale_log2 : {0.0, 40.0, 60.0, 80.0, 100.0, 120.0}) {
        try {
            run_case("boot_ckks", m2424::profiles::boot_ckks(), prescale_log2, 160.0);
        } catch (const std::exception& e) {
            auto message = std::string(e.what());
            std::replace(message.begin(), message.end(), ',', ';');
            std::printf("boot_ckks,%.0f,160,0,0,0,0,0,0,false,0,0,0,0,0,0,0,0,0,0,exception,%s,%s,FAIL\n",
                        prescale_log2,
                        message.c_str(),
                        message.c_str());
        }
    }
    return 0;
}
