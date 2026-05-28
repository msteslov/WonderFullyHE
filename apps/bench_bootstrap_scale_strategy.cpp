#include "m2424/bootstrap_scaling.hpp"
#include "m2424/diagonal_transform.hpp"
#include "m2424/eval_mod.hpp"
#include "m2424/profiles.hpp"
#include "m2424/seal_adapter.hpp"

#include <algorithm>
#include <cstddef>
#include <cmath>
#include <cstdio>
#include <vector>

namespace {

struct ProfileCase {
    const char* name{};
    m2424::CkksProfile profile{};
    double min_period_log2{};
    double max_period_log2{};
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

std::vector<int> active_coeff_modulus_bits(const m2424::CkksProfile& profile,
                                           const m2424::CipherInfo& info) {
    const std::size_t active_size = std::min(info.coeff_modulus_size, profile.coeff_modulus_bits.size());
    return {profile.coeff_modulus_bits.begin(),
            profile.coeff_modulus_bits.begin() + static_cast<std::ptrdiff_t>(active_size)};
}

void run_profile(const ProfileCase& profile_case) {
    constexpr std::size_t slots = 16;
    constexpr double amplitude = 1e-5;
    constexpr std::size_t level_drop = 2;
    constexpr double target_scale_log2 = 60.0;
    constexpr std::size_t min_chain_remaining = 3;
    const auto& profile = profile_case.profile;

    auto coeff_to_slot = m2424::DiagonalLinearTransform::from_matrix(
        m2424::canonical_embedding_matrix(slots));
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
    const auto active_bits = active_coeff_modulus_bits(profile, start_info);

    for (double period_log2 = profile_case.min_period_log2;
         period_log2 <= profile_case.max_period_log2;
         period_log2 += 2.0) {
        const double factor_log2 = -period_log2;
        const double expected_max_abs_after_normalization =
            max_abs_before_normalization * std::exp2(factor_log2);
        const bool evalmod_magnitude_ready =
            expected_max_abs_after_normalization <= m2424::EvalModPolynomial::approximation_bound;
        for (double max_plain_scale_log2 : {40.0, 50.0, 60.0, 80.0, 100.0, 120.0, 160.0, 200.0, 240.0}) {
            const auto plan = m2424::plan_bootstrap_scale_strategy(
                active_bits,
                start_info,
                factor_log2,
                max_plain_scale_log2,
                target_scale_log2,
                min_chain_remaining);
            const auto evalmod_capacity = m2424::plan_evalmod_first_multiply_capacity(
                active_bits, plan, 2.0);
            const bool p3_ready = evalmod_magnitude_ready && plan.feasible &&
                evalmod_capacity.first_multiply_ready;
            std::printf("%s,%.0f,%.0f,%zu,%zu,%.6e,%.6e,%.6e,%.6e,%zu,%zu,%zu,%zu,%zu,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%s,%s,%s,%s,%s\n",
                        profile_case.name,
                        period_log2,
                        max_plain_scale_log2,
                        start_info.chain_index,
                        start_info.coeff_modulus_size,
                        std::log2(start_info.scale),
                        start_info.coeff_modulus_log2,
                        max_abs_before_normalization,
                        expected_max_abs_after_normalization,
                        plan.max_consumable_levels,
                        plan.scalar_levels_needed,
                        plan.total_levels_needed,
                        plan.scale_squash_levels_needed,
                        plan.start_chain_index >= plan.total_levels_needed
                            ? plan.start_chain_index - plan.total_levels_needed
                            : 0,
                        plan.required_drop_log2,
                        plan.available_drop_log2,
                        plan.scale_after_scalar_log2,
                        plan.scale_after_squash_log2,
                        plan.target_scale_log2,
                        evalmod_capacity.remaining_coeff_modulus_log2,
                        evalmod_capacity.first_product_scale_log2,
                        evalmod_magnitude_ready ? "true" : "false",
                        plan.feasible ? "true" : "false",
                        evalmod_capacity.first_multiply_ready ? "true" : "false",
                        p3_ready ? "true" : "false",
                        plan.feasible ? evalmod_capacity.blocker.c_str() : plan.blocker.c_str());
        }
    }
}

} // namespace

int main() {
    std::printf("profile,period_log2,max_plain_scale_log2,start_chain_index,start_coeff_modulus_size,start_scale_log2,start_coeff_modulus_log2,max_abs_before_normalization,expected_max_abs_after_normalization,max_consumable_levels,scalar_levels_needed,total_levels_needed,scale_squash_levels_needed,chain_remaining_after_strategy,required_drop_log2,available_drop_log2,scale_after_scalar_log2,scale_after_squash_log2,target_scale_log2,remaining_coeff_modulus_log2_after_strategy,first_evalmod_product_scale_log2,evalmod_magnitude_ready,scale_strategy_feasible,evalmod_first_multiply_ready,p3_ready,blocker\n");
    run_profile({"boot_ckks", m2424::profiles::boot_ckks(), 220.0, 280.0});
    run_profile({"boot_deep_ckks", m2424::profiles::boot_deep_ckks(), 680.0, 800.0});
    return 0;
}
