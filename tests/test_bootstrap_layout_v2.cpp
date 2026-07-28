#include "m2424/bootstrap_layout_v2.hpp"

#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>

namespace {

std::size_t sum_bits(const std::vector<int>& bits) {
    return static_cast<std::size_t>(std::accumulate(bits.begin(), bits.end(), 0));
}

bool near(double lhs, double rhs, double tolerance) {
    return std::abs(lhs - rhs) <= tolerance;
}

} // namespace

int main() {
    bool ok = true;

    const auto layout = m2424::make_lattigo_like_bootstrap_layout_v2(15, 2, 1e-9, 3);

    ok = ok && layout.residual_profile.polyModulusDegree == 32768;
    ok = ok && layout.residual_profile.slots == 4;
    ok = ok && near(std::log2(layout.residual_profile.scale), 50.0, 1e-9);

    ok = ok && layout.slots_to_coeffs.name == "slots_to_coeffs";
    ok = ok && layout.slots_to_coeffs.log_scale >= 38.0;
    ok = ok && layout.slots_to_coeffs.log_scale <= 40.0;
    ok = ok && layout.slots_to_coeffs.levels == layout.slots_to_coeffs.log_q.size();

    ok = ok && layout.coeffs_to_slots.name == "coeffs_to_slots";
    ok = ok && layout.coeffs_to_slots.log_scale >= 55.0;
    ok = ok && layout.coeffs_to_slots.log_scale <= 57.0;

    ok = ok && layout.evalmod.name.find("degree30") != std::string::npos;
    ok = ok && layout.evalmod.name.find("double_angle3") != std::string::npos;
    ok = ok && near(layout.evalmod.log_scale, 60.0, 1e-9);
    ok = ok && layout.evalmod.levels >= 8;

    ok = ok && layout.key_switch.log_p == std::vector<int>({61, 61, 61});
    ok = ok && layout.key_switch.decomposition_count == 3;
    ok = ok && layout.key_switch.uses_special_primes;
    ok = ok && layout.key_switch.uses_ephemeral_sparse_secret;
    ok = ok && layout.key_switch.ephemeral_secret_weight == 32;

    const std::size_t expected_q = sum_bits(layout.slots_to_coeffs.log_q)
        + sum_bits(layout.modup.log_q)
        + sum_bits(layout.coeffs_to_slots.log_q)
        + sum_bits(layout.evalmod.log_q)
        + sum_bits(layout.output.log_q);
    ok = ok && layout.total_log_q == expected_q;
    ok = ok && layout.total_log_p == 183;
    ok = ok && layout.minimum_input_level == layout.slots_to_coeffs.levels
        + layout.modup.levels
        + layout.coeffs_to_slots.levels
        + layout.evalmod.levels;
    ok = ok && layout.expected_output_level == layout.output.levels - 1;
    ok = ok && layout.output.levels >= 6;

    bool invalid_threw = false;
    try {
        (void)m2424::make_lattigo_like_bootstrap_layout_v2(15, 15, 1e-9, 3);
    } catch (const std::invalid_argument&) {
        invalid_threw = true;
    }
    ok = ok && invalid_threw;

    invalid_threw = false;
    try {
        (void)m2424::make_lattigo_like_bootstrap_layout_v2(15, 2, -1.0, 3);
    } catch (const std::invalid_argument&) {
        invalid_threw = true;
    }
    ok = ok && invalid_threw;

    if (!ok) {
        std::cerr << "bootstrap layout v2 contract failed\n";
        return 1;
    }
    return 0;
}
