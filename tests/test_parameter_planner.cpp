#include "m2424/parameter_planner.hpp"

#include <cmath>
#include <cstdio>
#include <exception>
#include <stdexcept>
#include <string>

namespace {

void require(bool condition, const char* message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

void test_target_bits() {
    require(m2424::required_result_bits(1e-9) == 30, "1e-9 should require 30 result bits");
    require(m2424::required_result_bits(1e-6) == 20, "1e-6 should require 20 result bits");
}

void test_speed_plan_for_one_e_minus_nine() {
    m2424::CkksPlanningRequest request;
    request.target_error = 1e-9;
    request.multiplicative_depth = 2;
    request.slots = 4096;
    request.security_bits = 128;
    request.optimize_for = m2424::ParameterOptimizeFor::Speed;

    const auto result = m2424::plan_ckks_parameters(request);
    require(result.required_result_bits == 30, "planner result bits mismatch");
    require(result.calibrated_loss_bits == 14, "planner calibrated loss mismatch");
    require(result.selected_work_bits == 45, "speed plan should select 45-bit work moduli");
    require(result.selected_scale_log2 == 45, "speed plan should align scale and work bits");
    require(result.selected_work_levels == 2, "planner depth mismatch");
    require(result.profile.poly_modulus_degree == 8192, "planner should select minimal tc128 N for this request");
    require(result.profile.slots == 4096, "planner slots mismatch");
    require(result.profile.coeff_modulus_bits.size() == 4, "planner chain length mismatch");
    require(result.profile.coeff_modulus_bits[0] == 60, "planner first modulus mismatch");
    require(result.profile.coeff_modulus_bits[1] == 45, "planner work modulus mismatch");
    require(result.profile.coeff_modulus_bits[2] == 45, "planner work modulus mismatch");
    require(result.profile.coeff_modulus_bits[3] == 60, "planner last modulus mismatch");
    require(std::abs(std::log2(result.profile.scale) - 45.0) < 1e-9, "planner scale mismatch");
    require(result.security.passes_tc128, "planner security report should pass tc128");
}

void test_conservative_plan() {
    m2424::CkksPlanningRequest request;
    request.target_error = 1e-9;
    request.multiplicative_depth = 2;
    request.slots = 4096;
    request.optimize_for = m2424::ParameterOptimizeFor::Conservative;

    const auto result = m2424::plan_ckks_parameters(request);
    require(result.selected_work_bits == 50, "conservative plan should select at least 50-bit work moduli");
    require(result.profile.poly_modulus_degree == 16384, "50-bit depth-2 plan should move to minimal secure N");
    require(result.security.passes_tc128, "conservative plan should pass tc128");
}

void test_security_and_slots_drive_degree() {
    m2424::CkksPlanningRequest request;
    request.target_error = 1e-9;
    request.multiplicative_depth = 2;
    request.slots = 8192;
    request.security_bits = 128;

    const auto result = m2424::plan_ckks_parameters(request);
    require(result.profile.poly_modulus_degree == 16384, "slots should drive N to 16384");
    require(result.profile.slots == 8192, "planner should preserve requested slots");
}

void test_invalid_inputs() {
    bool threw = false;
    try {
        (void)m2424::required_result_bits(1.0);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    require(threw, "invalid target_error should throw");

    threw = false;
    try {
        m2424::CkksPlanningRequest request;
        request.multiplicative_depth = 0;
        (void)m2424::plan_ckks_parameters(request);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    require(threw, "zero depth should throw");
}

} // namespace

int main() {
    bool ok = true;
    try {
        test_target_bits();
        test_speed_plan_for_one_e_minus_nine();
        test_conservative_plan();
        test_security_and_slots_drive_degree();
        test_invalid_inputs();
    } catch (const std::exception& error) {
        ok = false;
        std::printf("[test_parameter_planner] FAIL: %s\n", error.what());
    }
    if (ok) {
        std::printf("[test_parameter_planner] PASS\n");
    }
    return ok ? 0 : 1;
}
