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
    require(result.estimated_precision_bits == 31, "planner estimated precision mismatch");
    require(result.estimated_abs_error_bound <= request.target_error, "planner should pass target error");
    require(result.passes_target_error, "planner pass flag mismatch");
    require(result.selected_work_levels == 2, "planner depth mismatch");
    require(result.profile.polyModulusDegree == 8192, "planner should select minimal tc128 N for this request");
    require(result.profile.slots == 4096, "planner slots mismatch");
    require(result.profile.coeffModulusBits.size() == 4, "planner chain length mismatch");
    require(result.profile.coeffModulusBits[0] == 60, "planner first modulus mismatch");
    require(result.profile.coeffModulusBits[1] == 45, "planner work modulus mismatch");
    require(result.profile.coeffModulusBits[2] == 45, "planner work modulus mismatch");
    require(result.profile.coeffModulusBits[3] == 60, "planner last modulus mismatch");
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
    require(result.estimated_precision_bits == 36, "conservative estimated precision mismatch");
    require(result.passes_target_error, "conservative plan should pass target error");
    require(result.profile.polyModulusDegree == 16384, "50-bit depth-2 plan should move to minimal secure N");
    require(result.security.passes_tc128, "conservative plan should pass tc128");
}

void test_security_and_slots_drive_degree() {
    m2424::CkksPlanningRequest request;
    request.target_error = 1e-9;
    request.multiplicative_depth = 2;
    request.slots = 8192;
    request.security_bits = 128;

    const auto result = m2424::plan_ckks_parameters(request);
    require(result.profile.polyModulusDegree == 16384, "slots should drive N to 16384");
    require(result.profile.slots == 8192, "planner should preserve requested slots");
}

void test_operation_budget_drives_loss_and_levels() {
    m2424::CkksOperationBudget budget;
    budget.additions = 3;
    budget.ciphertext_muls = 1;
    budget.plaintext_mul_rescales = 1;
    budget.mod_switches = 2;
    budget.rotations = 4;
    budget.linear_transforms = 1;

    require(m2424::estimated_level_budget(budget) == 3, "operation budget level estimate mismatch");
    require(m2424::calibrated_loss_bits(budget) == 23, "operation budget loss estimate mismatch");

    m2424::CkksPlanningRequest request;
    request.target_error = 1e-9;
    request.multiplicative_depth = 1;
    request.slots = 4096;
    request.use_operation_budget = true;
    request.operation_budget = budget;

    const auto result = m2424::plan_ckks_parameters(request);
    require(result.calibrated_loss_bits == 23, "planner should use explicit operation budget loss");
    require(result.selected_work_levels == 3, "planner should use explicit operation budget levels");
    require(result.selected_work_bits == 55, "operation budget should drive work bits");
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
        test_operation_budget_drives_loss_and_levels();
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
