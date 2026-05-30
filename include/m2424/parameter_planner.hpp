#pragma once

#include "m2424/seal_adapter.hpp"
#include "m2424/security_report.hpp"

#include <cstddef>
#include <string>

namespace m2424 {

enum class ParameterOptimizeFor {
    Speed,
    Conservative
};

enum class PlanningOperationProfile {
    BasicMulDepth,
    LinearTransform,
    BootstrapRefresh
};

struct CkksOperationBudget {
    std::size_t additions{};
    std::size_t plaintext_additions{};
    std::size_t ciphertext_muls{};
    std::size_t plaintext_mul_rescales{};
    std::size_t rescale_to_next{};
    std::size_t mod_switches{};
    std::size_t rotations{};
    std::size_t linear_transforms{};
    std::size_t evalmod_p3{};
    std::size_t bootstrap_refreshes{};
};

struct CkksPlanningRequest {
    double target_error{1e-9};
    std::size_t multiplicative_depth{1};
    std::size_t slots{1};
    int security_bits{128};
    ParameterOptimizeFor optimize_for{ParameterOptimizeFor::Speed};
    PlanningOperationProfile operation_profile{PlanningOperationProfile::BasicMulDepth};
    CkksOperationBudget operation_budget{};
    bool use_operation_budget{false};
    int calibrated_loss_bits_override{-1};
};

struct CkksPlanningResult {
    CkksProfile profile;
    double target_error{};
    int required_result_bits{};
    int calibrated_loss_bits{};
    int selected_work_bits{};
    int selected_scale_log2{};
    int estimated_precision_bits{};
    double estimated_abs_error_bound{};
    bool passes_target_error{};
    std::size_t selected_work_levels{};
    SecurityReport security;
};

int required_result_bits(double target_error);
int calibrated_loss_bits(PlanningOperationProfile profile);
int calibrated_loss_bits(const CkksOperationBudget& budget);
std::size_t estimated_level_budget(const CkksOperationBudget& budget);
CkksPlanningResult plan_ckks_parameters(const CkksPlanningRequest& request);

const char* to_string(ParameterOptimizeFor value) noexcept;
const char* to_string(PlanningOperationProfile value) noexcept;
std::string planning_result_summary(const CkksPlanningResult& result);

} // namespace m2424
