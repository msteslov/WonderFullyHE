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

struct CkksPlanningRequest {
    double target_error{1e-9};
    std::size_t multiplicative_depth{1};
    std::size_t slots{1};
    int security_bits{128};
    ParameterOptimizeFor optimize_for{ParameterOptimizeFor::Speed};
    PlanningOperationProfile operation_profile{PlanningOperationProfile::BasicMulDepth};
    int calibrated_loss_bits_override{-1};
};

struct CkksPlanningResult {
    CkksProfile profile;
    int required_result_bits{};
    int calibrated_loss_bits{};
    int selected_work_bits{};
    int selected_scale_log2{};
    std::size_t selected_work_levels{};
    SecurityReport security;
};

int required_result_bits(double target_error);
int calibrated_loss_bits(PlanningOperationProfile profile);
CkksPlanningResult plan_ckks_parameters(const CkksPlanningRequest& request);

const char* to_string(ParameterOptimizeFor value) noexcept;
const char* to_string(PlanningOperationProfile value) noexcept;
std::string planning_result_summary(const CkksPlanningResult& result);

} // namespace m2424
