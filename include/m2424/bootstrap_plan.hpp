#pragma once

#include "m2424/bootstrap_dft.hpp"
#include "m2424/bootstrap_scaling.hpp"
#include "m2424/eval_mod.hpp"
#include "m2424/parameter_planner.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace m2424 {

enum class BootstrapValueDomain {
    Coefficient,
    RaisedCoefficient,
    Slot,
    NormalizedSlot,
    ModularReducedSlot,
    RefreshedCoefficient
};

enum class BootstrapScalingStrategy {
    PlainMultiplyRescale,
    DecomposedPlainMultiplyRescale
};

enum class BootstrapPipelineGate {
    ModRaiseStructural,
    Scaling,
    EvalMod,
    FullRefresh
};

enum class BootstrapScaleDesignStatus {
    PeriodModelBlocked,
    ScaleStrategyBlocked,
    EvalModCapacityBlocked,
    ReadyForEvalModP3
};

enum class BootstrapRefreshPlanningStatus {
    ComputeFitsWithoutRefresh,
    RefreshRequired,
    RefreshPlanBlocked
};

struct BootstrapStageSpec {
    std::string name;
    BootstrapValueDomain input_domain{};
    BootstrapValueDomain output_domain{};
    bool structural_only{};
    bool value_preserving{};
    bool consumes_level{};
    std::string contract;
};

struct BootstrapPipelinePlan {
    std::size_t slots{};
    BootstrapCircuitOrder circuit_order{BootstrapCircuitOrder::ModRaiseFirst};
    BootstrapTransformBackend transform_backend{BootstrapTransformBackend::DenseDiagonal};
    BootstrapScalingStrategy scaling_strategy{BootstrapScalingStrategy::DecomposedPlainMultiplyRescale};
    BootstrapPeriodMode period_mode{BootstrapPeriodMode::ManualPowerOfTwo};
    double manual_period_log2{220.0};
    double plain_scale_log2{240.0};
    EvalModDegree evalmod_degree{EvalModDegree::P7};
    BootstrapPipelineGate active_gate{BootstrapPipelineGate::Scaling};
    std::vector<BootstrapStageSpec> stages;
};

struct BootstrapScaleDesign {
    BootstrapPeriodMode period_mode{BootstrapPeriodMode::ManualPowerOfTwo};
    BootstrapScalingStrategy normalization_strategy{BootstrapScalingStrategy::DecomposedPlainMultiplyRescale};
    EvalModDegree evalmod_degree{EvalModDegree::P3};
    double manual_period_log2{};
    double period_log2{};
    double coeff_to_slot_prescale_log2{};
    double plain_scale_log2{};
    double target_scale_log2{};
    std::size_t required_levels{};
    std::size_t chain_remaining_after_strategy{};
    bool magnitude_ok{};
    bool scale_strategy_ok{};
    bool evalmod_capacity_ok{};
    BootstrapScaleDesignStatus status{BootstrapScaleDesignStatus::ScaleStrategyBlocked};
    std::string blocker;
    BootstrapScaleStrategyPlan scale_plan;
    BootstrapEvalModCapacityPlan evalmod_capacity;
};

struct BootstrapPeriodFeasibilityWindow {
    double min_period_for_magnitude_log2{};
    double max_period_for_evalmod_capacity_log2{};
    double margin_log2{};
    bool possible{};
};

struct BootstrapRefreshPlanningRequest {
    CipherInfo current_info;
    CkksOperationBudget operation_budget;
    double target_error{1e-9};
    std::size_t slots{1};
    int security_bits{128};
    ParameterOptimizeFor optimize_for{ParameterOptimizeFor::Speed};
    std::size_t min_chain_remaining_after_compute{};
};

struct BootstrapRefreshPlanningResult {
    BootstrapRefreshPlanningStatus status{BootstrapRefreshPlanningStatus::RefreshPlanBlocked};
    std::size_t required_compute_levels{};
    std::size_t available_compute_levels{};
    bool needs_refresh{};
    bool parameter_plan_ok{};
    std::string blocker;
    CkksPlanningResult parameter_plan;
};

struct BootstrapRefreshScaleGateRequest {
    CkksProfile profile;
    CipherInfo before_mod_raise;
    CipherInfo after_mod_raise;
    CipherInfo slot_domain_info;
    double max_abs_before_normalization{};
    BootstrapPeriodMode period_mode{BootstrapPeriodMode::TotalCoeffModulus};
    double manual_period_log2{};
    BootstrapScalingStrategy normalization_strategy{BootstrapScalingStrategy::DecomposedPlainMultiplyRescale};
    double plain_scale_log2{40.0};
    double target_scale_log2{60.0};
    EvalModDegree evalmod_degree{EvalModDegree::P3};
    std::size_t min_chain_remaining{};
    double evalmod_capacity_margin_log2{2.0};
};

struct BootstrapRefreshScaleGateSearchRequest {
    CkksProfile profile;
    CipherInfo before_mod_raise;
    CipherInfo after_mod_raise;
    CipherInfo slot_domain_info;
    double max_abs_before_normalization{};
    std::vector<BootstrapPeriodMode> period_modes{BootstrapPeriodMode::TotalCoeffModulus};
    std::vector<double> manual_period_log2_values;
    std::vector<double> coeff_to_slot_prescale_log2_values{0.0};
    std::vector<double> plain_scale_log2_values{40.0, 50.0, 60.0, 80.0, 100.0, 120.0, 160.0, 200.0, 240.0};
    std::vector<double> target_scale_log2_values{30.0, 40.0, 50.0, 60.0};
    BootstrapScalingStrategy normalization_strategy{BootstrapScalingStrategy::DecomposedPlainMultiplyRescale};
    EvalModDegree evalmod_degree{EvalModDegree::P3};
    std::size_t min_chain_remaining{};
    double evalmod_capacity_margin_log2{2.0};
};

struct BootstrapRefreshScaleGateSearchResult {
    BootstrapScaleDesign best_design;
    std::size_t candidates{};
    std::size_t ready_candidates{};
    bool ready{};
};

const char* to_string(BootstrapValueDomain domain) noexcept;
const char* to_string(BootstrapTransformBackend backend) noexcept;
const char* to_string(BootstrapScalingStrategy strategy) noexcept;
const char* to_string(BootstrapPipelineGate gate) noexcept;
const char* to_string(BootstrapScaleDesignStatus status) noexcept;
const char* to_string(BootstrapRefreshPlanningStatus status) noexcept;

BootstrapPipelinePlan make_research_bootstrap_plan(std::size_t slots);
BootstrapPipelinePlan make_scalable_bootstrap_plan(std::size_t slots);
std::vector<int> bootstrap_plan_rotation_steps(const BootstrapPipelinePlan& plan);
std::vector<int> active_coeff_modulus_bits(const CkksProfile& profile, const CipherInfo& info);

BootstrapScaleDesign make_bootstrap_scale_design(BootstrapPeriodMode period_mode,
                                                 double manual_period_log2,
                                                 double period_log2,
                                                 BootstrapScalingStrategy normalization_strategy,
                                                 double plain_scale_log2,
                                                 double target_scale_log2,
                                                 EvalModDegree evalmod_degree,
                                                 const std::vector<int>& active_coeff_modulus_bits,
                                                 const CipherInfo& start_info,
                                                 double max_abs_before_normalization,
                                                 std::size_t min_chain_remaining,
                                                 double evalmod_capacity_margin_log2);

BootstrapPeriodFeasibilityWindow bootstrap_period_feasibility_window(
    double active_modulus_log2,
    double start_scale_log2,
    double target_scale_log2,
    double max_abs_before_normalization,
    double evalmod_bound,
    double evalmod_capacity_margin_log2);

BootstrapRefreshPlanningResult plan_bootstrap_refresh(const BootstrapRefreshPlanningRequest& request);
BootstrapScaleDesign plan_bootstrap_refresh_scale_gate(const BootstrapRefreshScaleGateRequest& request);
BootstrapRefreshScaleGateSearchResult search_bootstrap_refresh_scale_gate(
    const BootstrapRefreshScaleGateSearchRequest& request);

} // namespace m2424
