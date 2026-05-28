#pragma once

#include "m2424/bootstrap_dft.hpp"
#include "m2424/bootstrap_scaling.hpp"
#include "m2424/eval_mod.hpp"

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

const char* to_string(BootstrapValueDomain domain) noexcept;
const char* to_string(BootstrapTransformBackend backend) noexcept;
const char* to_string(BootstrapScalingStrategy strategy) noexcept;
const char* to_string(BootstrapPipelineGate gate) noexcept;
const char* to_string(BootstrapScaleDesignStatus status) noexcept;

BootstrapPipelinePlan make_research_bootstrap_plan(std::size_t slots);
BootstrapPipelinePlan make_scalable_bootstrap_plan(std::size_t slots);
std::vector<int> bootstrap_plan_rotation_steps(const BootstrapPipelinePlan& plan);

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

} // namespace m2424
