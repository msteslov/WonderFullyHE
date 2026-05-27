#pragma once

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

enum class BootstrapTransformBackend {
    DenseDiagonal,
    FftLike
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
    BootstrapTransformBackend transform_backend{BootstrapTransformBackend::DenseDiagonal};
    BootstrapScalingStrategy scaling_strategy{BootstrapScalingStrategy::DecomposedPlainMultiplyRescale};
    BootstrapPeriodMode period_mode{BootstrapPeriodMode::ManualPowerOfTwo};
    double manual_period_log2{220.0};
    double plain_scale_log2{240.0};
    EvalModDegree evalmod_degree{EvalModDegree::P7};
    BootstrapPipelineGate active_gate{BootstrapPipelineGate::Scaling};
    std::vector<BootstrapStageSpec> stages;
};

const char* to_string(BootstrapValueDomain domain) noexcept;
const char* to_string(BootstrapTransformBackend backend) noexcept;
const char* to_string(BootstrapScalingStrategy strategy) noexcept;
const char* to_string(BootstrapPipelineGate gate) noexcept;

BootstrapPipelinePlan make_research_bootstrap_plan(std::size_t slots);
BootstrapPipelinePlan make_scalable_bootstrap_plan(std::size_t slots);
std::vector<int> bootstrap_plan_rotation_steps(const BootstrapPipelinePlan& plan);

} // namespace m2424
