#include "m2424/bootstrap_plan.hpp"

#include "m2424/bootstrap_prototype.hpp"

#include <stdexcept>

namespace m2424 {

const char* to_string(BootstrapValueDomain domain) noexcept {
    switch (domain) {
    case BootstrapValueDomain::Coefficient:
        return "Coefficient";
    case BootstrapValueDomain::RaisedCoefficient:
        return "RaisedCoefficient";
    case BootstrapValueDomain::Slot:
        return "Slot";
    case BootstrapValueDomain::NormalizedSlot:
        return "NormalizedSlot";
    case BootstrapValueDomain::ModularReducedSlot:
        return "ModularReducedSlot";
    case BootstrapValueDomain::RefreshedCoefficient:
        return "RefreshedCoefficient";
    }
    return "unknown";
}

const char* to_string(BootstrapTransformBackend backend) noexcept {
    switch (backend) {
    case BootstrapTransformBackend::DenseDiagonal:
        return "DenseDiagonal";
    case BootstrapTransformBackend::FftLike:
        return "FftLike";
    }
    return "unknown";
}

const char* to_string(BootstrapScalingStrategy strategy) noexcept {
    switch (strategy) {
    case BootstrapScalingStrategy::PlainMultiplyRescale:
        return "PlainMultiplyRescale";
    case BootstrapScalingStrategy::DecomposedPlainMultiplyRescale:
        return "DecomposedPlainMultiplyRescale";
    }
    return "unknown";
}

const char* to_string(BootstrapPipelineGate gate) noexcept {
    switch (gate) {
    case BootstrapPipelineGate::ModRaiseStructural:
        return "ModRaiseStructural";
    case BootstrapPipelineGate::Scaling:
        return "Scaling";
    case BootstrapPipelineGate::EvalMod:
        return "EvalMod";
    case BootstrapPipelineGate::FullRefresh:
        return "FullRefresh";
    }
    return "unknown";
}

namespace {

std::vector<BootstrapStageSpec> standard_stages() {
    return {
        {
            "ModRaise",
            BootstrapValueDomain::Coefficient,
            BootstrapValueDomain::RaisedCoefficient,
            true,
            false,
            false,
            "Structural RNS-base extension only; decoded value is diagnostic, not a correctness criterion."
        },
        {
            "CoeffToSlot",
            BootstrapValueDomain::RaisedCoefficient,
            BootstrapValueDomain::Slot,
            false,
            false,
            true,
            "Move raised coefficient representation into slots before modular reduction."
        },
        {
            "Normalization",
            BootstrapValueDomain::Slot,
            BootstrapValueDomain::NormalizedSlot,
            false,
            false,
            true,
            "Apply bootstrap-period and amplitude scaling through plaintext multiplication; tiny scalars may be decomposed."
        },
        {
            "EvalMod",
            BootstrapValueDomain::NormalizedSlot,
            BootstrapValueDomain::ModularReducedSlot,
            false,
            false,
            true,
            "Approximate modular reduction on normalized slots; degree fallback is explicit, never silent."
        },
        {
            "Denormalization",
            BootstrapValueDomain::ModularReducedSlot,
            BootstrapValueDomain::Slot,
            false,
            false,
            true,
            "Restore slot magnitude before returning to coefficient representation."
        },
        {
            "SlotToCoeff",
            BootstrapValueDomain::Slot,
            BootstrapValueDomain::RefreshedCoefficient,
            false,
            true,
            true,
            "Return refreshed slot values to coefficient representation."
        }
    };
}

void validate_slots(std::size_t slots) {
    if (slots == 0 || (slots & (slots - 1)) != 0) {
        throw std::invalid_argument("bootstrap plan slots must be a non-zero power of two");
    }
}

} // namespace

BootstrapPipelinePlan make_research_bootstrap_plan(std::size_t slots) {
    validate_slots(slots);
    BootstrapPipelinePlan plan;
    plan.slots = slots;
    plan.transform_backend = BootstrapTransformBackend::DenseDiagonal;
    plan.scaling_strategy = BootstrapScalingStrategy::DecomposedPlainMultiplyRescale;
    plan.period_mode = BootstrapPeriodMode::ManualPowerOfTwo;
    plan.manual_period_log2 = 220.0;
    plan.plain_scale_log2 = 240.0;
    plan.evalmod_degree = EvalModDegree::P7;
    plan.active_gate = BootstrapPipelineGate::Scaling;
    plan.stages = standard_stages();
    return plan;
}

BootstrapPipelinePlan make_scalable_bootstrap_plan(std::size_t slots) {
    validate_slots(slots);
    auto plan = make_research_bootstrap_plan(slots);
    plan.transform_backend = BootstrapTransformBackend::FftLike;
    plan.active_gate = BootstrapPipelineGate::Scaling;
    return plan;
}

std::vector<int> bootstrap_plan_rotation_steps(const BootstrapPipelinePlan& plan) {
    validate_slots(plan.slots);
    if (plan.transform_backend == BootstrapTransformBackend::DenseDiagonal) {
        return BootstrapPrototype::required_rotation_steps(plan.slots);
    }
    throw std::logic_error("FFT-like bootstrap transform backend is specified but not implemented yet");
}

} // namespace m2424
