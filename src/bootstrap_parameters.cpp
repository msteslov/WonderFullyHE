#include "m2424/bootstrap_plan.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <stdexcept>

namespace m2424 {

const char* to_string(BootstrapParameterPlanningStatus status) noexcept {
    switch (status) {
    case BootstrapParameterPlanningStatus::Ready:
        return "ready";
    case BootstrapParameterPlanningStatus::BlockedByUnsupportedSlots:
        return "blocked_by_unsupported_slots";
    case BootstrapParameterPlanningStatus::BlockedBySecurityBudget:
        return "blocked_by_security_budget";
    case BootstrapParameterPlanningStatus::BlockedByEvalModCapacity:
        return "blocked_by_evalmod_capacity";
    case BootstrapParameterPlanningStatus::BlockedByScaleBudget:
        return "blocked_by_scale_budget";
    }
    return "unknown";
}

namespace {

bool supported_transform_slots(BootstrapTransformBackend backend, std::size_t slots) {
    if (backend == BootstrapTransformBackend::DenseDiagonal) {
        return slots != 0 && (slots & (slots - 1)) == 0;
    }
    return slots == 4 || slots == 8 || slots == 16;
}

BootstrapParameterPlanningStatus map_layout_status(BootstrapLayoutPlanningStatus status) {
    switch (status) {
    case BootstrapLayoutPlanningStatus::Ready:
        return BootstrapParameterPlanningStatus::Ready;
    case BootstrapLayoutPlanningStatus::BlockedBySecurityBudget:
        return BootstrapParameterPlanningStatus::BlockedBySecurityBudget;
    case BootstrapLayoutPlanningStatus::BlockedByEvalModCapacity:
        return BootstrapParameterPlanningStatus::BlockedByEvalModCapacity;
    case BootstrapLayoutPlanningStatus::BlockedByScaleBudget:
        return BootstrapParameterPlanningStatus::BlockedByScaleBudget;
    }
    return BootstrapParameterPlanningStatus::BlockedByScaleBudget;
}

} // namespace

BootstrapParametersResult plan_bootstrap_parameters(const BootstrapParametersRequest& request) {
    if (request.slots == 0 || (request.slots & (request.slots - 1)) != 0) {
        throw std::invalid_argument("bootstrap parameters slots must be a non-zero power of two");
    }
    if (request.poly_modulus_degree == 0 || (request.poly_modulus_degree & (request.poly_modulus_degree - 1)) != 0) {
        throw std::invalid_argument("bootstrap parameters poly_modulus_degree must be a non-zero power of two");
    }
    if (request.slots > request.poly_modulus_degree / 2) {
        throw std::invalid_argument("bootstrap parameters slots exceed CKKS slot count");
    }
    if (!std::isfinite(request.max_abs_after_coeff_to_slot_log2)) {
        throw std::invalid_argument("bootstrap parameters max_abs_after_coeff_to_slot_log2 must be finite");
    }
    if (!std::isfinite(request.transform_output_scale_log2) || request.transform_output_scale_log2 <= 0.0
        || !std::isfinite(request.normalization_plain_scale_log2) || request.normalization_plain_scale_log2 <= 0.0
        || !std::isfinite(request.target_scale_log2) || request.target_scale_log2 <= 0.0
        || !std::isfinite(request.evalmod_capacity_margin_log2) || request.evalmod_capacity_margin_log2 < 0.0) {
        throw std::invalid_argument("bootstrap parameters scale fields are invalid");
    }

    BootstrapParametersResult result;
    result.pipeline = request.circuit_order == BootstrapCircuitOrder::SlotsToCoeffsFirst
        ? make_scalable_bootstrap_plan(request.slots)
        : make_research_bootstrap_plan(request.slots);
    result.pipeline.circuit_order = request.circuit_order;
    result.pipeline.transform_backend = request.transform_backend;
    result.pipeline.plain_scale_log2 = request.normalization_plain_scale_log2;
    result.pipeline.evalmod_degree = request.mod1_model.type == BootstrapMod1Type::LegacySineP3
        ? EvalModDegree::P3
        : EvalModDegree::P3DoubleAngle;
    result.pipeline.active_gate = BootstrapPipelineGate::FullRefresh;

    if (!supported_transform_slots(request.transform_backend, request.slots)) {
        result.status = BootstrapParameterPlanningStatus::BlockedByUnsupportedSlots;
        result.blocker = "transform_backend_does_not_support_slots";
        return result;
    }

    result.rotation_steps = bootstrap_plan_rotation_steps(result.pipeline);

    BootstrapLayoutPlanningRequest layout_request;
    layout_request.slots = request.slots;
    layout_request.poly_modulus_degree = request.poly_modulus_degree;
    layout_request.security_level = request.security_level;
    layout_request.transform_backend = request.transform_backend;
    layout_request.max_abs_after_coeff_to_slot_log2 = request.max_abs_after_coeff_to_slot_log2;
    layout_request.transform_output_scale_log2 = request.transform_output_scale_log2;
    layout_request.normalization_plain_scale_log2 = request.normalization_plain_scale_log2;
    layout_request.target_scale_log2 = request.target_scale_log2;
    layout_request.evalmod_capacity_margin_log2 = request.evalmod_capacity_margin_log2;
    layout_request.residual_levels = request.residual_levels;
    layout_request.first_mod_bits = request.first_mod_bits;
    layout_request.middle_mod_bits = request.middle_mod_bits;
    layout_request.last_mod_bits = request.last_mod_bits;
    layout_request.mod1_model = request.mod1_model;

    result.layout = plan_bootstrap_layout(layout_request);
    result.status = map_layout_status(result.layout.status);
    result.blocker = result.layout.blocker;
    result.profile = result.layout.profile;
    result.consumed_levels = result.layout.total_levels >= request.residual_levels
        ? result.layout.total_levels - request.residual_levels
        : result.layout.total_levels;
    result.required_input_chain_index = result.layout.total_levels;
    result.expected_output_chain_index = request.residual_levels;
    result.expected_output_scale_log2 = request.target_scale_log2;

    if (result.status == BootstrapParameterPlanningStatus::Ready) {
        if (result.rotation_steps.empty()) {
            result.status = BootstrapParameterPlanningStatus::BlockedByUnsupportedSlots;
            result.blocker = "bootstrap_rotation_steps_empty";
        } else {
            result.blocker = "none";
        }
    }
    return result;
}

} // namespace m2424
