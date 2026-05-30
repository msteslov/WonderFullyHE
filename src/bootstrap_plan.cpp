#include "m2424/bootstrap_plan.hpp"

#include "m2424/bootstrap_prototype.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>

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

const char* to_string(BootstrapScaleDesignStatus status) noexcept {
    switch (status) {
    case BootstrapScaleDesignStatus::PeriodModelBlocked:
        return "period_model_blocked";
    case BootstrapScaleDesignStatus::ScaleStrategyBlocked:
        return "scale_strategy_blocked";
    case BootstrapScaleDesignStatus::EvalModCapacityBlocked:
        return "evalmod_capacity_blocked";
    case BootstrapScaleDesignStatus::ReadyForEvalModP3:
        return "ready_for_evalmod_p3";
    }
    return "unknown";
}

const char* to_string(BootstrapRefreshPlanningStatus status) noexcept {
    switch (status) {
    case BootstrapRefreshPlanningStatus::ComputeFitsWithoutRefresh:
        return "compute_fits_without_refresh";
    case BootstrapRefreshPlanningStatus::RefreshRequired:
        return "refresh_required";
    case BootstrapRefreshPlanningStatus::RefreshPlanBlocked:
        return "refresh_plan_blocked";
    }
    return "unknown";
}

std::vector<int> active_coeff_modulus_bits(const CkksProfile& profile, const CipherInfo& info) {
    const std::size_t active_size = std::min(info.coeff_modulus_size, profile.coeff_modulus_bits.size());
    return {profile.coeff_modulus_bits.begin(),
            profile.coeff_modulus_bits.begin() + static_cast<std::ptrdiff_t>(active_size)};
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
    plan.circuit_order = BootstrapCircuitOrder::ModRaiseFirst;
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
    plan.circuit_order = BootstrapCircuitOrder::SlotsToCoeffsFirst;
    plan.transform_backend = BootstrapTransformBackend::FftLike;
    plan.active_gate = BootstrapPipelineGate::Scaling;
    return plan;
}

std::vector<int> bootstrap_plan_rotation_steps(const BootstrapPipelinePlan& plan) {
    validate_slots(plan.slots);
    if (plan.transform_backend == BootstrapTransformBackend::DenseDiagonal) {
        return BootstrapPrototype::required_rotation_steps(plan.slots);
    }
    auto coeff_to_slot = make_bootstrap_dft_plan(plan.slots,
                                                 BootstrapDftType::HomomorphicDecode,
                                                 plan.plain_scale_log2);
    auto slot_to_coeff = make_bootstrap_dft_plan(plan.slots,
                                                 BootstrapDftType::HomomorphicEncode,
                                                 plan.plain_scale_log2);
    auto steps = coeff_to_slot.rotation_steps();
    auto inverse_steps = slot_to_coeff.rotation_steps();
    steps.insert(steps.end(), inverse_steps.begin(), inverse_steps.end());
    std::sort(steps.begin(), steps.end());
    steps.erase(std::unique(steps.begin(), steps.end()), steps.end());
    return steps;
}

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
                                                 double evalmod_capacity_margin_log2) {
    if (!std::isfinite(period_log2)) {
        throw std::invalid_argument("bootstrap scale design period log2 must be finite");
    }
    if (!std::isfinite(max_abs_before_normalization) || max_abs_before_normalization < 0.0) {
        throw std::invalid_argument("max_abs_before_normalization must be finite and non-negative");
    }

    BootstrapScaleDesign design;
    design.period_mode = period_mode;
    design.manual_period_log2 = manual_period_log2;
    design.period_log2 = period_log2;
    design.normalization_strategy = normalization_strategy;
    design.plain_scale_log2 = plain_scale_log2;
    design.target_scale_log2 = target_scale_log2;
    design.evalmod_degree = evalmod_degree;

    const double factor_log2 = -period_log2;
    const double expected_max_abs_after_normalization =
        max_abs_before_normalization * std::exp2(factor_log2);
    design.magnitude_ok =
        expected_max_abs_after_normalization <= EvalModPolynomial::approximation_bound;
    if (!design.magnitude_ok) {
        design.status = BootstrapScaleDesignStatus::PeriodModelBlocked;
        design.blocker = to_string(design.status);
        return design;
    }

    design.scale_plan = plan_bootstrap_scale_strategy(active_coeff_modulus_bits,
                                                      start_info,
                                                      factor_log2,
                                                      plain_scale_log2,
                                                      target_scale_log2,
                                                      min_chain_remaining);
    design.scale_strategy_ok = design.scale_plan.feasible;
    design.required_levels = design.scale_plan.total_levels_needed;
    design.chain_remaining_after_strategy = design.scale_plan.start_chain_index >= design.required_levels
        ? design.scale_plan.start_chain_index - design.required_levels
        : 0;
    if (!design.scale_strategy_ok) {
        design.status = BootstrapScaleDesignStatus::ScaleStrategyBlocked;
        design.blocker = design.scale_plan.blocker;
        return design;
    }
    if (evalmod_degree != EvalModDegree::P3) {
        design.status = BootstrapScaleDesignStatus::EvalModCapacityBlocked;
        design.blocker = "only_p3_design_gate_implemented";
        return design;
    }

    design.evalmod_capacity = plan_evalmod_first_multiply_capacity(active_coeff_modulus_bits,
                                                                   design.scale_plan,
                                                                   evalmod_capacity_margin_log2);
    design.evalmod_capacity_ok = design.evalmod_capacity.first_multiply_ready;
    if (!design.evalmod_capacity_ok) {
        design.status = BootstrapScaleDesignStatus::EvalModCapacityBlocked;
        design.blocker = design.evalmod_capacity.blocker;
        return design;
    }

    design.status = BootstrapScaleDesignStatus::ReadyForEvalModP3;
    design.blocker = "none";
    return design;
}

BootstrapPeriodFeasibilityWindow bootstrap_period_feasibility_window(
    double active_modulus_log2,
    double start_scale_log2,
    double target_scale_log2,
    double max_abs_before_normalization,
    double evalmod_bound,
    double evalmod_capacity_margin_log2) {
    if (!std::isfinite(active_modulus_log2)) {
        throw std::invalid_argument("active modulus log2 must be finite");
    }
    if (!std::isfinite(start_scale_log2) || start_scale_log2 <= 0.0) {
        throw std::invalid_argument("start scale log2 must be positive and finite");
    }
    if (!std::isfinite(target_scale_log2) || target_scale_log2 <= 0.0) {
        throw std::invalid_argument("target scale log2 must be positive and finite");
    }
    if (!std::isfinite(max_abs_before_normalization) || max_abs_before_normalization < 0.0) {
        throw std::invalid_argument("max_abs_before_normalization must be finite and non-negative");
    }
    if (!std::isfinite(evalmod_bound) || evalmod_bound <= 0.0) {
        throw std::invalid_argument("evalmod bound must be positive and finite");
    }
    if (!std::isfinite(evalmod_capacity_margin_log2) || evalmod_capacity_margin_log2 < 0.0) {
        throw std::invalid_argument("evalmod capacity margin log2 must be finite and non-negative");
    }

    BootstrapPeriodFeasibilityWindow window;
    if (max_abs_before_normalization == 0.0) {
        window.min_period_for_magnitude_log2 = 0.0;
    } else {
        window.min_period_for_magnitude_log2 =
            std::max(0.0, std::log2(max_abs_before_normalization / evalmod_bound));
    }
    window.max_period_for_evalmod_capacity_log2 =
        active_modulus_log2 - start_scale_log2 - target_scale_log2 - evalmod_capacity_margin_log2;
    window.margin_log2 =
        window.max_period_for_evalmod_capacity_log2 - window.min_period_for_magnitude_log2;
    window.possible = window.margin_log2 >= 0.0;
    return window;
}

BootstrapRefreshPlanningResult plan_bootstrap_refresh(const BootstrapRefreshPlanningRequest& request) {
    if (request.slots == 0) {
        throw std::invalid_argument("bootstrap refresh planner slots must be positive");
    }
    if (request.current_info.chain_index < request.min_chain_remaining_after_compute) {
        throw std::invalid_argument("min_chain_remaining_after_compute exceeds current chain index");
    }

    BootstrapRefreshPlanningResult result;
    result.required_compute_levels = estimated_level_budget(request.operation_budget);
    result.available_compute_levels =
        request.current_info.chain_index - request.min_chain_remaining_after_compute;

    CkksPlanningRequest parameter_request;
    parameter_request.target_error = request.target_error;
    parameter_request.multiplicative_depth = std::max<std::size_t>(1, result.required_compute_levels);
    parameter_request.slots = request.slots;
    parameter_request.security_bits = request.security_bits;
    parameter_request.optimize_for = request.optimize_for;
    parameter_request.use_operation_budget = true;
    parameter_request.operation_budget = request.operation_budget;

    try {
        result.parameter_plan = plan_ckks_parameters(parameter_request);
        result.parameter_plan_ok = result.parameter_plan.passes_target_error;
    } catch (const std::exception& error) {
        result.status = BootstrapRefreshPlanningStatus::RefreshPlanBlocked;
        result.needs_refresh = true;
        result.blocker = error.what();
        return result;
    }

    if (!result.parameter_plan_ok) {
        result.status = BootstrapRefreshPlanningStatus::RefreshPlanBlocked;
        result.needs_refresh = true;
        result.blocker = "parameter_plan_failed_target_error";
        return result;
    }

    if (result.required_compute_levels <= result.available_compute_levels) {
        result.status = BootstrapRefreshPlanningStatus::ComputeFitsWithoutRefresh;
        result.needs_refresh = false;
        result.blocker = "none";
        return result;
    }

    result.status = BootstrapRefreshPlanningStatus::RefreshRequired;
    result.needs_refresh = true;
    result.blocker = "insufficient_chain_for_next_budget";
    return result;
}

BootstrapScaleDesign plan_bootstrap_refresh_scale_gate(const BootstrapRefreshScaleGateRequest& request) {
    const double period_log2 = bootstrap_period_log2(
        request.period_mode,
        request.manual_period_log2,
        request.profile.coeff_modulus_bits,
        request.before_mod_raise,
        request.after_mod_raise);
    return make_bootstrap_scale_design(
        request.period_mode,
        request.manual_period_log2,
        period_log2,
        request.normalization_strategy,
        request.plain_scale_log2,
        request.target_scale_log2,
        request.evalmod_degree,
        active_coeff_modulus_bits(request.profile, request.slot_domain_info),
        request.slot_domain_info,
        request.max_abs_before_normalization,
        request.min_chain_remaining,
        request.evalmod_capacity_margin_log2);
}

BootstrapRefreshScaleGateSearchResult search_bootstrap_refresh_scale_gate(
    const BootstrapRefreshScaleGateSearchRequest& request) {
    if (request.period_modes.empty()) {
        throw std::invalid_argument("bootstrap scale gate search period modes must not be empty");
    }
    if (request.plain_scale_log2_values.empty()) {
        throw std::invalid_argument("bootstrap scale gate search plain scale values must not be empty");
    }
    if (request.target_scale_log2_values.empty()) {
        throw std::invalid_argument("bootstrap scale gate search target scale values must not be empty");
    }

    const auto active_bits = active_coeff_modulus_bits(request.profile, request.slot_domain_info);
    BootstrapRefreshScaleGateSearchResult result;
    double best_score = -std::numeric_limits<double>::infinity();

    const auto consider = [&](BootstrapPeriodMode period_mode,
                              double manual_period_log2,
                              double period_log2,
                              double plain_scale_log2,
                              double target_scale_log2) {
        auto design = make_bootstrap_scale_design(
            period_mode,
            manual_period_log2,
            period_log2,
            request.normalization_strategy,
            plain_scale_log2,
            target_scale_log2,
            request.evalmod_degree,
            active_bits,
            request.slot_domain_info,
            request.max_abs_before_normalization,
            request.min_chain_remaining,
            request.evalmod_capacity_margin_log2);
        ++result.candidates;

        const bool ready = design.status == BootstrapScaleDesignStatus::ReadyForEvalModP3;
        if (ready) {
            ++result.ready_candidates;
        }

        double score = 0.0;
        if (ready) {
            score += 1'000'000.0;
            score += static_cast<double>(design.chain_remaining_after_strategy) * 1'000.0;
            score += design.evalmod_capacity.margin_log2;
            score -= std::abs(design.target_scale_log2 - 40.0);
            score -= 0.01 * design.plain_scale_log2;
        } else {
            score += design.magnitude_ok ? 10'000.0 : 0.0;
            score += design.scale_strategy_ok ? 1'000.0 : 0.0;
            score += design.evalmod_capacity_ok ? 100.0 : 0.0;
            score -= static_cast<double>(design.required_levels) * 10.0;
            score += static_cast<double>(design.chain_remaining_after_strategy);
            score -= 0.01 * design.plain_scale_log2;
        }

        if (result.candidates == 1 || score > best_score) {
            best_score = score;
            result.best_design = std::move(design);
            result.ready = ready;
        } else if (ready) {
            result.ready = true;
        }
    };

    for (BootstrapPeriodMode period_mode : request.period_modes) {
        if (period_mode == BootstrapPeriodMode::ManualPowerOfTwo) {
            if (request.manual_period_log2_values.empty()) {
                throw std::invalid_argument("manual period search requires manual period values");
            }
            for (double manual_period_log2 : request.manual_period_log2_values) {
                for (double plain_scale_log2 : request.plain_scale_log2_values) {
                    for (double target_scale_log2 : request.target_scale_log2_values) {
                        consider(period_mode,
                                 manual_period_log2,
                                 manual_period_log2,
                                 plain_scale_log2,
                                 target_scale_log2);
                    }
                }
            }
            continue;
        }

        const double period_log2 = bootstrap_period_log2(
            period_mode,
            0.0,
            request.profile.coeff_modulus_bits,
            request.before_mod_raise,
            request.after_mod_raise);
        for (double plain_scale_log2 : request.plain_scale_log2_values) {
            for (double target_scale_log2 : request.target_scale_log2_values) {
                consider(period_mode,
                         0.0,
                         period_log2,
                         plain_scale_log2,
                         target_scale_log2);
            }
        }
    }

    return result;
}

} // namespace m2424
