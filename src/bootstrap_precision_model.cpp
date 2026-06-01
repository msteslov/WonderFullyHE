#include "m2424/bootstrap_precision_model.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace m2424 {
namespace {

constexpr const char* kFloorBlocker =
    "DFT precision is key-switch/rotation-noise-floor limited; increasing plaintext scale alone is insufficient.";

void validate_positive_finite(double value, const char* name) {
    if (!std::isfinite(value) || value <= 0.0) {
        throw std::invalid_argument(std::string(name) + " must be positive and finite");
    }
}

std::vector<DftPrecisionMeasurement> sorted_measurements(std::vector<DftPrecisionMeasurement> values) {
    if (values.size() < 2) {
        throw std::invalid_argument("DFT precision fit requires at least two measurements");
    }
    std::sort(values.begin(), values.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.transform_plain_scale_log2 < rhs.transform_plain_scale_log2;
    });
    return values;
}

double minimum_roundtrip_error(const std::vector<DftPrecisionMeasurement>& measurements) {
    double result = std::numeric_limits<double>::infinity();
    for (const auto& measurement : measurements) {
        result = std::min(result, measurement.roundtrip_error);
    }
    return result;
}

double best_scale_for_roundtrip(const std::vector<DftPrecisionMeasurement>& measurements) {
    return std::min_element(measurements.begin(), measurements.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.roundtrip_error < rhs.roundtrip_error;
    })->transform_plain_scale_log2;
}

} // namespace

BootstrapStageErrorBudget make_bootstrap_error_budget(double target_total_error,
                                                      std::size_t cycles) {
    validate_positive_finite(target_total_error, "target_total_error");
    if (cycles == 0) {
        throw std::invalid_argument("cycles must be positive");
    }

    BootstrapStageErrorBudget budget;
    budget.target_total_error = target_total_error;
    budget.cycles = cycles;
    budget.per_cycle_budget = target_total_error / (static_cast<double>(cycles) * 2.0);
    budget.dft_roundtrip_budget = budget.per_cycle_budget * 0.5;
    budget.evalmod_budget = budget.per_cycle_budget * 0.25;
    budget.modraise_budget = budget.per_cycle_budget * 0.25;
    return budget;
}

DftPrecisionFit fit_dft_precision_floor(const std::vector<DftPrecisionMeasurement>& measurements,
                                        double target_roundtrip_error) {
    validate_positive_finite(target_roundtrip_error, "target_roundtrip_error");
    const auto sorted = sorted_measurements(measurements);

    const auto& low = sorted.front();
    const auto& high = sorted.back();
    const double observed_best = minimum_roundtrip_error(sorted);
    const double best_scale = best_scale_for_roundtrip(sorted);

    double noise_floor = high.roundtrip_error;
    if (sorted.size() >= 2) {
        const auto& prev = sorted[sorted.size() - 2];
        const double high_scale_improvement = prev.roundtrip_error - high.roundtrip_error;
        const double relative_improvement =
            high_scale_improvement / std::max(prev.roundtrip_error, std::numeric_limits<double>::min());
        if (high_scale_improvement > 0.0 && relative_improvement < 0.35) {
            noise_floor = std::min(prev.roundtrip_error, high.roundtrip_error);
        } else {
            noise_floor = observed_best;
        }
    }

    const double low_quant = std::exp2(-low.transform_plain_scale_log2);
    const double high_quant = std::exp2(-high.transform_plain_scale_log2);
    const double denominator = low_quant - high_quant;
    double coefficient = 0.0;
    if (denominator > 0.0) {
        coefficient = std::max(0.0, (low.roundtrip_error - high.roundtrip_error) / denominator);
    }

    DftPrecisionFit fit;
    fit.quantization_coefficient = coefficient;
    fit.noise_floor = noise_floor;
    fit.predicted_best_error = std::max(noise_floor, observed_best);
    fit.best_plain_scale_log2 = best_scale;
    fit.floor_blocks_target = fit.predicted_best_error > target_roundtrip_error;
    if (fit.floor_blocks_target) {
        fit.blocker = kFloorBlocker;
    }
    return fit;
}

BootstrapDftCost estimate_bootstrap_dft_cost(std::size_t slots,
                                             BootstrapDftType type,
                                             double plain_scale_log2) {
    auto plan = make_bootstrap_dft_plan(slots, type, plain_scale_log2);
    BootstrapDftCost cost;
    cost.slots = slots;
    cost.type = type;
    cost.layer_count = plan.layers.size();
    for (const auto& layer : plan.layers) {
        cost.diagonal_term_count += layer.transform.terms().size();
        auto steps = layer.transform.rotation_steps();
        cost.rotation_count += steps.size();
        cost.plaintext_multiplication_count += layer.transform.terms().size();
        cost.rescale_count += 1;
        cost.rotation_steps.insert(cost.rotation_steps.end(), steps.begin(), steps.end());
    }
    std::sort(cost.rotation_steps.begin(), cost.rotation_steps.end());
    cost.rotation_steps.erase(std::unique(cost.rotation_steps.begin(), cost.rotation_steps.end()),
                              cost.rotation_steps.end());
    return cost;
}

BootstrapPrecisionPlanningResult plan_bootstrap_precision(
    const BootstrapPrecisionPlanningRequest& request,
    const std::vector<DftPrecisionMeasurement>& calibration) {
    BootstrapPrecisionPlanningResult result;
    result.budget = make_bootstrap_error_budget(request.target_total_error, request.cycles);
    if (request.slots == 0) {
        throw std::invalid_argument("slots must be positive");
    }
    if (request.candidate_profiles.empty()) {
        throw std::invalid_argument("candidate_profiles must not be empty");
    }
    if (request.candidate_transform_scales.empty()) {
        throw std::invalid_argument("candidate_transform_scales must not be empty");
    }

    result.dft_fit = fit_dft_precision_floor(calibration, result.budget.dft_roundtrip_budget);
    result.selected_profile = request.candidate_profiles.front();
    result.selected_transform_plain_scale_log2 = result.dft_fit.best_plain_scale_log2;
    result.feasible = !result.dft_fit.floor_blocks_target;
    if (!result.feasible) {
        result.blocker = result.dft_fit.blocker;
    }
    return result;
}

} // namespace m2424
