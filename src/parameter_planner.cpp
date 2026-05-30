#include "m2424/parameter_planner.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace m2424 {
namespace {

SecurityLevel security_level_from_bits(int bits) {
    switch (bits) {
    case 128:
        return SecurityLevel::TC128;
    case 192:
        return SecurityLevel::TC192;
    case 256:
        return SecurityLevel::TC256;
    default:
        throw std::invalid_argument("security_bits must be one of 128, 192, 256");
    }
}

bool passes_requested_security(const SecurityReport& report, SecurityLevel level) {
    switch (level) {
    case SecurityLevel::TC128:
        return report.passes_tc128;
    case SecurityLevel::TC192:
        return report.passes_tc192;
    case SecurityLevel::TC256:
        return report.passes_tc256;
    case SecurityLevel::None:
        return false;
    }
    return false;
}

int round_up_to_multiple_of_five(int value) {
    return ((value + 4) / 5) * 5;
}

std::size_t min_poly_degree_for_slots(std::size_t slots) {
    if (slots == 0) {
        throw std::invalid_argument("slots must be positive");
    }
    std::size_t degree = 4096;
    while (degree / 2 < slots) {
        degree *= 2;
        if (degree > 32768) {
            throw std::invalid_argument("requested slots exceed supported planner range");
        }
    }
    return degree;
}

CkksProfile make_profile(std::size_t poly_degree, int work_bits, std::size_t depth, std::size_t slots) {
    std::vector<int> bits;
    bits.reserve(depth + 2);
    bits.push_back(60);
    for (std::size_t i = 0; i < depth; ++i) {
        bits.push_back(work_bits);
    }
    bits.push_back(60);
    return CkksProfile{
        poly_degree,
        std::move(bits),
        std::exp2(static_cast<double>(work_bits)),
        slots
    };
}

} // namespace

const char* to_string(ParameterOptimizeFor value) noexcept {
    switch (value) {
    case ParameterOptimizeFor::Speed:
        return "speed";
    case ParameterOptimizeFor::Conservative:
        return "conservative";
    }
    return "unknown";
}

const char* to_string(PlanningOperationProfile value) noexcept {
    switch (value) {
    case PlanningOperationProfile::BasicMulDepth:
        return "basic_mul_depth";
    case PlanningOperationProfile::LinearTransform:
        return "linear_transform";
    case PlanningOperationProfile::BootstrapRefresh:
        return "bootstrap_refresh";
    }
    return "unknown";
}

int required_result_bits(double target_error) {
    if (!std::isfinite(target_error) || target_error <= 0.0 || target_error >= 1.0) {
        throw std::invalid_argument("target_error must be finite and in (0, 1)");
    }
    return static_cast<int>(std::ceil(-std::log2(target_error)));
}

int calibrated_loss_bits(PlanningOperationProfile profile) {
    switch (profile) {
    case PlanningOperationProfile::BasicMulDepth:
        return 14;
    case PlanningOperationProfile::LinearTransform:
        return 18;
    case PlanningOperationProfile::BootstrapRefresh:
        return 24;
    }
    return 14;
}

int calibrated_loss_bits(const CkksOperationBudget& budget) {
    int loss = 8;
    loss += static_cast<int>((budget.additions + 7) / 8);
    loss += static_cast<int>((budget.plaintext_additions + 7) / 8);
    loss += static_cast<int>(3 * budget.ciphertext_muls);
    loss += static_cast<int>(2 * budget.plaintext_mul_rescales);
    loss += static_cast<int>(2 * budget.rescale_to_next);
    loss += static_cast<int>((budget.mod_switches + 3) / 4);
    loss += static_cast<int>((budget.rotations + 1) / 2);
    loss += static_cast<int>(6 * budget.linear_transforms);
    loss += static_cast<int>(10 * budget.evalmod_p3);
    loss += static_cast<int>(24 * budget.bootstrap_refreshes);
    return loss;
}

std::size_t estimated_level_budget(const CkksOperationBudget& budget) {
    return budget.ciphertext_muls
        + budget.plaintext_mul_rescales
        + budget.rescale_to_next
        + budget.linear_transforms
        + 3 * budget.evalmod_p3;
}

CkksPlanningResult plan_ckks_parameters(const CkksPlanningRequest& request) {
    if (request.multiplicative_depth == 0) {
        throw std::invalid_argument("multiplicative_depth must be positive");
    }

    const auto requested_security = security_level_from_bits(request.security_bits);
    const int result_bits = required_result_bits(request.target_error);
    const int loss_bits = request.calibrated_loss_bits_override >= 0
        ? request.calibrated_loss_bits_override
        : request.use_operation_budget
            ? calibrated_loss_bits(request.operation_budget)
            : calibrated_loss_bits(request.operation_profile);
    if (loss_bits < 0) {
        throw std::invalid_argument("calibrated loss bits must be non-negative");
    }

    int work_bits = round_up_to_multiple_of_five(result_bits + loss_bits);
    if (request.optimize_for == ParameterOptimizeFor::Conservative) {
        work_bits = work_bits < 50 ? 50 : work_bits;
    }
    if (work_bits > 60) {
        throw std::invalid_argument("requested target requires work_bits > 60; increase strategy is unsupported");
    }
    const int estimated_precision_bits = work_bits - loss_bits;
    const double estimated_abs_error_bound = std::exp2(-static_cast<double>(estimated_precision_bits));
    const bool passes_target_error = estimated_abs_error_bound <= request.target_error;

    const std::size_t work_levels = request.use_operation_budget
        ? std::max(request.multiplicative_depth, estimated_level_budget(request.operation_budget))
        : request.multiplicative_depth;
    if (work_levels == 0) {
        throw std::invalid_argument("estimated work levels must be positive");
    }

    const std::size_t min_degree = min_poly_degree_for_slots(request.slots);
    for (std::size_t degree : {4096UL, 8192UL, 16384UL, 32768UL}) {
        if (degree < min_degree) {
            continue;
        }
        auto profile = make_profile(degree, work_bits, work_levels, request.slots);
        auto security = analyze_security("planned_ckks", profile);
        if (passes_requested_security(security, requested_security)) {
            return CkksPlanningResult{
                profile,
                request.target_error,
                result_bits,
                loss_bits,
                work_bits,
                work_bits,
                estimated_precision_bits,
                estimated_abs_error_bound,
                passes_target_error,
                work_levels,
                security
            };
        }
    }

    throw std::invalid_argument("no supported poly_modulus_degree satisfies requested depth, slots, and security");
}

std::string planning_result_summary(const CkksPlanningResult& result) {
    std::ostringstream out;
    out << "target_result_bits=" << result.required_result_bits
        << ",calibrated_loss_bits=" << result.calibrated_loss_bits
        << ",work_bits=" << result.selected_work_bits
        << ",scale_log2=" << result.selected_scale_log2
        << ",estimated_precision_bits=" << result.estimated_precision_bits
        << ",estimated_abs_error_bound=" << result.estimated_abs_error_bound
        << ",passes_target_error=" << (result.passes_target_error ? "true" : "false")
        << ",work_levels=" << result.selected_work_levels
        << ",poly_modulus_degree=" << result.profile.poly_modulus_degree
        << ",slots=" << result.profile.slots
        << ",effective_security=" << to_string(result.security.effective_level);
    return out.str();
}

} // namespace m2424
