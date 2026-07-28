#pragma once

#include "m2424/bootstrap_dft.hpp"
#include "m2424/seal_adapter.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace m2424 {

struct BootstrapStageErrorBudget {
    double target_total_error{};
    std::size_t cycles{};
    double per_cycle_budget{};
    double dft_roundtrip_budget{};
    double evalmod_budget{};
    double modraise_budget{};
};

struct BootstrapErrorRecurrence {
    double target_total_error{};
    std::size_t cycles{};
    double amplification{};
    double per_cycle_budget{};
    double dft_roundtrip_budget{};
    double rotation_budget{};
};

struct CalibratedRotationNoiseModel {
    double measured_scale_log2{};
    double measured_rotation_error{};
    double target_rotation_error{};
    double safety_bits{};
};

struct EvalModSmallSignalDecision {
    double max_abs_input{};
    double evalmod_budget{};
    double cubic_coefficient_abs{};
    double cubic_bound{};
    bool linear_path_allowed{};
};

struct DftPrecisionMeasurement {
    std::string profile_name;
    std::size_t slots{};
    double ciphertext_scale_log2{};
    double transform_plain_scale_log2{};
    std::size_t rotation_count{};
    std::size_t input_chain_index{};
    std::size_t output_chain_index{};
    std::size_t consumed_levels{};
    double baseline_error{};
    double slot_to_coeff_error{};
    double roundtrip_error{};
};

struct DftPrecisionFit {
    double quantization_coefficient{};
    double noise_floor{};
    double predicted_best_error{};
    double best_plain_scale_log2{};
    bool floor_blocks_target{};
    std::string blocker;
};

struct BootstrapDftCost {
    std::size_t slots{};
    BootstrapDftType type{BootstrapDftType::HomomorphicDecode};
    std::size_t layer_count{};
    std::size_t diagonal_term_count{};
    std::size_t rotation_count{};
    std::size_t plaintext_multiplication_count{};
    std::size_t rescale_count{};
    std::vector<int> rotationSteps;
};

struct BootstrapPrecisionPlanningRequest {
    double target_total_error{};
    std::size_t cycles{};
    std::size_t slots{};
    std::vector<CkksProfile> candidate_profiles;
    std::vector<double> candidate_transform_scales;
};

struct BootstrapPrecisionPlanningResult {
    bool feasible{};
    CkksProfile selected_profile{};
    double selected_transform_plain_scale_log2{};
    BootstrapStageErrorBudget budget{};
    DftPrecisionFit dft_fit{};
    std::string blocker;
};

BootstrapStageErrorBudget make_bootstrap_error_budget(double target_total_error,
                                                      std::size_t cycles);

BootstrapErrorRecurrence make_bootstrap_error_recurrence(double target_total_error,
                                                        std::size_t cycles,
                                                        double amplification);

double required_ciphertext_scale_log2(const CalibratedRotationNoiseModel& model);

EvalModSmallSignalDecision decide_evalmod_small_signal(double max_abs_input,
                                                       double evalmod_budget);

DftPrecisionFit fit_dft_precision_floor(const std::vector<DftPrecisionMeasurement>& measurements,
                                        double target_roundtrip_error);

BootstrapDftCost estimate_bootstrap_dft_cost(std::size_t slots,
                                             BootstrapDftType type,
                                             double plain_scale_log2 = 60.0);

BootstrapPrecisionPlanningResult plan_bootstrap_precision(
    const BootstrapPrecisionPlanningRequest& request,
    const std::vector<DftPrecisionMeasurement>& calibration);

} // namespace m2424
