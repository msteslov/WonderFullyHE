#pragma once

#include "m2424/bootstrap_prototype.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace m2424 {

// A declarative description of one interchangeable bootstrap prototype variant.
struct BootstrapExperimentConfig {
    std::string name{"unnamed"};
    std::size_t slots{4};
    double tolerance{1e-3};
    double normalization_factor{1.0};
    BootstrapCircuitOrder circuit_order{BootstrapCircuitOrder::SlotsToCoeffsFirst};
    BootstrapTransformBackend transform_backend{BootstrapTransformBackend::FftLike};
    EvalModDegree evalmod_degree{EvalModDegree::P3};
    EvalModEvaluationPolicy evalmod_policy{EvalModEvaluationPolicy::Polynomial};
    BootstrapNormalizationMode normalization_mode{BootstrapNormalizationMode::PlainMultiplyRescale};
    BootstrapDenormalizationPosition denormalization_position{
        BootstrapDenormalizationPosition::AfterSlotToCoeff};
    BootstrapPeriodMode period_mode{BootstrapPeriodMode::TotalCoeffModulus};
    double manual_period_log2{};
    double plain_scale_log2{40.0};
    double output_correction_factor{1.0};
    bool post_refresh_mod_raise_enabled{};
    std::size_t stc_first_target_chain_index{2};
    double stc_first_period_offset_log2{3.0};
};

enum class BootstrapExperimentOutcome {
    Passed,
    Blocked,
    Failed
};

struct BootstrapExperimentResult {
    BootstrapExperimentConfig config;
    std::vector<int> rotation_steps;
    BootstrapExperimentOutcome outcome{BootstrapExperimentOutcome::Failed};
    std::string blocker;
    BootstrapPrototypeReport report;
};

const char* to_string(BootstrapExperimentOutcome outcome) noexcept;

// Validates only configuration-level constraints and returns the exact Galois
// rotation steps needed by the selected transform implementation.
std::vector<int> bootstrap_experiment_rotation_steps(const BootstrapExperimentConfig& config);

// Runs one configured variant. The caller owns adapter setup and must create the
// returned rotation keys before invoking this function. expected may be null for
// an execution-only experiment; otherwise every stage is checked against it.
BootstrapExperimentResult run_bootstrap_experiment(SealAdapter& adapter,
                                                    const Cipher& input,
                                                    const ComplexVector* expected,
                                                    const BootstrapExperimentConfig& config);

std::vector<BootstrapExperimentResult> run_bootstrap_experiments(
    SealAdapter& adapter,
    const Cipher& input,
    const ComplexVector* expected,
    const std::vector<BootstrapExperimentConfig>& configs);

} // namespace m2424
