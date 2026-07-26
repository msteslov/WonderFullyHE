#include "m2424/bootstrap_experiment.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace m2424 {
namespace {

bool is_power_of_two(std::size_t value) {
    return value != 0 && (value & (value - 1)) == 0;
}

void validate_config(const BootstrapExperimentConfig& config) {
    if (config.name.empty()) {
        throw std::invalid_argument("bootstrap experiment name must not be empty");
    }
    if (!is_power_of_two(config.slots)) {
        throw std::invalid_argument("bootstrap experiment slots must be a non-zero power of two");
    }
    if (!std::isfinite(config.tolerance) || config.tolerance < 0.0) {
        throw std::invalid_argument("bootstrap experiment tolerance must be a non-negative finite value");
    }
    if (!std::isfinite(config.normalization_factor) || config.normalization_factor <= 0.0) {
        throw std::invalid_argument("bootstrap experiment normalization factor must be positive and finite");
    }
    if (!std::isfinite(config.plain_scale_log2) || config.plain_scale_log2 <= 0.0) {
        throw std::invalid_argument("bootstrap experiment plain scale log2 must be positive and finite");
    }
    if (!std::isfinite(config.manual_period_log2) || config.manual_period_log2 < 0.0) {
        throw std::invalid_argument("bootstrap experiment manual period log2 must be non-negative and finite");
    }
    if (!std::isfinite(config.stc_first_period_offset_log2)) {
        throw std::invalid_argument("bootstrap experiment STC-first period offset must be finite");
    }
    if (config.circuit_order == BootstrapCircuitOrder::SlotsToCoeffsFirst
        && config.evalmod_degree != EvalModDegree::P3
        && config.evalmod_degree != EvalModDegree::P3DoubleAngle) {
        throw std::invalid_argument("SlotsToCoeffsFirst experiment supports only P3/P3DoubleAngle");
    }
}

void append_steps(std::vector<int>& destination, const std::vector<int>& source) {
    destination.insert(destination.end(), source.begin(), source.end());
}

void normalize_steps(std::vector<int>& steps) {
    std::sort(steps.begin(), steps.end());
    steps.erase(std::unique(steps.begin(), steps.end()), steps.end());
}

BootstrapPrototype make_prototype(SealAdapter& adapter, const BootstrapExperimentConfig& config) {
    BootstrapPrototype prototype(adapter, config.slots, config.tolerance, config.normalization_factor);
    prototype.set_circuit_order(config.circuit_order);
    prototype.set_transform_backend(config.transform_backend);
    prototype.set_evalmod_degree(config.evalmod_degree);
    prototype.set_evalmod_policy(config.evalmod_policy);
    prototype.set_normalization_mode(config.normalization_mode);
    prototype.set_denormalization_position(config.denormalization_position);
    prototype.set_period_mode(config.period_mode);
    prototype.set_manual_period_log2(config.manual_period_log2);
    prototype.set_plain_scale_log2(config.plain_scale_log2);
    prototype.set_output_correction_factor(config.output_correction_factor);
    prototype.set_post_refresh_mod_raise_enabled(config.post_refresh_mod_raise_enabled);
    prototype.set_stc_first_target_chain_index(config.stc_first_target_chain_index);
    prototype.set_stc_first_period_offset_log2(config.stc_first_period_offset_log2);
    return prototype;
}

BootstrapExperimentOutcome classify(const BootstrapPrototypeReport& report, std::string& blocker) {
    if (report.stages.empty()) {
        blocker = "empty_report";
        return BootstrapExperimentOutcome::Failed;
    }
    for (const auto& stage : report.stages) {
        if (stage.status == "BLOCKED") {
            blocker = stage.name;
            return BootstrapExperimentOutcome::Blocked;
        }
        if (stage.status == "FAIL") {
            blocker = stage.name;
            return BootstrapExperimentOutcome::Failed;
        }
    }
    if (report.checked && !report.preserve_value_criterion) {
        blocker = "preserve_value_criterion";
        return BootstrapExperimentOutcome::Failed;
    }
    if (!report.restore_level_criterion) {
        blocker = "restore_level_criterion";
        return BootstrapExperimentOutcome::Failed;
    }
    blocker = "none";
    return BootstrapExperimentOutcome::Passed;
}

} // namespace

const char* to_string(BootstrapExperimentOutcome outcome) noexcept {
    switch (outcome) {
    case BootstrapExperimentOutcome::Passed:
        return "passed";
    case BootstrapExperimentOutcome::Blocked:
        return "blocked";
    case BootstrapExperimentOutcome::Failed:
        return "failed";
    }
    return "unknown";
}

std::vector<int> bootstrap_experiment_rotation_steps(const BootstrapExperimentConfig& config) {
    validate_config(config);
    if (config.transform_backend == BootstrapTransformBackend::DenseDiagonal) {
        return BootstrapPrototype::required_rotation_steps(config.slots);
    }

    std::vector<int> steps;
    if (config.transform_backend == BootstrapTransformBackend::SmallSlots4Butterfly) {
        append_steps(steps, make_small_slots4_butterfly_stc_plan(config.plain_scale_log2).rotation_steps());
        append_steps(steps, make_small_slots4_butterfly_cts_plan(config.plain_scale_log2).rotation_steps());
    } else {
        append_steps(steps, make_bootstrap_dft_plan(
            config.slots, BootstrapDftType::HomomorphicEncode, config.plain_scale_log2).rotation_steps());
        append_steps(steps, make_bootstrap_dft_plan(
            config.slots, BootstrapDftType::HomomorphicDecode, config.plain_scale_log2).rotation_steps());
    }
    normalize_steps(steps);
    return steps;
}

BootstrapExperimentResult run_bootstrap_experiment(SealAdapter& adapter,
                                                    const Cipher& input,
                                                    const ComplexVector* expected,
                                                    const BootstrapExperimentConfig& config) {
    validate_config(config);
    BootstrapExperimentResult result;
    result.config = config;
    result.rotation_steps = bootstrap_experiment_rotation_steps(config);

    try {
        auto prototype = make_prototype(adapter, config);
        result.report = expected == nullptr
            ? prototype.refresh_cipher_fast(input)
            : prototype.refresh_cipher_checked(input, *expected);
        result.outcome = classify(result.report, result.blocker);
    } catch (const std::exception& error) {
        result.outcome = BootstrapExperimentOutcome::Failed;
        result.blocker = error.what();
    }
    return result;
}

std::vector<BootstrapExperimentResult> run_bootstrap_experiments(
    SealAdapter& adapter,
    const Cipher& input,
    const ComplexVector* expected,
    const std::vector<BootstrapExperimentConfig>& configs) {
    if (configs.empty()) {
        throw std::invalid_argument("bootstrap experiment suite must not be empty");
    }
    std::vector<BootstrapExperimentResult> results;
    results.reserve(configs.size());
    for (const auto& config : configs) {
        results.push_back(run_bootstrap_experiment(adapter, input, expected, config));
    }
    return results;
}

} // namespace m2424
