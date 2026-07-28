#include "m2424/bootstrap_pipeline.hpp"
#include "m2424/accuracy.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <utility>

namespace m2424 {

BootstrapVectorReferenceOracle::BootstrapVectorReferenceOracle(std::vector<double> expected)
    : expected_(std::move(expected)) {
    if (expected_.empty()) {
        throw std::invalid_argument("bootstrap reference vector must not be empty");
    }
}

BootstrapReferenceReport BootstrapVectorReferenceOracle::validate(SealAdapter& adapter,
                                                                   const Cipher& ciphertext,
                                                                   double targetError) const {
    if (targetError < 0.0 || !std::isfinite(targetError)) {
        throw std::invalid_argument("bootstrap reference tolerance must be finite and non-negative");
    }
    std::vector<double> actual = adapter.decode(adapter.decrypt(ciphertext));
    if (actual.size() < expected_.size()) {
        return {true, false, 0.0, "decoded ciphertext has fewer slots than reference"};
    }
    actual.resize(expected_.size());
    const AccuracyReport accuracy = compare(expected_, actual, targetError);
    return {true, accuracy.ok, accuracy.max_abs_error, accuracy.ok ? "ok" : "reference mismatch"};
}

BootstrapPipeline::BootstrapPipeline(std::vector<std::unique_ptr<BootstrapStage>> stages)
    : stages_(std::move(stages)) {
    if (stages_.empty()) {
        throw std::invalid_argument("bootstrap pipeline must contain at least one stage");
    }
    if (std::any_of(stages_.begin(), stages_.end(), [](const auto& stage) { return !stage; })) {
        throw std::invalid_argument("bootstrap pipeline must not contain null stages");
    }
}

BootstrapKeyRequirements BootstrapPipeline::keyRequirements() const {
    BootstrapKeyRequirements result;
    for (const auto& stage : stages_) {
        const auto requirements = stage->keyRequirements();
        result.requiresRelin = result.requiresRelin || requirements.requiresRelin;
        if (std::any_of(requirements.rotationSteps.begin(), requirements.rotationSteps.end(),
                        [](int step) { return step == 0; })) {
            throw std::invalid_argument("bootstrap stage must not request zero rotation");
        }
        result.rotationSteps.insert(result.rotationSteps.end(),
                                    requirements.rotationSteps.begin(),
                                    requirements.rotationSteps.end());
    }
    std::sort(result.rotationSteps.begin(), result.rotationSteps.end());
    result.rotationSteps.erase(std::unique(result.rotationSteps.begin(), result.rotationSteps.end()),
                               result.rotationSteps.end());
    return result;
}

std::size_t BootstrapPipeline::estimatedLevelCost() const noexcept {
    std::size_t result = 0;
    for (const auto& stage : stages_) {
        const std::size_t cost = stage->estimatedLevelCost();
        if (cost > std::numeric_limits<std::size_t>::max() - result) {
            return std::numeric_limits<std::size_t>::max();
        }
        result += cost;
    }
    return result;
}

BootstrapPreparationPlan BootstrapPipeline::planPreparation(const SealAdapter& adapter,
                                                             const Cipher& input,
                                                             const BootstrapExecutionConfig& config) const {
    BootstrapPreparationPlan plan;
    plan.keyRequirements = keyRequirements();
    plan.estimatedLevelsConsumed = estimatedLevelCost();

    const CipherInfo inputInfo = adapter.info(input);
    plan.keysAvailable = (!plan.keyRequirements.requiresRelin || adapter.hasRelinKeys())
        && adapter.hasRotationKeys(plan.keyRequirements.rotationSteps);
    plan.levelsAvailable = inputInfo.chainIndex >= config.minRemainingLevels
        && inputInfo.chainIndex - config.minRemainingLevels >= plan.estimatedLevelsConsumed;
    plan.ready = plan.keysAvailable && plan.levelsAvailable;
    if (!plan.keysAvailable) {
        plan.blocker = "required evaluation keys are not available";
    } else if (!plan.levelsAvailable) {
        plan.blocker = "insufficient modulus-chain levels for declared stage costs";
    }
    return plan;
}

BootstrapPipelineResult BootstrapPipeline::run(SealAdapter& adapter,
                                                const Cipher& input,
                                                const BootstrapMeasurementConfig& measurements) const {
    Cipher current = input;
    BootstrapPipelineResult result;
    result.stages.reserve(stages_.size());

    for (const auto& stage : stages_) {
        const CipherInfo before = adapter.info(current);
        const auto startedAt = std::chrono::steady_clock::now();
        BootstrapStageResult stageResult = stage->execute(adapter, current);
        const auto finishedAt = std::chrono::steady_clock::now();

        BootstrapStageReport report;
        report.kind = stage->kind();
        report.name = std::string(stage->name());
        report.input = before;
        report.output = adapter.info(stageResult.ciphertext);
        report.chainIndexDelta = static_cast<std::int64_t>(report.output.chainIndex)
            - static_cast<std::int64_t>(report.input.chainIndex);
        if (measurements.measureDuration) {
            report.durationMs = std::chrono::duration<double, std::milli>(finishedAt - startedAt).count();
        }
        result.stages.push_back(std::move(report));
        current = std::move(stageResult.ciphertext);
    }

    result.ciphertext = std::move(current);
    return result;
}

BootstrapPipelineResult BootstrapPipeline::run(SealAdapter& adapter,
                                                const Cipher& input,
                                                const BootstrapExecutionConfig& config,
                                                const BootstrapReferenceOracle* oracle) const {
    const BootstrapPreparationPlan plan = planPreparation(adapter, input, config);
    if (!plan.ready) {
        throw std::runtime_error("bootstrap pipeline is not ready: " + plan.blocker);
    }
    BootstrapPipelineResult result = run(adapter, input, config.measurements);
    if (oracle) {
        result.reference = oracle->validate(adapter, result.ciphertext, config.targetError);
    }
    return result;
}

void writeBootstrapCsv(std::ostream& output, const BootstrapPipelineResult& result) {
    output << "stage,kind,input_chain_index,output_chain_index,chain_index_delta,input_scale,output_scale,duration_ms\n";
    output << std::setprecision(17);
    for (const auto& stage : result.stages) {
        output << stage.name << ',' << static_cast<int>(stage.kind) << ','
               << stage.input.chainIndex << ',' << stage.output.chainIndex << ','
               << stage.chainIndexDelta << ',' << stage.input.scale << ','
               << stage.output.scale << ',' << stage.durationMs << '\n';
    }
}

} // namespace m2424
