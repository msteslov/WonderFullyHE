#include "m2424/bootstrap_pipeline.hpp"

#include <algorithm>
#include <chrono>
#include <stdexcept>
#include <utility>

namespace m2424 {

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

} // namespace m2424
