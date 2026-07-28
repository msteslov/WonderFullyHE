#include "m2424/accuracy.hpp"
#include "m2424/bootstrap_pipeline.hpp"
#include "m2424/profiles.hpp"

#include <algorithm>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

class IdentityStage final : public m2424::BootstrapStage {
public:
    m2424::BootstrapStageKind kind() const noexcept override {
        return m2424::BootstrapStageKind::ModUp;
    }

    std::string_view name() const noexcept override {
        return "identity-reference";
    }

    m2424::BootstrapKeyRequirements keyRequirements() const override {
        return {false, {1, 2}};
    }

    m2424::BootstrapStageResult execute(m2424::SealAdapter&, const m2424::Cipher& input) const override {
        return {input};
    }
};

class RotateStage final : public m2424::BootstrapStage {
public:
    m2424::BootstrapStageKind kind() const noexcept override {
        return m2424::BootstrapStageKind::CoeffToSlot;
    }

    std::string_view name() const noexcept override {
        return "rotate-reference";
    }

    m2424::BootstrapKeyRequirements keyRequirements() const override {
        return {true, {1, 2}};
    }

    m2424::BootstrapStageResult execute(m2424::SealAdapter& adapter, const m2424::Cipher& input) const override {
        return {adapter.rotate(input, 1)};
    }
};

class InvalidRotationStage final : public m2424::BootstrapStage {
public:
    m2424::BootstrapStageKind kind() const noexcept override {
        return m2424::BootstrapStageKind::EvalMod;
    }

    std::string_view name() const noexcept override {
        return "invalid-rotation";
    }

    m2424::BootstrapKeyRequirements keyRequirements() const override {
        return {false, {0}};
    }

    m2424::BootstrapStageResult execute(m2424::SealAdapter&, const m2424::Cipher& input) const override {
        return {input};
    }
};

std::vector<double> head(const std::vector<double>& values, std::size_t size) {
    return {values.begin(), values.begin() + static_cast<std::ptrdiff_t>(size)};
}

bool rejectsZeroRotation() {
    std::vector<std::unique_ptr<m2424::BootstrapStage>> stages;
    stages.push_back(std::make_unique<InvalidRotationStage>());
    try {
        const m2424::BootstrapPipeline pipeline(std::move(stages));
        (void)pipeline.keyRequirements();
    } catch (const std::invalid_argument&) {
        return true;
    }
    return false;
}

} // namespace

int main() {
    auto adapter = m2424::SealAdapter::create(m2424::profiles::basic_ckks());
    adapter.generateKeys({1, 2}, true);

    const std::vector<double> input{1.0, 2.0, 3.0, 4.0};
    std::vector<std::unique_ptr<m2424::BootstrapStage>> stages;
    stages.push_back(std::make_unique<IdentityStage>());
    stages.push_back(std::make_unique<RotateStage>());
    const m2424::BootstrapPipeline pipeline(std::move(stages));

    const auto requirements = pipeline.keyRequirements();
    const auto encrypted = adapter.encrypt(adapter.encode(input));
    const auto result = pipeline.run(adapter, encrypted);
    const auto actual = head(adapter.decode(adapter.decrypt(result.ciphertext)), input.size());
    const std::vector<double> expected{2.0, 3.0, 4.0, 0.0};

    const bool keysOk = requirements.requiresRelin && requirements.rotationSteps == std::vector<int>({1, 2});
    const bool stagesOk = result.stages.size() == 2
        && result.stages[0].name == "identity-reference"
        && result.stages[1].name == "rotate-reference"
        && result.stages[0].chainIndexDelta == 0
        && result.stages[1].chainIndexDelta == 0;
    const bool accuracyOk = m2424::compare(expected, actual, 1e-5).ok;
    const bool ok = keysOk && stagesOk && accuracyOk && rejectsZeroRotation();
    std::printf("[test_bootstrap_pipeline] %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
