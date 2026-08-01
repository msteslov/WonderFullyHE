#include "m2424/bootstrap_candidates.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace m2424 {
namespace {

constexpr BootstrapErrorBudget kErrorBudget{
    1e-10,
    1e-10,
    2e-10,
    2e-10,
    2e-10,
    2e-10
};

BootstrapCandidate makeCandidate(std::string id,
                                 std::string description,
                                 std::size_t activeSlots,
                                 std::size_t polyModulusDegree,
                                 double targetScaleLog2,
                                 BootstrapResourceLimits resources) {
    return BootstrapCandidate{
        std::move(id),
        std::move(description),
        activeSlots,
        polyModulusDegree,
        targetScaleLog2,
        1.0,
        128,
        resources,
        kErrorBudget
    };
}

bool isFiniteNonNegative(const std::optional<double>& value) {
    return value.has_value() && std::isfinite(*value) && *value >= 0.0;
}

double totalBudget(const BootstrapErrorBudget& budget) {
    return budget.encode + budget.modRaise + budget.coeffToSlot
        + budget.evalMod + budget.slotToCoeff + budget.reserve;
}

} // namespace

const std::vector<BootstrapCandidate>& bootstrapCandidates() {
    static const std::vector<BootstrapCandidate> candidates{
        makeCandidate("validation_4_s55", "Математическая валидация компонентов на 4 слотах.",
                      4, 8192, 55.0, {1'000.0, 1'024, 256}),
        makeCandidate("validation_16_s55", "Валидация роста ошибок на 16 слотах.",
                      16, 8192, 55.0, {1'000.0, 1'024, 256}),
        makeCandidate("validation_64_s55", "Валидация малой пакетной обработки на 64 слотах.",
                      64, 8192, 55.0, {2'000.0, 1'024, 256}),
        makeCandidate("memory_1024_s55", "Ограниченный памятью вариант для 1024 слотов.",
                      1024, 8192, 55.0, {60'000.0, 4'096, 512}),
        makeCandidate("throughput_1024_s55", "Быстрый вариант для 1024 слотов.",
                      1024, 8192, 55.0, {10'000.0, 8'192, 2'048}),
        makeCandidate("memory_4096_s55", "Ограниченный памятью вариант для 4096 слотов.",
                      4096, 8192, 55.0, {120'000.0, 8'192, 1'024}),
        makeCandidate("throughput_4096_s55", "Базовый вариант пакетной обработки на 4096 слотах.",
                      4096, 8192, 55.0, {15'000.0, 12'288, 3'072}),
        makeCandidate("balanced_8192_s55", "Сбалансированный вариант для 8192 слотов.",
                      8192, 16384, 55.0, {30'000.0, 16'384, 4'096}),
        makeCandidate("precision_8192_s59", "Вариант с повышенным scale для 8192 слотов.",
                      8192, 16384, 59.0, {60'000.0, 16'384, 4'096}),
        makeCandidate("memory_8192_s55", "Экономный по ключам вариант для 8192 слотов.",
                      8192, 16384, 55.0, {180'000.0, 12'288, 1'536}),
        makeCandidate("throughput_16384_s55", "Высокая пакетная производительность на 16384 слотах.",
                      16384, 32768, 55.0, {60'000.0, 24'576, 6'144}),
        makeCandidate("precision_16384_s59", "Максимальный запас по scale для 16384 слотов.",
                      16384, 32768, 59.0, {120'000.0, 32'768, 8'192})
    };
    return candidates;
}

const BootstrapCandidate& bootstrapCandidateById(const std::string& id) {
    const auto& candidates = bootstrapCandidates();
    const auto it = std::find_if(candidates.begin(), candidates.end(), [&](const BootstrapCandidate& candidate) {
        return candidate.id == id;
    });
    if (it == candidates.end()) {
        throw std::invalid_argument("unknown bootstrap candidate: " + id);
    }
    return *it;
}

bool isBootstrapCandidateValid(const BootstrapCandidate& candidate) {
    if (candidate.activeSlots == 0 || candidate.polyModulusDegree < 2 * candidate.activeSlots
        || candidate.targetScaleLog2 <= 0.0 || !std::isfinite(candidate.targetScaleLog2)
        || candidate.inputAbsBound <= 0.0 || !std::isfinite(candidate.inputAbsBound)
        || candidate.minimumSecurityBits < 128 || candidate.resources.maxLatencyMs <= 0.0
        || !std::isfinite(candidate.resources.maxLatencyMs) || candidate.resources.maxRamMiB == 0
        || candidate.resources.maxEvaluationKeyMiB == 0) {
        return false;
    }
    return totalBudget(candidate.errorBudget) <= kTargetAbsoluteError;
}

BootstrapFeasibilityForecast forecastBootstrapFeasibility(const BootstrapCandidate& candidate,
                                                          const BootstrapComponentBounds& bounds) {
    if (!isBootstrapCandidateValid(candidate)) {
        return {false, false, 0.0, "invalid_candidate"};
    }
    if (!isFiniteNonNegative(bounds.encode) || !isFiniteNonNegative(bounds.modRaise)
        || !isFiniteNonNegative(bounds.coeffToSlot) || !isFiniteNonNegative(bounds.evalMod)
        || !isFiniteNonNegative(bounds.slotToCoeff)) {
        return {false, false, 0.0, "missing_confirmed_component_bound"};
    }

    const double predicted = *bounds.encode + *bounds.modRaise + *bounds.coeffToSlot
        + *bounds.evalMod + *bounds.slotToCoeff + candidate.errorBudget.reserve;
    const bool withinComponentBudgets = *bounds.encode <= candidate.errorBudget.encode
        && *bounds.modRaise <= candidate.errorBudget.modRaise
        && *bounds.coeffToSlot <= candidate.errorBudget.coeffToSlot
        && *bounds.evalMod <= candidate.errorBudget.evalMod
        && *bounds.slotToCoeff <= candidate.errorBudget.slotToCoeff;
    return {true, withinComponentBudgets && predicted <= kTargetAbsoluteError, predicted,
            withinComponentBudgets && predicted <= kTargetAbsoluteError ? "" : "error_budget_exceeded"};
}

double propagatedBootstrapError(const BootstrapPropagationBounds& bounds) {
    const double values[]{bounds.coeffToSlotNormalized, bounds.approximation,
                          bounds.polynomialArithmetic, bounds.scaleMismatch,
                          bounds.periodMismatch, bounds.evalModDerivativeMax,
                          bounds.slotToCoeffOperatorNorm, bounds.slotToCoeffAdditive,
                          bounds.finalAdditive};
    for (double value : values) {
        if (!std::isfinite(value) || value < 0.0) {
            throw std::invalid_argument("invalid bootstrap propagation bound");
        }
    }
    return bounds.slotToCoeffOperatorNorm
        * (bounds.evalModDerivativeMax * bounds.coeffToSlotNormalized
           + bounds.approximation + bounds.polynomialArithmetic
           + bounds.scaleMismatch + bounds.periodMismatch)
        + bounds.slotToCoeffAdditive + bounds.finalAdditive;
}

} // namespace m2424
