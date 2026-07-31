#include "m2424/coeff_to_slot.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <numeric>
#include <vector>

namespace {

double maxHalfError(const std::vector<double>& coefficients,
                    std::size_t offset,
                    const std::vector<std::complex<double>>& actual) {
    double result = 0.0;
    for (std::size_t index = 0; index < actual.size(); ++index) {
        result = std::max(result, std::abs(actual[index] - coefficients[index + offset]));
    }
    return result;
}

} // namespace

int main() {
    constexpr std::size_t degree = 16384;
    constexpr std::size_t slots = degree / 2;
    constexpr double errorLimit = 2e-10;
    const m2424::CkksProfile profile{
        degree, {60, 60, 60, 60, 60, 60, 60}, std::exp2(59.5), slots
    };

    m2424::CoeffToSlot transform(degree);
    const auto requirements = transform.requirements();
    const auto metrics = transform.plan().metrics();
    const auto ranked = m2424::CoeffToSlotPlan::estimateFactorizations(degree, 4, 5);
    const auto stageRotations = transform.plan().stageRotationSteps();
    const auto& radices = transform.plan().factorization().radices;
    auto keySteps = requirements.rotationSteps;
    keySteps.push_back(0);

    auto adapter = m2424::SealAdapter::create(profile);
    adapter.generateKeys(keySteps, false);
    std::vector<double> values(slots);
    for (std::size_t index = 0; index < values.size(); ++index) {
        values[index] = static_cast<double>(static_cast<int>(index % 17) - 8) / 32.0;
    }
    const auto encrypted = adapter.encrypt(adapter.encode(values));
    auto lowered = encrypted;
    const std::size_t sourceDropCount = adapter.info(encrypted).chainIndex;
    for (std::size_t level = 0; level < sourceDropCount; ++level) {
        lowered = adapter.rescaleToNext(adapter.multiplyPlain(
            lowered, adapter.encodeScalarAtScaleFor(1.0, std::exp2(60.0), lowered)));
    }
    auto raised = adapter.modRaiseToTop(lowered);
    const auto raisedInfo = adapter.info(raised);
    const auto coefficients = adapter.decryptRaisedCoefficientsAtRaisedModulus(raised);
    const auto sourceCoefficients = adapter.decryptRaisedCoefficientsAtSourceModulus(raised);
    bool raisedTermObserved = false;
    for (std::size_t index = 0; index < coefficients.size(); ++index) {
        raisedTermObserved = raisedTermObserved
            || std::abs(coefficients[index] - sourceCoefficients[index]) > 1e-6;
    }

    const m2424::CoeffToSlotContract contract{
        "validation_full_s59_5", slots, degree, 59.5, 59.5, 0.25, errorLimit
    };
    const auto preflight = m2424::preflightCoeffToSlot(adapter, raised, contract, requirements);
    auto prepared = transform.prepare(adapter, raised, contract);
    const bool preparedBeforeApply =
        transform.plan().isPreparedFor(prepared, adapter, raised, contract);
    m2424::CoeffToSlot wrongDegreeTransform(degree / 2);
    bool rejectsWrongPlanDegreePrepare = false;
    try {
        (void)wrongDegreeTransform.prepare(adapter, raised, contract);
    } catch (const std::invalid_argument&) {
        rejectsWrongPlanDegreePrepare = true;
    }
    const bool rejectsWrongPlanDegreePreparedState =
        !wrongDegreeTransform.plan().isPreparedFor(prepared, adapter, raised, contract);
    bool rejectsWrongPlanDegreeApply = false;
    auto wrongDegreeRaised = adapter.modRaiseToTop(lowered);
    try {
        (void)wrongDegreeTransform.apply(
            adapter, std::move(wrongDegreeRaised), contract, prepared);
    } catch (const std::invalid_argument&) {
        rejectsWrongPlanDegreeApply = true;
    }
    auto wrongContract = contract;
    wrongContract.inputScaleLog2 = 55.0;
    bool rejectsWrongContract = false;
    try {
        (void)transform.apply(adapter, std::move(raised), wrongContract, prepared);
    } catch (const std::invalid_argument&) {
        rejectsWrongContract = true;
    }
    const auto result = transform.apply(adapter, std::move(raised), contract, prepared);
    const auto firstInfo = adapter.info(result.slotCipherFirst);
    const auto secondInfo = adapter.info(result.slotCipherSecond);
    const auto first = adapter.decodeComplex(adapter.decrypt(result.slotCipherFirst));
    const auto second = adapter.decodeComplex(adapter.decrypt(result.slotCipherSecond));
    const double firstError = maxHalfError(coefficients, 0, first);
    const double secondError = maxHalfError(coefficients, slots, second);
    double oracleMax = 0.0;
    double actualMax = 0.0;
    for (double value : coefficients) oracleMax = std::max(oracleMax, std::abs(value));
    for (const auto& value : first) actualMax = std::max(actualMax, std::abs(value));
    for (const auto& value : second) actualMax = std::max(actualMax, std::abs(value));

    const bool ok = preflight.ready
        && requirements.requiresConjugation
        && rejectsWrongContract
        && rejectsWrongPlanDegreePrepare
        && rejectsWrongPlanDegreePreparedState
        && rejectsWrongPlanDegreeApply
        && preparedBeforeApply
        && transform.plan().butterflyStageCount() == 13
        && !ranked.empty()
        && metrics.depth == requirements.minRemainingLevels
        && metrics.rescalesPerApply == 2 * metrics.depth
        && metrics.uniqueEvaluationKeys == requirements.rotationSteps.size() + 1
        && metrics.ordinaryGiantRotationsPerApply == 0
        && metrics.hoistedAutomorphismsPerApply == metrics.rotationsPerApply
        && metrics.finalModDownsPerApply == 2 * metrics.depth
        && metrics.storedComplexValues > 0
        && stageRotations.size() == transform.plan().depth()
        && radices.size() == transform.plan().depth()
        && std::accumulate(radices.begin(), radices.end(), std::size_t{0})
            == transform.plan().rawStageCount()
        && raisedTermObserved
        && adapter.hasConjugationKey()
        && firstError <= errorLimit
        && secondError <= errorLimit
        && firstInfo.chainIndex + requirements.minRemainingLevels == raisedInfo.chainIndex
        && secondInfo.chainIndex == firstInfo.chainIndex
        && std::abs(std::log2(firstInfo.scale) - 59.5) <= 0.25
        && std::abs(std::log2(secondInfo.scale) - 59.5) <= 0.25;

    auto replacement = m2424::SealAdapter::create({
        degree, {60, 60, 60, 60, 60, 59, 60}, std::exp2(59.5), slots
    });
    adapter = std::move(replacement);
    const bool staleContextRejected =
        !transform.plan().isPreparedFor(prepared, adapter, raised, contract);
    const bool finalOk = ok && staleContextRejected;
    std::printf("[test_coeff_to_slot] first=%.3e second=%.3e oracle=%.3e actual=%.3e raw=%zu levels=%zu keys=%zu qI=%s %s\n",
                firstError, secondError, oracleMax, actualMax, transform.plan().rawStageCount(),
                requirements.minRemainingLevels, requirements.rotationSteps.size() + 1,
                raisedTermObserved ? "yes" : "no", finalOk ? "PASS" : "FAIL");
    return finalOk ? 0 : 1;
}
