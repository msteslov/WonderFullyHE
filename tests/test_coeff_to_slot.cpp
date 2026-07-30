#include "m2424/coeff_to_slot.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
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
        degree, {60, 60, 60, 60}, std::exp2(59.0), slots
    };

    m2424::CoeffToSlot transform(degree);
    const auto requirements = transform.requirements();
    auto keySteps = requirements.rotationSteps;
    keySteps.push_back(0);

    auto adapter = m2424::SealAdapter::create(profile);
    adapter.generateKeys(keySteps, false);
    std::vector<double> values(slots);
    for (std::size_t index = 0; index < values.size(); ++index) {
        values[index] = static_cast<double>(static_cast<int>(index % 17) - 8) / 32.0;
    }
    const auto encrypted = adapter.encrypt(adapter.encode(values));
    const auto loweredOnce = adapter.rescaleToNext(adapter.multiplyPlain(
        encrypted, adapter.encodeScalarAtScaleFor(1.0, std::exp2(60.0), encrypted)));
    const auto lowered = adapter.rescaleToNext(adapter.multiplyPlain(
        loweredOnce, adapter.encodeScalarAtScaleFor(1.0, std::exp2(60.0), loweredOnce)));
    auto raised = adapter.modRaiseToTop(lowered);
    const auto raisedInfo = adapter.info(raised);
    const auto coefficients = adapter.decryptRaisedCoefficientsAtRaisedModulus(raised);

    const m2424::CoeffToSlotContract contract{
        "validation_full_s59", slots, degree, 59.0, 59.0, 0.25, errorLimit
    };
    const auto preflight = m2424::preflightCoeffToSlot(adapter, raised, contract, requirements);
    const auto result = transform.apply(adapter, std::move(raised));
    const auto firstInfo = adapter.info(result.slotCipherFirst);
    const auto secondInfo = adapter.info(result.slotCipherSecond);
    const auto first = adapter.decodeComplex(adapter.decrypt(result.slotCipherFirst));
    const auto second = adapter.decodeComplex(adapter.decrypt(result.slotCipherSecond));
    const double firstError = maxHalfError(coefficients, 0, first);
    const double secondError = maxHalfError(coefficients, slots, second);

    const bool ok = preflight.ready
        && requirements.requiresConjugation
        && adapter.hasConjugationKey()
        && firstError <= errorLimit
        && secondError <= errorLimit
        && firstInfo.chainIndex + requirements.minRemainingLevels == raisedInfo.chainIndex
        && secondInfo.chainIndex == firstInfo.chainIndex
        && std::abs(std::log2(firstInfo.scale) - 59.0) <= 0.25
        && std::abs(std::log2(secondInfo.scale) - 59.0) <= 0.25;
    std::printf("[test_coeff_to_slot] first=%.3e second=%.3e levels=%zu keys=%zu %s\n",
                firstError, secondError, requirements.minRemainingLevels,
                requirements.rotationSteps.size() + 1, ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
