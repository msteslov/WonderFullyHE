#include "m2424/coeff_to_slot_contract.hpp"

#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <vector>

namespace {

m2424::CkksProfile calibrationProfile(double scaleLog2) {
    return {16384, {60, static_cast<int>(scaleLog2), static_cast<int>(scaleLog2), 60},
            std::exp2(scaleLog2), 8192};
}

bool rejectsZeroRotation(const m2424::SealAdapter& adapter,
                         const m2424::RaisedCipher& input,
                         const m2424::CoeffToSlotContract& contract) {
    try {
        (void)m2424::preflightCoeffToSlot(adapter, input, contract, {1, {0}});
    } catch (const std::invalid_argument&) {
        return true;
    }
    return false;
}

} // namespace

int main() {
    const auto& candidate = m2424::bootstrapCandidateById("precision_8192_s59");
    const auto contract = m2424::makeCoeffToSlotContract(candidate);
    const m2424::CoeffToSlotPlanRequirements requirements{1, {1, 2, 4}, true};

    auto adapter = m2424::SealAdapter::create(calibrationProfile(59.0));
    auto keySteps = requirements.rotationSteps;
    keySteps.push_back(0);
    adapter.generateKeys(keySteps, true);
    const auto input = adapter.encrypt(adapter.encode({0.25, -0.5, 0.75, -1.0}));
    const auto lowered = adapter.rescaleToNext(
        adapter.multiplyPlain(input, adapter.encodeScalarFor(1.0, input)));
    const auto raised = adapter.modRaiseToTop(lowered);
    const auto ready = m2424::preflightCoeffToSlot(adapter, raised, contract, requirements);

    auto missingKeysAdapter = m2424::SealAdapter::create(calibrationProfile(59.0));
    missingKeysAdapter.generateKeys(false, false);
    const auto missingKeysInput = missingKeysAdapter.encrypt(missingKeysAdapter.encode({0.25, -0.5}));
    const auto missingKeysLowered = missingKeysAdapter.rescaleToNext(
        missingKeysAdapter.multiplyPlain(
            missingKeysInput, missingKeysAdapter.encodeScalarFor(1.0, missingKeysInput)));
    const auto missingKeysRaised = missingKeysAdapter.modRaiseToTop(missingKeysLowered);
    const auto missingKeys = m2424::preflightCoeffToSlot(
        missingKeysAdapter, missingKeysRaised, contract, requirements);

    auto missingConjugationAdapter = m2424::SealAdapter::create(calibrationProfile(59.0));
    missingConjugationAdapter.generateKeys(requirements.rotationSteps, true);
    const auto missingConjugationInput =
        missingConjugationAdapter.encrypt(missingConjugationAdapter.encode({0.25, -0.5}));
    const auto missingConjugationLowered = missingConjugationAdapter.rescaleToNext(
        missingConjugationAdapter.multiplyPlain(
            missingConjugationInput,
            missingConjugationAdapter.encodeScalarFor(1.0, missingConjugationInput)));
    const auto missingConjugationRaised =
        missingConjugationAdapter.modRaiseToTop(missingConjugationLowered);
    const auto missingConjugation = m2424::preflightCoeffToSlot(
        missingConjugationAdapter, missingConjugationRaised, contract, requirements);

    auto scaleMismatchAdapter = m2424::SealAdapter::create(calibrationProfile(55.0));
    scaleMismatchAdapter.generateKeys(keySteps, true);
    const auto scaleMismatchInput = scaleMismatchAdapter.encrypt(scaleMismatchAdapter.encode({0.25, -0.5}));
    const auto scaleMismatchLowered = scaleMismatchAdapter.rescaleToNext(
        scaleMismatchAdapter.multiplyPlain(
            scaleMismatchInput, scaleMismatchAdapter.encodeScalarFor(1.0, scaleMismatchInput)));
    const auto scaleMismatchRaised = scaleMismatchAdapter.modRaiseToTop(scaleMismatchLowered);
    const auto scaleMismatch = m2424::preflightCoeffToSlot(
        scaleMismatchAdapter, scaleMismatchRaised, contract, requirements);

    const bool ok = contract.candidateId == candidate.id
        && contract.maxAbsError == candidate.errorBudget.coeffToSlot
        && ready.ready && ready.blocker.empty()
        && !missingKeys.ready && missingKeys.blocker == "missing_rotation_keys"
        && !missingConjugation.ready && missingConjugation.blocker == "missing_conjugation_key"
        && !scaleMismatch.ready && scaleMismatch.blocker == "input_scale_mismatch"
        && m2424::isCoeffToSlotPlanRequirementsValid(requirements)
        && rejectsZeroRotation(adapter, raised, contract);

    std::printf("[test_coeff_to_slot_contract] budget=%.1e %s\n",
                contract.maxAbsError, ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
