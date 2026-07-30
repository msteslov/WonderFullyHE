#include "m2424/coeff_to_slot_contract.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace m2424 {
namespace {

bool isFinitePositive(double value) {
    return std::isfinite(value) && value > 0.0;
}

bool containsZeroRotation(const std::vector<int>& rotationSteps) {
    return std::find(rotationSteps.begin(), rotationSteps.end(), 0) != rotationSteps.end();
}

} // namespace

CoeffToSlotContract makeCoeffToSlotContract(const BootstrapCandidate& candidate) {
    if (!isBootstrapCandidateValid(candidate)) {
        throw std::invalid_argument("cannot create CoeffToSlot contract from invalid bootstrap candidate");
    }
    return {
        candidate.id,
        candidate.activeSlots,
        candidate.polyModulusDegree,
        candidate.targetScaleLog2,
        candidate.targetScaleLog2,
        0.25,
        candidate.errorBudget.coeffToSlot
    };
}

bool isCoeffToSlotPlanRequirementsValid(const CoeffToSlotPlanRequirements& requirements) {
    return !containsZeroRotation(requirements.rotationSteps);
}

CoeffToSlotPreflight preflightCoeffToSlot(const SealAdapter& adapter,
                                          const RaisedCipher& input,
                                          const CoeffToSlotContract& contract,
                                          const CoeffToSlotPlanRequirements& requirements) {
    if (!isFinitePositive(contract.inputScaleLog2) || !isFinitePositive(contract.outputScaleLog2)
        || !isFinitePositive(contract.inputScaleToleranceLog2) || !isFinitePositive(contract.maxAbsError)
        || contract.activeSlots == 0 || contract.polyModulusDegree < 2 * contract.activeSlots) {
        throw std::invalid_argument("invalid CoeffToSlot contract");
    }
    if (!isCoeffToSlotPlanRequirementsValid(requirements)) {
        throw std::invalid_argument("CoeffToSlot plan must not request rotation step zero");
    }

    const auto info = adapter.info(input);
    const bool physicalSlotsAvailable = adapter.slotCount() >= contract.activeSlots;
    const bool inputScaleMatches = std::abs(std::log2(info.scale) - contract.inputScaleLog2)
        <= contract.inputScaleToleranceLog2;
    const bool levelsAvailable = info.chainIndex >= requirements.minRemainingLevels;
    const bool rotationKeysAvailable = adapter.hasRotationKeys(requirements.rotationSteps);
    const bool conjugationKeyAvailable = !requirements.requiresConjugation || adapter.hasConjugationKey();
    const bool ready = physicalSlotsAvailable && inputScaleMatches && levelsAvailable
        && rotationKeysAvailable && conjugationKeyAvailable;

    std::string blocker;
    if (!physicalSlotsAvailable) {
        blocker = "insufficient_physical_slots";
    } else if (!inputScaleMatches) {
        blocker = "input_scale_mismatch";
    } else if (!levelsAvailable) {
        blocker = "insufficient_remaining_levels";
    } else if (!rotationKeysAvailable) {
        blocker = "missing_rotation_keys";
    } else if (!conjugationKeyAvailable) {
        blocker = "missing_conjugation_key";
    }
    return {physicalSlotsAvailable, inputScaleMatches, levelsAvailable, rotationKeysAvailable,
            conjugationKeyAvailable, ready, blocker};
}

} // namespace m2424
