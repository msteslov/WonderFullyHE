#pragma once

#include "m2424/bootstrap_candidates.hpp"
#include "m2424/seal_adapter.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace m2424 {

/**
 * Неподвижные требования CoeffToSlot для одного кандидата bootstrap.
 *
 * Входной ciphertext представляет коэффициенты полинома, нормализованные делением на
 * 2^inputScaleLog2. Будущий ModUp обязан сохранить именно эту семантику. Преобразование
 * возвращает комплексные слоты в порядке Microsoft SEAL с целевым outputScaleLog2.
 */
struct CoeffToSlotContract {
    std::string candidateId;
    std::size_t activeSlots{};
    std::size_t polyModulusDegree{};
    double inputScaleLog2{};
    double outputScaleLog2{};
    double inputScaleToleranceLog2{};
    double maxAbsError{};
};

/// Требования конкретного FFT или BSGS-плана к ключам и доступным уровням ciphertext.
struct CoeffToSlotPlanRequirements {
    std::size_t minRemainingLevels{};
    std::vector<int> rotationSteps;
};

/// Результат проверки совместимости ciphertext с контрактом и выбранным планом CoeffToSlot.
struct CoeffToSlotPreflight {
    bool physicalSlotsAvailable{};
    bool inputScaleMatches{};
    bool levelsAvailable{};
    bool rotationKeysAvailable{};
    bool ready{};
    std::string blocker;
};

/// Создаёт неизменяемый контракт CoeffToSlot из зарегистрированного кандидата bootstrap.
CoeffToSlotContract makeCoeffToSlotContract(const BootstrapCandidate& candidate);

/// Проверяет корректность требований конкретного плана до обращения к ciphertext.
bool isCoeffToSlotPlanRequirementsValid(const CoeffToSlotPlanRequirements& requirements);

/// Проверяет, что ciphertext можно передать выбранному CoeffToSlot-плану без неявных преобразований scale.
CoeffToSlotPreflight preflightCoeffToSlot(const SealAdapter& adapter,
                                          const Cipher& input,
                                          const CoeffToSlotContract& contract,
                                          const CoeffToSlotPlanRequirements& requirements);

} // namespace m2424
