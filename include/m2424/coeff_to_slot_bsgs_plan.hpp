#pragma once

#include "m2424/canonical_embedding_reference.hpp"
#include "m2424/coeff_to_slot_contract.hpp"

#include <cstddef>
#include <vector>

namespace m2424 {

/**
 * BSGS-кандидат CoeffToSlot для forward DFT-матрицы малого размера.
 *
 * План раскладывает диагонали на baby и giant steps. Baby rotations вычисляются один раз,
 * затем переиспользуются во всех группах. План расходует один CKKS-уровень независимо от
 * числа giant-групп, поскольку plaintext-умножения выполняются от ciphertext одного уровня.
 *
 * До реализации ModUp зашифрованная проверка использует синтетический slot fixture.
 * Следовательно, план подтверждает исполнение DFT-матрицы, но не весь bootstrap CoeffToSlot.
 */
class CoeffToSlotBsgsPlan {
public:
    /// Строит BSGS-план размера transformSlots с заданной длиной baby step.
    CoeffToSlotBsgsPlan(std::size_t transformSlots, std::size_t physicalSlotCount, std::size_t babyStep);

    /// Возвращает расход уровней плана.
    std::size_t requiredLevels() const noexcept;
    /// Возвращает все rotation steps для baby и giant rotations.
    std::vector<int> rotationSteps() const;
    /// Вычисляет тот же forward DFT в plaintext для проверки факторизации.
    ComplexVector applyPlain(const ComplexVector& input) const;
    /// Кодирует неизменяемые BSGS-диагонали для уровня и scale входного ciphertext.
    void prepare(SealAdapter& adapter, const Cipher& input) const;
    /// Возвращает размер подготовленных plaintext-диагоналей в сериализованном виде.
    std::size_t preparedPlaintextBytes(SealAdapter& adapter) const;
    /// Применяет BSGS-план к ciphertext.
    Cipher apply(SealAdapter& adapter, const Cipher& input) const;
    /// Формирует требования BSGS-плана для preflightCoeffToSlot.
    CoeffToSlotPlanRequirements requirements() const;

private:
    struct WeightedBabyTerm {
        int babyRotation{};
        ComplexVector diagonal;
    };
    struct GiantGroup {
        int giantRotation{};
        std::vector<WeightedBabyTerm> terms;
    };

    std::size_t transformSlots_{};
    std::size_t physicalSlotCount_{};
    std::size_t babyStep_{};
    std::vector<GiantGroup> groups_;
    mutable const SealAdapter* preparedAdapter_{};
    mutable std::size_t preparedChainIndex_{};
    mutable double preparedInputScale_{};
    mutable std::vector<std::vector<Plain>> preparedDiagonals_;
};

} // namespace m2424
