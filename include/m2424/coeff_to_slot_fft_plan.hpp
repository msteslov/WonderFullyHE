#pragma once

#include "m2424/coeff_to_slot_contract.hpp"
#include "m2424/homomorphic_linear_transform.hpp"

#include <cstddef>
#include <vector>

namespace m2424 {

/**
 * Полный radix-2 DIT FFT-план для кандидата CoeffToSlot.
 *
 * План включает входную bit-reversal перестановку и все butterfly-слои. Он предназначен
 * для малых размеров 4, 8 и 16 как первый полностью проверяемый encrypted-кандидат.
 * Каждая стадия расходует один CKKS-уровень через HomomorphicLinearTransform. До появления
 * ModUp план проверяется на synthetic slot fixture и не считается полным bootstrap CoeffToSlot.
 */
class CoeffToSlotFftPlan {
public:
    /// Строит FFT-план логического размера transformSlots в физическом CKKS-векторе.
    CoeffToSlotFftPlan(std::size_t transformSlots, std::size_t physicalSlotCount);

    /// Возвращает число уровней, необходимое для входной перестановки и butterfly-слоёв.
    std::size_t requiredLevels() const noexcept;
    /// Возвращает набор Galois rotation steps всего FFT-плана.
    std::vector<int> rotationSteps() const;
    /// Преобразует входные plaintext-слоты в forward DFT для проверки математики плана.
    ComplexVector applyPlain(const ComplexVector& input) const;
    /// Применяет все стадии плана к ciphertext; вход должен пройти CoeffToSlot preflight.
    Cipher apply(SealAdapter& adapter, const Cipher& input) const;
    /// Формирует требования плана для preflightCoeffToSlot.
    CoeffToSlotPlanRequirements requirements() const;

private:
    std::size_t transformSlots_{};
    std::vector<HomomorphicLinearTransform> layers_;
};

} // namespace m2424
