#pragma once

#include "m2424/coeff_to_slot_contract.hpp"
#include "m2424/homomorphic_linear_transform.hpp"
#include "m2424/seal_adapter.hpp"

#include <cstddef>
#include <memory>
#include <vector>

namespace m2424 {

/// Два slot ciphertext, содержащие первую и вторую половины коэффициентов RaisedCipher.
struct CoeffToSlotResult {
    Cipher slotCipherFirst;
    Cipher slotCipherSecond;
};

/**
 * Canonical bootstrap CoeffToSlot для полного CKKS packing.
 *
 * План гомоморфно вычисляет две половины обратного canonical embedding через
 * диагональные матрицы, Galois automorphisms и conjugated sums. Входом может
 * быть только RaisedCipher.
 */
class CoeffToSlot {
public:
    explicit CoeffToSlot(std::size_t polyModulusDegree);
    ~CoeffToSlot();
    CoeffToSlot(CoeffToSlot&&) noexcept;
    CoeffToSlot& operator=(CoeffToSlot&&) noexcept;
    CoeffToSlot(const CoeffToSlot&) = delete;
    CoeffToSlot& operator=(const CoeffToSlot&) = delete;

    /// Требуется один rescale-уровень, все ненулевые rotations и conjugation key.
    CoeffToSlotPlanRequirements requirements() const;

    /// Потребляет RaisedCipher и возвращает две coefficient-side slot упаковки.
    CoeffToSlotResult apply(SealAdapter& adapter, RaisedCipher&& input) const;

private:
    struct Impl;
    std::unique_ptr<Impl> pimpl_;
};

} // namespace m2424
