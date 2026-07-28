#pragma once

#include "m2424/canonical_embedding_reference.hpp"
#include "m2424/seal_adapter.hpp"

#include <cstddef>
#include <vector>

namespace m2424 {

/// Одна диагональ линейного преобразования: diagonal[i] * rotate(input, rotation)[i].
struct HomomorphicDiagonalTerm {
    int rotation{};
    ComplexVector diagonal;
};

/**
 * Исполняет заданное линейное преобразование CKKS-слотов через ротации и plaintext-диагонали.
 *
 * Класс не выбирает стратегию CoeffToSlot и не меняет scale вручную. FFT и BSGS должны
 * строить свои планы из этих диагональных terms и проходить отдельную проверку точности.
 */
class HomomorphicLinearTransform {
public:
    /// Создаёт преобразование для физической ёмкости CKKS-слотов.
    HomomorphicLinearTransform(std::size_t physicalSlotCount, std::vector<HomomorphicDiagonalTerm> terms);

    /// Возвращает шаги ротации, для которых требуются Galois keys.
    std::vector<int> rotationSteps() const;

    /// Применяет линейное преобразование с одним общим уровнем: каждая диагональ умножается и rescale-ится.
    Cipher apply(SealAdapter& adapter, const Cipher& input) const;

private:
    std::size_t physicalSlotCount_{};
    std::vector<HomomorphicDiagonalTerm> terms_;
};

} // namespace m2424
