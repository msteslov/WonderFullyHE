#pragma once

#include <complex>
#include <cstddef>
#include <vector>

namespace m2424 {

using Complex = std::complex<double>;
using ComplexVector = std::vector<Complex>;

/**
 * Возвращает порядок нечётных степеней корня, соответствующий порядку CKKS-слотов Microsoft SEAL.
 *
 * Для кольца X^N + 1 SEAL последовательно использует степени p[i], где p[0] = 1 и
 * p[i + 1] = 3 * p[i] mod 2N. Bit-reversal в исходном коде SEAL относится только
 * к внутреннему расположению данных FFT и не меняет математический порядок результата.
 */
std::vector<std::size_t> canonicalEmbeddingRootExponents(std::size_t polyModulusDegree);

/**
 * Вычисляет canonical embedding нормализованного коэффициентного вектора в CKKS-слоты.
 *
 * Вход содержит N вещественных коэффициентов полинома из R[X]/(X^N + 1), уже делённых
 * на входной CKKS scale. Функция является O(N^2) plaintext-эталоном и не предназначена
 * для гомоморфного или production-выполнения.
 */
ComplexVector coeffToSlotReference(const std::vector<double>& normalizedCoefficients);

/**
 * Выполняет обратный canonical embedding и возвращает нормализованные вещественные коэффициенты.
 *
 * Вход должен содержать ровно N/2 комплексных слотов в порядке Microsoft SEAL.
 */
std::vector<double> slotToCoeffReference(const ComplexVector& slots);

} // namespace m2424
