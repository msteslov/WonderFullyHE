#include "m2424/canonical_embedding_reference.hpp"

#include <cmath>
#include <stdexcept>

namespace m2424 {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

bool isPowerOfTwo(std::size_t value) {
    return value >= 2 && (value & (value - 1)) == 0;
}

void validatePolyModulusDegree(std::size_t polyModulusDegree) {
    if (!isPowerOfTwo(polyModulusDegree)) {
        throw std::invalid_argument("polyModulusDegree must be a power of two greater than one");
    }
}

Complex rootOfUnity(std::size_t exponent, std::size_t order) {
    const double angle = (2.0 * kPi * static_cast<double>(exponent)) / static_cast<double>(order);
    return {std::cos(angle), std::sin(angle)};
}

} // namespace

std::vector<std::size_t> canonicalEmbeddingRootExponents(std::size_t polyModulusDegree) {
    validatePolyModulusDegree(polyModulusDegree);

    const std::size_t slotCount = polyModulusDegree / 2;
    const std::size_t rootOrder = 2 * polyModulusDegree;
    std::vector<std::size_t> exponents;
    exponents.reserve(slotCount);

    std::size_t position = 1;
    for (std::size_t index = 0; index < slotCount; ++index) {
        exponents.push_back(position);
        position = (position * 3) % rootOrder;
    }
    return exponents;
}

ComplexVector coeffToSlotReference(const std::vector<double>& normalizedCoefficients) {
    const std::size_t polyModulusDegree = normalizedCoefficients.size();
    validatePolyModulusDegree(polyModulusDegree);

    const auto exponents = canonicalEmbeddingRootExponents(polyModulusDegree);
    const std::size_t rootOrder = 2 * polyModulusDegree;
    ComplexVector slots;
    slots.reserve(exponents.size());

    for (const std::size_t exponent : exponents) {
        Complex value{};
        for (std::size_t coefficient = 0; coefficient < polyModulusDegree; ++coefficient) {
            value += normalizedCoefficients[coefficient] * rootOfUnity(exponent * coefficient, rootOrder);
        }
        slots.push_back(value);
    }
    return slots;
}

std::vector<double> slotToCoeffReference(const ComplexVector& slots) {
    const std::size_t slotCount = slots.size();
    const std::size_t polyModulusDegree = 2 * slotCount;
    validatePolyModulusDegree(polyModulusDegree);

    const auto exponents = canonicalEmbeddingRootExponents(polyModulusDegree);
    const std::size_t rootOrder = 2 * polyModulusDegree;
    std::vector<double> coefficients(polyModulusDegree, 0.0);

    for (std::size_t coefficient = 0; coefficient < polyModulusDegree; ++coefficient) {
        Complex sum{};
        for (std::size_t index = 0; index < slotCount; ++index) {
            const std::size_t exponent = exponents[index];
            sum += slots[index] * rootOfUnity(rootOrder - ((exponent * coefficient) % rootOrder), rootOrder);
            sum += std::conj(slots[index]) * rootOfUnity(exponent * coefficient, rootOrder);
        }
        coefficients[coefficient] = (sum / static_cast<double>(polyModulusDegree)).real();
    }
    return coefficients;
}

} // namespace m2424
