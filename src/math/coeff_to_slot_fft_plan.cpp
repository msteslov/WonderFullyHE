#include "m2424/coeff_to_slot_fft_plan.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <utility>

namespace m2424 {
namespace {

using ComplexMatrix = std::vector<ComplexVector>;
constexpr double kPi = 3.141592653589793238462643383279502884;

bool isPowerOfTwo(std::size_t value) {
    return value >= 2 && (value & (value - 1)) == 0;
}

std::size_t reverseBits(std::size_t value, std::size_t bitCount) {
    std::size_t result = 0;
    for (std::size_t bit = 0; bit < bitCount; ++bit) {
        result = (result << 1) | (value & 1U);
        value >>= 1U;
    }
    return result;
}

std::size_t log2Exact(std::size_t value) {
    std::size_t result = 0;
    while (value > 1) {
        value >>= 1U;
        ++result;
    }
    return result;
}

ComplexMatrix emptyMatrix(std::size_t size) {
    return ComplexMatrix(size, ComplexVector(size));
}

ComplexMatrix bitReversalMatrix(std::size_t size) {
    auto matrix = emptyMatrix(size);
    const auto bits = log2Exact(size);
    for (std::size_t row = 0; row < size; ++row) {
        matrix[row][reverseBits(row, bits)] = {1.0, 0.0};
    }
    return matrix;
}

ComplexMatrix butterflyMatrix(std::size_t size, std::size_t width) {
    auto matrix = emptyMatrix(size);
    const std::size_t half = width / 2;
    for (std::size_t block = 0; block < size; block += width) {
        for (std::size_t offset = 0; offset < half; ++offset) {
            const std::size_t low = block + offset;
            const std::size_t high = low + half;
            const double angle = -2.0 * kPi * static_cast<double>(offset) / static_cast<double>(width);
            const Complex twiddle{std::cos(angle), std::sin(angle)};
            matrix[low][low] = {1.0, 0.0};
            matrix[low][high] = twiddle;
            matrix[high][low] = {1.0, 0.0};
            matrix[high][high] = -twiddle;
        }
    }
    return matrix;
}

HomomorphicLinearTransform matrixToTransform(const ComplexMatrix& matrix, std::size_t physicalSlotCount) {
    const std::size_t size = matrix.size();
    std::map<int, ComplexVector> diagonals;
    for (std::size_t row = 0; row < size; ++row) {
        if (matrix[row].size() != size) {
            throw std::invalid_argument("FFT layer matrix must be square");
        }
        for (std::size_t column = 0; column < size; ++column) {
            if (std::abs(matrix[row][column]) == 0.0) {
                continue;
            }
            const int rotation = static_cast<int>(column) - static_cast<int>(row);
            auto [it, inserted] = diagonals.emplace(rotation, ComplexVector(physicalSlotCount));
            it->second[row] += matrix[row][column];
        }
    }

    std::vector<HomomorphicDiagonalTerm> terms;
    terms.reserve(diagonals.size());
    for (auto& [rotation, diagonal] : diagonals) {
        terms.push_back({rotation, std::move(diagonal)});
    }
    return {physicalSlotCount, std::move(terms)};
}

ComplexVector applyButterflyPlain(ComplexVector values, std::size_t width) {
    const std::size_t half = width / 2;
    for (std::size_t block = 0; block < values.size(); block += width) {
        for (std::size_t offset = 0; offset < half; ++offset) {
            const std::size_t low = block + offset;
            const std::size_t high = low + half;
            const double angle = -2.0 * kPi * static_cast<double>(offset) / static_cast<double>(width);
            const Complex twiddle{std::cos(angle), std::sin(angle)};
            const Complex lhs = values[low];
            const Complex rhs = twiddle * values[high];
            values[low] = lhs + rhs;
            values[high] = lhs - rhs;
        }
    }
    return values;
}

} // namespace

CoeffToSlotFftPlan::CoeffToSlotFftPlan(std::size_t transformSlots, std::size_t physicalSlotCount)
    : transformSlots_(transformSlots) {
    if (!isPowerOfTwo(transformSlots_) || transformSlots_ > physicalSlotCount) {
        throw std::invalid_argument("FFT transform slots must be a power of two within physical slot count");
    }
    layers_.push_back(matrixToTransform(bitReversalMatrix(transformSlots_), physicalSlotCount));
    for (std::size_t width = 2; width <= transformSlots_; width *= 2) {
        layers_.push_back(matrixToTransform(butterflyMatrix(transformSlots_, width), physicalSlotCount));
    }
}

std::size_t CoeffToSlotFftPlan::requiredLevels() const noexcept {
    return layers_.size();
}

std::vector<int> CoeffToSlotFftPlan::rotationSteps() const {
    std::vector<int> steps;
    for (const auto& layer : layers_) {
        const auto layerSteps = layer.rotationSteps();
        steps.insert(steps.end(), layerSteps.begin(), layerSteps.end());
    }
    std::sort(steps.begin(), steps.end());
    steps.erase(std::unique(steps.begin(), steps.end()), steps.end());
    return steps;
}

ComplexVector CoeffToSlotFftPlan::applyPlain(const ComplexVector& input) const {
    if (input.size() != transformSlots_) {
        throw std::invalid_argument("FFT plaintext input size must match transform slots");
    }
    ComplexVector result(transformSlots_);
    const auto bits = log2Exact(transformSlots_);
    for (std::size_t index = 0; index < transformSlots_; ++index) {
        result[index] = input[reverseBits(index, bits)];
    }
    for (std::size_t width = 2; width <= transformSlots_; width *= 2) {
        result = applyButterflyPlain(std::move(result), width);
    }
    return result;
}

Cipher CoeffToSlotFftPlan::apply(SealAdapter& adapter, const Cipher& input) const {
    Cipher result = input;
    for (const auto& layer : layers_) {
        result = layer.apply(adapter, result);
    }
    return result;
}

CoeffToSlotPlanRequirements CoeffToSlotFftPlan::requirements() const {
    return {requiredLevels(), rotationSteps()};
}

} // namespace m2424
