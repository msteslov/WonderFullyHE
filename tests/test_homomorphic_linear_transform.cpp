#include "m2424/accuracy.hpp"
#include "m2424/homomorphic_linear_transform.hpp"
#include "m2424/profiles.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <vector>

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

std::vector<m2424::HomomorphicDiagonalTerm> fourPointDftTerms(std::size_t physicalSlots) {
    constexpr std::size_t kSize = 4;
    std::map<int, m2424::ComplexVector> diagonals;
    for (std::size_t row = 0; row < kSize; ++row) {
        for (std::size_t column = 0; column < kSize; ++column) {
            const int rotation = static_cast<int>(column) - static_cast<int>(row);
            auto [it, inserted] = diagonals.emplace(rotation, m2424::ComplexVector(physicalSlots));
            const double angle = -2.0 * kPi * static_cast<double>(row * column) / static_cast<double>(kSize);
            it->second[row] += m2424::Complex{std::cos(angle), std::sin(angle)};
        }
    }

    std::vector<m2424::HomomorphicDiagonalTerm> terms;
    terms.reserve(diagonals.size());
    for (auto& [rotation, diagonal] : diagonals) {
        terms.push_back({rotation, std::move(diagonal)});
    }
    return terms;
}

m2424::ComplexVector fourPointDftReference(const m2424::ComplexVector& input) {
    m2424::ComplexVector result;
    result.reserve(input.size());
    for (std::size_t row = 0; row < input.size(); ++row) {
        m2424::Complex value{};
        for (std::size_t column = 0; column < input.size(); ++column) {
            const double angle = -2.0 * kPi * static_cast<double>(row * column) / static_cast<double>(input.size());
            value += input[column] * m2424::Complex{std::cos(angle), std::sin(angle)};
        }
        result.push_back(value);
    }
    return result;
}

} // namespace

int main() {
    auto adapter = m2424::SealAdapter::create(m2424::profiles::high_precision_ckks());
    const auto terms = fourPointDftTerms(adapter.slotCount());
    const m2424::HomomorphicLinearTransform transform(adapter.slotCount(), terms);
    adapter.generateKeys(transform.rotationSteps(), false);

    const m2424::ComplexVector input{{0.25, -0.5}, {-0.75, 0.125}, {0.5, 0.25}, {-0.125, -0.375}};
    const auto encrypted = adapter.encrypt(adapter.encodeComplex(input));
    const auto result = transform.apply(adapter, encrypted);
    const auto inputInfo = adapter.info(encrypted);
    const auto resultInfo = adapter.info(result);
    const auto decoded = adapter.decodeComplex(adapter.decrypt(result));
    const auto expected = fourPointDftReference(input);

    double maxError = 0.0;
    for (std::size_t index = 0; index < input.size(); ++index) {
        maxError = std::max(maxError, std::abs(expected[index] - decoded[index]));
    }
    const bool ok = maxError <= m2424::kTargetAbsoluteError
        && resultInfo.chainIndex + 1 == inputInfo.chainIndex
        && std::abs(std::log2(resultInfo.scale) - std::log2(inputInfo.scale)) <= 0.25
        && transform.rotationSteps() == std::vector<int>({-3, -2, -1, 1, 2, 3});
    std::printf("[test_homomorphic_linear_transform] max_abs_error=%.3e %s\n",
                maxError, ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
