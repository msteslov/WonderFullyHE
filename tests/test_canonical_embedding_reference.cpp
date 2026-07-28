#include "m2424/canonical_embedding_reference.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <vector>

namespace {

constexpr double kReferenceTolerance = 1e-12;
constexpr double kPi = 3.141592653589793238462643383279502884;

bool close(double lhs, double rhs) {
    return std::abs(lhs - rhs) <= kReferenceTolerance;
}

bool close(m2424::Complex lhs, m2424::Complex rhs) {
    return std::abs(lhs - rhs) <= kReferenceTolerance;
}

bool rejectsInvalidDegree() {
    try {
        (void)m2424::canonicalEmbeddingRootExponents(12);
    } catch (const std::invalid_argument&) {
        return true;
    }
    return false;
}

} // namespace

int main() {
    bool ok = true;

    const auto exponents = m2424::canonicalEmbeddingRootExponents(8);
    ok = ok && exponents == std::vector<std::size_t>({1, 3, 9, 11});

    std::vector<double> constant(8, 0.0);
    constant[0] = 1.0;
    const auto constantSlots = m2424::coeffToSlotReference(constant);
    for (const auto value : constantSlots) {
        ok = ok && close(value, {1.0, 0.0});
    }

    std::vector<double> linear(8, 0.0);
    linear[1] = 1.0;
    const auto linearSlots = m2424::coeffToSlotReference(linear);
    for (std::size_t index = 0; index < linearSlots.size(); ++index) {
        const double angle = (2.0 * kPi * static_cast<double>(exponents[index])) / 16.0;
        ok = ok && close(linearSlots[index], {std::cos(angle), std::sin(angle)});
    }

    const std::vector<double> coefficients{
        -1.0, 0.75, -0.5, 0.25, 0.0, -0.125, 0.625, -0.875,
        0.375, -0.25, 0.5, -0.75, 1.0, -0.625, 0.125, 0.0
    };
    const auto roundtrip = m2424::slotToCoeffReference(m2424::coeffToSlotReference(coefficients));
    double maxError = 0.0;
    for (std::size_t index = 0; index < coefficients.size(); ++index) {
        maxError = std::max(maxError, std::abs(coefficients[index] - roundtrip[index]));
    }
    ok = ok && maxError <= kReferenceTolerance && rejectsInvalidDegree();

    std::printf("[test_canonical_embedding_reference] max_abs_error=%.3e %s\n",
                maxError, ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
