#include "m2424/canonical_embedding_reference.hpp"
#include "m2424/coeff_to_slot.hpp"

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

double factorizedError(const std::vector<double>& coefficients, std::size_t depth) {
    m2424::CoeffToSlotPlan plan(coefficients.size(), depth);
    const auto transformed = plan.applyPlain(m2424::coeffToSlotReference(coefficients));
    const std::size_t slots = coefficients.size() / 2;
    double error = 0.0;
    for (std::size_t index = 0; index < slots; ++index) {
        error = std::max(error, std::abs(transformed.first[index] - coefficients[index]));
        error = std::max(error, std::abs(transformed.second[index] - coefficients[index + slots]));
    }
    return error;
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

    std::vector<double> basis(32);
    basis[0] = 1.0;
    const double basisError = factorizedError(basis, 3);
    basis.assign(32, 0.0);
    basis[31] = -1.0;
    const double boundaryError = factorizedError(basis, 3);
    std::vector<double> deterministic(32);
    for (std::size_t index = 0; index < deterministic.size(); ++index) {
        deterministic[index] =
            static_cast<double>(static_cast<int>((index * 29) % 41) - 20) / 32.0;
    }
    double factorizedRandomError = 0.0;
    for (std::size_t depth = 1; depth <= 5; ++depth) {
        factorizedRandomError =
            std::max(factorizedRandomError, factorizedError(deterministic, depth));
    }
    ok = ok && basisError <= kReferenceTolerance
        && boundaryError <= kReferenceTolerance
        && factorizedRandomError <= kReferenceTolerance;

    std::printf("[test_canonical_embedding_reference] roundtrip=%.3e factorized=%.3e %s\n",
                maxError, std::max({basisError, boundaryError, factorizedRandomError}),
                ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
