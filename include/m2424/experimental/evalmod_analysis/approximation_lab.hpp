#pragma once

#include "m2424/experimental/evalmod_analysis/domain_analysis.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace m2424::experimental {

enum class PolynomialBasis { Monomial, Chebyshev, Composite };

struct EvalModPolynomial {
    PolynomialBasis basis{PolynomialBasis::Monomial};
    std::vector<std::string> decimalCoefficients;
};

struct EvalModGridDiagnostic {
    double approximationMaxError{};
    double realDerivativeMax{};
    double complexDerivativeMax{};
    double complexBoundaryErrorMax{};
    double arithmeticErrorEstimate{};
    std::size_t evaluations{};
    std::size_t precisionBits{};
};

/// MPFR grid diagnostic; это не interval certificate.
EvalModGridDiagnostic diagnoseEvalModPolynomialOnGrid(const EvalModPolynomial& polynomial,
                                                      const EvalModDomain& domain,
                                                      std::size_t samplesPerEdge,
                                                      const std::string& complexRadiusDecimal,
                                                      double arithmeticErrorEstimate,
                                                      std::size_t precisionBits = 256,
                                                      std::size_t maxEvaluations = 1'000'000);

} // namespace m2424::experimental
