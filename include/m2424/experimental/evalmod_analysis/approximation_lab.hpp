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

struct EvalModIntervalCertificate {
    std::string approximationErrorUpperBound;
    std::string derivativeUpperBound;
    std::string complexErrorUpperBound;
    double approximationErrorUpperBoundDouble{};
    double derivativeUpperBoundDouble{};
    double complexErrorUpperBoundDouble{};
    bool proved{};
};

/// Converts standard Chebyshev coefficients (sum c_k T_k(x)) to an exact
/// executable monomial polynomial. Decimal inputs are transformed as rationals.
EvalModPolynomial convertChebyshevToMonomial(const EvalModPolynomial& polynomial);

/// MPFR grid diagnostic; это не interval certificate.
EvalModGridDiagnostic diagnoseEvalModPolynomialOnGrid(const EvalModPolynomial& polynomial,
                                                      const EvalModDomain& domain,
                                                      std::size_t samplesPerEdge,
                                                      const std::string& complexRadiusDecimal,
                                                      double arithmeticErrorEstimate,
                                                      std::size_t precisionBits = 256,
                                                      std::size_t maxEvaluations = 1'000'000);

EvalModIntervalCertificate certifyEvalModPolynomialIntervals(
    const EvalModPolynomial& polynomial, const EvalModDomain& domain,
    const std::string& complexRadiusDecimal, std::size_t subdivisions = 4096,
    std::size_t precisionBits = 384);

} // namespace m2424::experimental
