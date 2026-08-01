#pragma once

#include "m2424/evalmod_domain_analysis.hpp"

#include <vector>

namespace m2424 {

struct EvalModErrorCertificate {
    double approximationMaxError{};
    double derivativeMax{};
    double complexNeighborhoodError{};
    double arithmeticErrorEstimate{};
    double totalPredictedError{};
};

/// Сертифицирует monomial polynomial на объединении [k-rho,k+rho], |k|<=K.
/// Значения между интервалами намеренно не входят в real-domain maximum.
EvalModErrorCertificate certifyEvalModPolynomial(const std::vector<double>& coefficients,
                                                 const EvalModDomain& domain,
                                                 std::size_t samplesPerInterval,
                                                 double complexNeighborhoodRadius,
                                                 double arithmeticErrorEstimate);

} // namespace m2424
