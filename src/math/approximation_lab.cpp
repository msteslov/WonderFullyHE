#include "m2424/approximation_lab.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <stdexcept>

namespace m2424 {
namespace {

template <class T>
T evaluate(const std::vector<double>& coefficients, T x) {
    T result{};
    for (auto it = coefficients.rbegin(); it != coefficients.rend(); ++it) result = result * x + *it;
    return result;
}

double evaluateDerivative(const std::vector<double>& coefficients, double x) {
    double result = 0.0;
    for (std::size_t i = coefficients.size(); i-- > 1;) result = result * x + i * coefficients[i];
    return result;
}

} // namespace

EvalModErrorCertificate certifyEvalModPolynomial(const std::vector<double>& coefficients,
                                                 const EvalModDomain& domain,
                                                 std::size_t samplesPerInterval,
                                                 double complexNeighborhoodRadius,
                                                 double arithmeticErrorEstimate) {
    if (coefficients.empty() || !isEvalModDomainValid(domain) || samplesPerInterval < 2
        || !std::isfinite(complexNeighborhoodRadius) || complexNeighborhoodRadius < 0.0
        || !std::isfinite(arithmeticErrorEstimate) || arithmeticErrorEstimate < 0.0) {
        throw std::invalid_argument("invalid approximation certificate input");
    }
    for (double coefficient : coefficients) {
        if (!std::isfinite(coefficient)) throw std::invalid_argument("non-finite polynomial coefficient");
    }

    EvalModErrorCertificate result;
    for (std::int64_t integer = -static_cast<std::int64_t>(domain.integerBound);
         integer <= static_cast<std::int64_t>(domain.integerBound); ++integer) {
        for (std::size_t sample = 0; sample < samplesPerInterval; ++sample) {
            const double residual = -domain.normalizedResidualBound
                + 2.0 * domain.normalizedResidualBound * static_cast<double>(sample)
                    / static_cast<double>(samplesPerInterval - 1);
            const double x = static_cast<double>(integer) + residual;
            result.approximationMaxError = std::max(result.approximationMaxError,
                                                     std::abs(evaluate(coefficients, x) - residual));
            result.derivativeMax = std::max(result.derivativeMax,
                                             std::abs(evaluateDerivative(coefficients, x)));
            for (double sign : {-1.0, 1.0}) {
                const std::complex<double> z{x, sign * complexNeighborhoodRadius};
                const std::complex<double> expected{residual, sign * complexNeighborhoodRadius};
                result.complexNeighborhoodError = std::max(result.complexNeighborhoodError,
                                                             std::abs(evaluate(coefficients, z) - expected));
            }
        }
    }
    result.arithmeticErrorEstimate = arithmeticErrorEstimate;
    result.totalPredictedError = result.approximationMaxError + arithmeticErrorEstimate;
    return result;
}

} // namespace m2424
