#include "m2424/evalmod_domain_analysis.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace m2424 {

EvalModDomain analyzeEvalModDomain(const EvalModCiphertextModel& model) {
    if (model.coefficientCount == 0 || !std::isfinite(model.deterministicIntegerOffset)
        || model.deterministicIntegerOffset < 0.0 || !std::isfinite(model.integerNoiseStddev)
        || model.integerNoiseStddev < 0.0 || !std::isfinite(model.normalizedMessageAbsBound)
        || model.normalizedMessageAbsBound < 0.0
        || !std::isfinite(model.normalizedEncodingErrorAbsBound)
        || model.normalizedEncodingErrorAbsBound < 0.0
        || !std::isfinite(model.failureProbabilityLog2)
        || model.failureProbabilityLog2 >= 0.0) {
        throw std::invalid_argument("invalid EvalMod ciphertext model");
    }

    // P(max_i |X_i| > t) <= 2*n*exp(-t^2/(2*sigma^2)).
    const double logFailure = model.failureProbabilityLog2 * std::log(2.0);
    const double logUnionFactor = std::log(2.0 * static_cast<double>(model.coefficientCount));
    const double tail = model.integerNoiseStddev
        * std::sqrt(2.0 * (logUnionFactor - logFailure));
    const double integerBound = std::ceil(model.deterministicIntegerOffset + tail);
    if (!std::isfinite(integerBound)
        || integerBound > static_cast<double>(std::numeric_limits<std::size_t>::max())) {
        throw std::overflow_error("EvalMod integer bound does not fit size_t");
    }

    const double rho = model.normalizedMessageAbsBound + model.normalizedEncodingErrorAbsBound;
    EvalModDomain result{static_cast<std::size_t>(integerBound), rho, 0.5 - rho,
                         model.failureProbabilityLog2};
    if (!isEvalModDomainValid(result)) {
        throw std::domain_error("EvalMod residual reaches a rounding discontinuity");
    }
    return result;
}

bool isEvalModDomainValid(const EvalModDomain& domain) {
    return std::isfinite(domain.normalizedResidualBound)
        && domain.normalizedResidualBound >= 0.0 && domain.normalizedResidualBound < 0.5
        && std::isfinite(domain.discontinuityMargin) && domain.discontinuityMargin > 0.0
        && std::abs(domain.discontinuityMargin - (0.5 - domain.normalizedResidualBound)) < 1e-12
        && std::isfinite(domain.failureProbabilityLog2) && domain.failureProbabilityLog2 < 0.0;
}

} // namespace m2424
