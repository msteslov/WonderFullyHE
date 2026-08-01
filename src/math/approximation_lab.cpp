#include "m2424/experimental/evalmod_analysis/approximation_lab.hpp"

#include <mpfr.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>

namespace m2424::experimental {
namespace {

class Real {
public:
    explicit Real(mpfr_prec_t precision) { mpfr_init2(value_, precision); }
    Real(const Real& other) { mpfr_init2(value_, mpfr_get_prec(other.value_)); mpfr_set(value_, other.value_, MPFR_RNDN); }
    ~Real() { mpfr_clear(value_); }
    mpfr_ptr get() { return value_; }
    mpfr_srcptr get() const { return value_; }
private:
    mpfr_t value_;
};

struct ComplexValue { Real real; Real imag; explicit ComplexValue(mpfr_prec_t p) : real(p), imag(p) {} };

void evaluate(const std::vector<Real>& coefficients, const Real& xr, const Real& xi,
              ComplexValue& value, ComplexValue& derivative) {
    const mpfr_prec_t precision = mpfr_get_prec(xr.get());
    Real nextReal(precision), nextImag(precision), a(precision), b(precision);
    mpfr_set_zero(value.real.get(), 0); mpfr_set_zero(value.imag.get(), 0);
    mpfr_set_zero(derivative.real.get(), 0); mpfr_set_zero(derivative.imag.get(), 0);
    for (std::size_t index = coefficients.size(); index-- > 0;) {
        mpfr_mul(a.get(), derivative.real.get(), xr.get(), MPFR_RNDN);
        mpfr_mul(b.get(), derivative.imag.get(), xi.get(), MPFR_RNDN);
        mpfr_sub(nextReal.get(), a.get(), b.get(), MPFR_RNDN);
        mpfr_mul(a.get(), derivative.real.get(), xi.get(), MPFR_RNDN);
        mpfr_mul(b.get(), derivative.imag.get(), xr.get(), MPFR_RNDN);
        mpfr_add(nextImag.get(), a.get(), b.get(), MPFR_RNDN);
        mpfr_add(nextReal.get(), nextReal.get(), value.real.get(), MPFR_RNDN);
        mpfr_add(nextImag.get(), nextImag.get(), value.imag.get(), MPFR_RNDN);
        mpfr_set(derivative.real.get(), nextReal.get(), MPFR_RNDN);
        mpfr_set(derivative.imag.get(), nextImag.get(), MPFR_RNDN);

        mpfr_mul(a.get(), value.real.get(), xr.get(), MPFR_RNDN);
        mpfr_mul(b.get(), value.imag.get(), xi.get(), MPFR_RNDN);
        mpfr_sub(nextReal.get(), a.get(), b.get(), MPFR_RNDN);
        mpfr_mul(a.get(), value.real.get(), xi.get(), MPFR_RNDN);
        mpfr_mul(b.get(), value.imag.get(), xr.get(), MPFR_RNDN);
        mpfr_add(nextImag.get(), a.get(), b.get(), MPFR_RNDN);
        mpfr_add(nextReal.get(), nextReal.get(), coefficients[index].get(), MPFR_RNDN);
        mpfr_set(value.real.get(), nextReal.get(), MPFR_RNDN);
        mpfr_set(value.imag.get(), nextImag.get(), MPFR_RNDN);
    }
}

double magnitudeUp(const Real& real, const Real& imag) {
    Real magnitude(mpfr_get_prec(real.get()));
    mpfr_hypot(magnitude.get(), real.get(), imag.get(), MPFR_RNDU);
    return mpfr_get_d(magnitude.get(), MPFR_RNDU);
}

} // namespace

EvalModGridDiagnostic diagnoseEvalModPolynomialOnGrid(const EvalModPolynomial& polynomial,
                                                      const EvalModDomain& domain,
                                                      std::size_t samplesPerEdge,
                                                      const std::string& radiusDecimal,
                                                      double arithmeticErrorEstimate,
                                                      std::size_t precisionBits,
                                                      std::size_t maxEvaluations) {
    if (polynomial.basis != PolynomialBasis::Monomial || polynomial.decimalCoefficients.empty()
        || !isEvalModDomainValid(domain) || samplesPerEdge < 2 || precisionBits < 64
        || maxEvaluations == 0 || !std::isfinite(arithmeticErrorEstimate)
        || arithmeticErrorEstimate < 0.0
        || domain.integerBound > static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max() - 1)) {
        throw std::invalid_argument("invalid EvalMod grid diagnostic input");
    }
    const std::size_t intervals = domain.integerBound * 2 + 1;
    if (samplesPerEdge > std::numeric_limits<std::size_t>::max() / 5
        || intervals > maxEvaluations / (samplesPerEdge * 5)) {
        throw std::length_error("EvalMod grid diagnostic evaluation limit exceeded");
    }
    const auto precision = static_cast<mpfr_prec_t>(precisionBits);
    std::vector<Real> coefficients;
    coefficients.reserve(polynomial.decimalCoefficients.size());
    for (const auto& decimal : polynomial.decimalCoefficients) {
        coefficients.emplace_back(precision);
        if (mpfr_set_str(coefficients.back().get(), decimal.c_str(), 10, MPFR_RNDN) != 0) {
            throw std::invalid_argument("invalid decimal polynomial coefficient");
        }
    }
    Real radius(precision);
    if (mpfr_set_str(radius.get(), radiusDecimal.c_str(), 10, MPFR_RNDU) != 0 || mpfr_sgn(radius.get()) < 0) {
        throw std::invalid_argument("invalid complex radius");
    }

    EvalModGridDiagnostic result;
    result.arithmeticErrorEstimate = arithmeticErrorEstimate;
    result.precisionBits = precisionBits;
    Real xr(precision), xi(precision), expectedReal(precision), errorReal(precision), errorImag(precision);
    ComplexValue value(precision), derivative(precision);
    const auto observe = [&](double real, double imag, std::int64_t integer, bool realAxis) {
        mpfr_set_d(xr.get(), real, MPFR_RNDN); mpfr_set_d(xi.get(), imag, MPFR_RNDN);
        evaluate(coefficients, xr, xi, value, derivative);
        mpfr_set_si(expectedReal.get(), integer, MPFR_RNDN);
        mpfr_sub(errorReal.get(), value.real.get(), xr.get(), MPFR_RNDN);
        mpfr_add(errorReal.get(), errorReal.get(), expectedReal.get(), MPFR_RNDN);
        mpfr_sub(errorImag.get(), value.imag.get(), xi.get(), MPFR_RNDN);
        const double error = magnitudeUp(errorReal, errorImag);
        const double derivativeMagnitude = magnitudeUp(derivative.real, derivative.imag);
        if (realAxis) {
            result.approximationMaxError = std::max(result.approximationMaxError, error);
            result.realDerivativeMax = std::max(result.realDerivativeMax, derivativeMagnitude);
        } else {
            result.complexBoundaryErrorMax = std::max(result.complexBoundaryErrorMax, error);
            result.complexDerivativeMax = std::max(result.complexDerivativeMax, derivativeMagnitude);
        }
        ++result.evaluations;
    };

    const double rho = domain.normalizedResidualBound;
    const double radiusValue = mpfr_get_d(radius.get(), MPFR_RNDU);
    for (std::int64_t integer = -static_cast<std::int64_t>(domain.integerBound);;
         ++integer) {
        for (std::size_t sample = 0; sample < samplesPerEdge; ++sample) {
            const double t = static_cast<double>(sample) / static_cast<double>(samplesPerEdge - 1);
            const double real = static_cast<double>(integer) - rho + 2.0 * rho * t;
            observe(real, 0.0, integer, true);
            observe(real, -radiusValue, integer, false);
            observe(real, radiusValue, integer, false);
            const double imag = -radiusValue + 2.0 * radiusValue * t;
            observe(static_cast<double>(integer) - rho, imag, integer, false);
            observe(static_cast<double>(integer) + rho, imag, integer, false);
        }
        if (integer == static_cast<std::int64_t>(domain.integerBound)) break;
    }
    return result;
}

} // namespace m2424::experimental
