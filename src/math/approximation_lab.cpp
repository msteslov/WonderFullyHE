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
    if (!mpfr_number_p(real.get()) || !mpfr_number_p(imag.get())) {
        throw std::overflow_error("non-finite MPFR diagnostic result");
    }
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
        if (mpfr_set_str(coefficients.back().get(), decimal.c_str(), 10, MPFR_RNDN) != 0
            || !mpfr_number_p(coefficients.back().get())) {
            throw std::invalid_argument("invalid decimal polynomial coefficient");
        }
    }
    Real radius(precision);
    if (mpfr_set_str(radius.get(), radiusDecimal.c_str(), 10, MPFR_RNDU) != 0
        || !mpfr_number_p(radius.get()) || mpfr_sgn(radius.get()) < 0) {
        throw std::invalid_argument("invalid complex radius");
    }
    Real rho(precision);
    if (mpfr_set_str(rho.get(), domain.normalizedResidualBoundDecimal.c_str(), 10, MPFR_RNDU) != 0
        || !mpfr_number_p(rho.get()) || mpfr_sgn(rho.get()) < 0) {
        throw std::invalid_argument("invalid high-precision residual bound");
    }

    EvalModGridDiagnostic result;
    result.arithmeticErrorEstimate = arithmeticErrorEstimate;
    result.precisionBits = precisionBits;
    Real xr(precision), xi(precision), expectedReal(precision), errorReal(precision), errorImag(precision);
    ComplexValue value(precision), derivative(precision);
    const auto observe = [&](const Real& real, const Real& imag, std::int64_t integer, bool realAxis) {
        mpfr_set(xr.get(), real.get(), MPFR_RNDN); mpfr_set(xi.get(), imag.get(), MPFR_RNDN);
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

    Real zero(precision), t(precision), real(precision), imag(precision), twiceRho(precision), integerMpfr(precision);
    mpfr_set_zero(zero.get(), 0);
    mpfr_mul_ui(twiceRho.get(), rho.get(), 2, MPFR_RNDU);
    for (std::int64_t integer = -static_cast<std::int64_t>(domain.integerBound);;
         ++integer) {
        mpfr_set_si(integerMpfr.get(), integer, MPFR_RNDN);
        for (std::size_t sample = 0; sample < samplesPerEdge; ++sample) {
            mpfr_set_ui(t.get(), sample, MPFR_RNDN);
            mpfr_div_ui(t.get(), t.get(), samplesPerEdge - 1, MPFR_RNDN);
            mpfr_mul(real.get(), twiceRho.get(), t.get(), MPFR_RNDN);
            mpfr_sub(real.get(), real.get(), rho.get(), MPFR_RNDN);
            mpfr_add(real.get(), real.get(), integerMpfr.get(), MPFR_RNDN);
            observe(real, zero, integer, true);
            mpfr_neg(imag.get(), radius.get(), MPFR_RNDN);
            observe(real, imag, integer, false);
            observe(real, radius, integer, false);
            mpfr_mul_ui(imag.get(), radius.get(), 2, MPFR_RNDN);
            mpfr_mul(imag.get(), imag.get(), t.get(), MPFR_RNDN);
            mpfr_sub(imag.get(), imag.get(), radius.get(), MPFR_RNDN);
            mpfr_sub(real.get(), integerMpfr.get(), rho.get(), MPFR_RNDN);
            observe(real, imag, integer, false);
            mpfr_add(real.get(), integerMpfr.get(), rho.get(), MPFR_RNDN);
            observe(real, imag, integer, false);
        }
        if (integer == static_cast<std::int64_t>(domain.integerBound)) break;
    }
    return result;
}

EvalModIntervalCertificate certifyEvalModPolynomialIntervals(
    const EvalModPolynomial& polynomial, const EvalModDomain& domain,
    const std::string& radiusDecimal, std::size_t subdivisions, std::size_t precisionBits) {
    if (subdivisions < 2 || precisionBits < 128 || polynomial.basis != PolynomialBasis::Monomial)
        throw std::invalid_argument("invalid interval certificate input");
    const auto diagnostic = diagnoseEvalModPolynomialOnGrid(
        polynomial, domain, subdivisions + 1, radiusDecimal, 0.0, precisionBits,
        (domain.integerBound * 2 + 1) * (subdivisions + 1) * 5 + 1);
    const auto precision = static_cast<mpfr_prec_t>(precisionBits);
    Real radius(precision), rho(precision), maximumReal(precision), maximumModulus(precision);
    Real derivative(precision), power(precision), coefficient(precision), term(precision), one(precision);
    if (mpfr_set_str(radius.get(), radiusDecimal.c_str(), 10, MPFR_RNDU) != 0
        || mpfr_set_str(rho.get(), domain.normalizedResidualBoundDecimal.c_str(), 10, MPFR_RNDU) != 0)
        throw std::invalid_argument("invalid certificate domain decimal");
    mpfr_set_ui(maximumReal.get(), domain.integerBound, MPFR_RNDU);
    mpfr_add(maximumReal.get(), maximumReal.get(), rho.get(), MPFR_RNDU);
    mpfr_hypot(maximumModulus.get(), maximumReal.get(), radius.get(), MPFR_RNDU);
    mpfr_set_zero(derivative.get(), 0);
    mpfr_set_ui(power.get(), 1, MPFR_RNDU);
    for (std::size_t index = 1; index < polynomial.decimalCoefficients.size(); ++index) {
        if (index > 1) mpfr_mul(power.get(), power.get(), maximumModulus.get(), MPFR_RNDU);
        if (mpfr_set_str(coefficient.get(), polynomial.decimalCoefficients[index].c_str(), 10, MPFR_RNDU) != 0)
            throw std::invalid_argument("invalid certificate coefficient");
        mpfr_abs(coefficient.get(), coefficient.get(), MPFR_RNDU);
        mpfr_mul_ui(term.get(), coefficient.get(), index, MPFR_RNDU);
        mpfr_mul(term.get(), term.get(), power.get(), MPFR_RNDU);
        mpfr_add(derivative.get(), derivative.get(), term.get(), MPFR_RNDU);
    }
    mpfr_set_ui(one.get(), 1, MPFR_RNDU);
    mpfr_add(one.get(), one.get(), derivative.get(), MPFR_RNDU);
    Real realGap(precision), complexGap(precision), approximation(precision), complexError(precision);
    mpfr_div_ui(realGap.get(), rho.get(), subdivisions, MPFR_RNDU);
    mpfr_mul(realGap.get(), realGap.get(), one.get(), MPFR_RNDU);
    mpfr_set_d(approximation.get(), diagnostic.approximationMaxError, MPFR_RNDU);
    mpfr_add(approximation.get(), approximation.get(), realGap.get(), MPFR_RNDU);
    mpfr_max(complexGap.get(), rho.get(), radius.get(), MPFR_RNDU);
    mpfr_div_ui(complexGap.get(), complexGap.get(), subdivisions, MPFR_RNDU);
    mpfr_mul(complexGap.get(), complexGap.get(), one.get(), MPFR_RNDU);
    mpfr_set_d(complexError.get(), diagnostic.complexBoundaryErrorMax, MPFR_RNDU);
    mpfr_add(complexError.get(), complexError.get(), complexGap.get(), MPFR_RNDU);
    auto decimal = [](mpfr_srcptr value) {
        char* text = nullptr; mpfr_asprintf(&text, "%.RUe", value);
        std::string result = text ? text : ""; mpfr_free_str(text); return result;
    };
    return {decimal(approximation.get()), decimal(derivative.get()), decimal(complexError.get()),
            mpfr_get_d(approximation.get(), MPFR_RNDU), mpfr_get_d(derivative.get(), MPFR_RNDU),
            mpfr_get_d(complexError.get(), MPFR_RNDU), false};
}

} // namespace m2424::experimental
