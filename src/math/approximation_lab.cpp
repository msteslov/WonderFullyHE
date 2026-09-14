#include "m2424/experimental/evalmod_analysis/approximation_lab.hpp"

#include <mpfr.h>
#include "m2424/experimental/evalmod_analysis/evalround_execution.hpp"
#include "m2424/experimental/evalmod_analysis/exact_decimal.hpp"
#include "evalround_interval_internal.hpp"

#include <algorithm>
#include <array>
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

struct Interval {
    Real lo, hi;
    explicit Interval(mpfr_prec_t p) : lo(p), hi(p) {}
};
void multiplyInterval(const Interval& left, const Interval& right, Interval& out) {
    const auto precision=mpfr_get_prec(left.lo.get());
    std::array<Real, 4> productsLo{Real(precision), Real(precision), Real(precision), Real(precision)};
    std::array<Real, 4> productsHi{Real(precision), Real(precision), Real(precision), Real(precision)};
    mpfr_mul(productsLo[0].get(), left.lo.get(), right.lo.get(), MPFR_RNDD);
    mpfr_mul(productsLo[1].get(), left.lo.get(), right.hi.get(), MPFR_RNDD);
    mpfr_mul(productsLo[2].get(), left.hi.get(), right.lo.get(), MPFR_RNDD);
    mpfr_mul(productsLo[3].get(), left.hi.get(), right.hi.get(), MPFR_RNDD);
    mpfr_mul(productsHi[0].get(), left.lo.get(), right.lo.get(), MPFR_RNDU);
    mpfr_mul(productsHi[1].get(), left.lo.get(), right.hi.get(), MPFR_RNDU);
    mpfr_mul(productsHi[2].get(), left.hi.get(), right.lo.get(), MPFR_RNDU);
    mpfr_mul(productsHi[3].get(), left.hi.get(), right.hi.get(), MPFR_RNDU);
    mpfr_set(out.lo.get(), productsLo[0].get(), MPFR_RNDD);
    mpfr_set(out.hi.get(), productsHi[0].get(), MPFR_RNDU);
    for (int i = 1; i < 4; ++i) {
        mpfr_min(out.lo.get(), out.lo.get(), productsLo[i].get(), MPFR_RNDD);
        mpfr_max(out.hi.get(), out.hi.get(), productsHi[i].get(), MPFR_RNDU);
    }
}
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

std::vector<mpq_class> shiftExactCoefficients(
    const std::vector<mpq_class>& coefficients, std::int64_t center) {
    std::vector<mpq_class> shifted{coefficients.back()};
    const mpq_class exactCenter(static_cast<long>(center));
    for (std::size_t index = coefficients.size() - 1; index-- > 0;) {
        std::vector<mpq_class> next(shifted.size() + 1);
        for (std::size_t power = 0; power < shifted.size(); ++power) {
            next[power] += exactCenter * shifted[power];
            next[power + 1] += shifted[power];
        }
        next[0] += coefficients[index];
        shifted = std::move(next);
    }
    for (auto& coefficient : shifted) coefficient.canonicalize();
    return shifted;
}

} // namespace

std::vector<mpq_class> detail::shiftPolynomialToIntegerCenterExact(
    const EvalModPolynomial& polynomial, std::int64_t center) {
    if (polynomial.basis != PolynomialBasis::Monomial
        || polynomial.decimalCoefficients.empty()
        || polynomial.decimalCoefficients.size() > 257)
        throw std::invalid_argument("invalid exact polynomial shift input");
    std::vector<mpq_class> coefficients;
    coefficients.reserve(polynomial.decimalCoefficients.size());
    for (const auto& text : polynomial.decimalCoefficients)
        coefficients.push_back(parseExactDecimal(text));
    return shiftExactCoefficients(coefficients, center);
}

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
    if (subdivisions < 2 || precisionBits < 128 || polynomial.basis != PolynomialBasis::Monomial
        || polynomial.decimalCoefficients.empty() || !isEvalModDomainValid(domain)
        || domain.integerBound > static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max()))
        throw std::invalid_argument("invalid interval certificate input");
    const auto precision = static_cast<mpfr_prec_t>(precisionBits);
    Real radius(precision), rho(precision), maximumReal(precision), maximumModulus(precision);
    Real derivative(precision), power(precision), coefficient(precision), term(precision);
    Real polynomialBound(precision);
    if (mpfr_set_str(radius.get(), radiusDecimal.c_str(), 10, MPFR_RNDU) != 0
        || mpfr_set_str(rho.get(), domain.normalizedResidualBoundDecimal.c_str(), 10, MPFR_RNDU) != 0)
        throw std::invalid_argument("invalid certificate domain decimal");
    mpfr_set_ui(maximumReal.get(), domain.integerBound, MPFR_RNDU);
    mpfr_add(maximumReal.get(), maximumReal.get(), rho.get(), MPFR_RNDU);
    mpfr_hypot(maximumModulus.get(), maximumReal.get(), radius.get(), MPFR_RNDU);
    mpfr_set_zero(derivative.get(), 0);
    mpfr_set_zero(polynomialBound.get(), 0);
    for (std::size_t index = polynomial.decimalCoefficients.size(); index-- > 0;) {
        mpfr_mul(polynomialBound.get(), polynomialBound.get(), maximumModulus.get(), MPFR_RNDU);
        if (mpfr_set_str(coefficient.get(), polynomial.decimalCoefficients[index].c_str(), 10,
                         MPFR_RNDA) != 0)
            throw std::invalid_argument("invalid certificate coefficient");
        mpfr_abs(coefficient.get(), coefficient.get(), MPFR_RNDA);
        mpfr_add(polynomialBound.get(), polynomialBound.get(), coefficient.get(), MPFR_RNDU);
    }
    mpfr_set_ui(power.get(), 1, MPFR_RNDU);
    for (std::size_t index = 1; index < polynomial.decimalCoefficients.size(); ++index) {
        if (index > 1) mpfr_mul(power.get(), power.get(), maximumModulus.get(), MPFR_RNDU);
        if (mpfr_set_str(coefficient.get(), polynomial.decimalCoefficients[index].c_str(), 10,
                         MPFR_RNDA) != 0)
            throw std::invalid_argument("invalid certificate coefficient");
        mpfr_abs(coefficient.get(), coefficient.get(), MPFR_RNDA);
        mpfr_mul_ui(term.get(), coefficient.get(), index, MPFR_RNDU);
        mpfr_mul(term.get(), term.get(), power.get(), MPFR_RNDU);
        mpfr_add(derivative.get(), derivative.get(), term.get(), MPFR_RNDU);
    }
    std::vector<Interval> coefficientIntervals;
    coefficientIntervals.reserve(polynomial.decimalCoefficients.size());
    for (const auto& decimalCoefficient : polynomial.decimalCoefficients) {
        coefficientIntervals.emplace_back(precision);
        if (mpfr_set_str(coefficientIntervals.back().lo.get(), decimalCoefficient.c_str(), 10,
                         MPFR_RNDD) != 0
            || mpfr_set_str(coefficientIntervals.back().hi.get(), decimalCoefficient.c_str(), 10,
                            MPFR_RNDU) != 0)
            throw std::invalid_argument("invalid interval coefficient");
    }

    Real approximation(precision), complexError(precision), targetModulus(precision);
    mpfr_set_zero(approximation.get(), 0);
    for (std::int64_t integer = -static_cast<std::int64_t>(domain.integerBound);
         integer <= static_cast<std::int64_t>(domain.integerBound); ++integer) {
        for (std::size_t cell = 0; cell < subdivisions; ++cell) {
            Interval x(precision), valueInterval(precision), productInterval(precision);
            Real fraction(precision), scaled(precision);
            mpfr_set_ui(fraction.get(), 2 * cell + 1, MPFR_RNDD);
            mpfr_div_ui(fraction.get(), fraction.get(), subdivisions, MPFR_RNDD);
            mpfr_sub_ui(fraction.get(), fraction.get(), 1, MPFR_RNDD);
            mpfr_mul(scaled.get(), rho.get(), fraction.get(), MPFR_RNDD);
            mpfr_add_si(x.lo.get(), scaled.get(), integer, MPFR_RNDD);
            mpfr_set_ui(fraction.get(), 2 * cell + 1, MPFR_RNDU);
            mpfr_div_ui(fraction.get(), fraction.get(), subdivisions, MPFR_RNDU);
            mpfr_sub_ui(fraction.get(), fraction.get(), 1, MPFR_RNDU);
            mpfr_mul(scaled.get(), rho.get(), fraction.get(), MPFR_RNDU);
            mpfr_add_si(x.hi.get(), scaled.get(), integer, MPFR_RNDU);
            mpfr_set_zero(valueInterval.lo.get(), 0);
            mpfr_set_zero(valueInterval.hi.get(), 0);
            for (std::size_t index = coefficientIntervals.size(); index-- > 0;) {
                multiplyInterval(valueInterval, x, productInterval);
                mpfr_add(valueInterval.lo.get(), productInterval.lo.get(),
                         coefficientIntervals[index].lo.get(), MPFR_RNDD);
                mpfr_add(valueInterval.hi.get(), productInterval.hi.get(),
                         coefficientIntervals[index].hi.get(), MPFR_RNDU);
            }
            mpfr_sub(valueInterval.lo.get(), valueInterval.lo.get(), x.hi.get(), MPFR_RNDD);
            mpfr_add_si(valueInterval.lo.get(), valueInterval.lo.get(), integer, MPFR_RNDD);
            mpfr_sub(valueInterval.hi.get(), valueInterval.hi.get(), x.lo.get(), MPFR_RNDU);
            mpfr_add_si(valueInterval.hi.get(), valueInterval.hi.get(), integer, MPFR_RNDU);
            Real endpoint(precision);
            mpfr_abs(endpoint.get(), valueInterval.lo.get(), MPFR_RNDU);
            mpfr_max(approximation.get(), approximation.get(), endpoint.get(), MPFR_RNDU);
            mpfr_abs(endpoint.get(), valueInterval.hi.get(), MPFR_RNDU);
            mpfr_max(approximation.get(), approximation.get(), endpoint.get(), MPFR_RNDU);
        }
    }
    Real halfCell(precision), lipschitz(precision);
    mpfr_div_ui(halfCell.get(), rho.get(), subdivisions, MPFR_RNDU);
    mpfr_add_ui(lipschitz.get(), derivative.get(), 1, MPFR_RNDU);
    mpfr_mul(halfCell.get(), halfCell.get(), lipschitz.get(), MPFR_RNDU);
    mpfr_add(approximation.get(), approximation.get(), halfCell.get(), MPFR_RNDU);
    mpfr_hypot(targetModulus.get(), rho.get(), radius.get(), MPFR_RNDU);
    mpfr_add(complexError.get(), polynomialBound.get(), targetModulus.get(), MPFR_RNDU);
    auto decimal = [](mpfr_srcptr value) {
        char* text = nullptr; mpfr_asprintf(&text, "%.RUe", value);
        std::string result = text ? text : ""; mpfr_free_str(text); return result;
    };
    return {decimal(approximation.get()), decimal(derivative.get()), decimal(complexError.get()),
            mpfr_get_d(approximation.get(), MPFR_RNDU), mpfr_get_d(derivative.get(), MPFR_RNDU),
            mpfr_get_d(complexError.get(), MPFR_RNDU),
            mpfr_number_p(approximation.get()) && mpfr_number_p(derivative.get())
                && mpfr_number_p(complexError.get())};
}

EvalRoundDigitPolynomial certifyEvalRoundDigitPolynomial(const EvalRoundProblem& problem,std::size_t digit,
    const EvalModPolynomial& polynomial,std::size_t subdivisions) {
    if(problem.K>4096||!std::isfinite(problem.rho)||problem.rho<0||problem.rho>=.5||
       digit>=evalRoundDigitCount(problem.K,EvalRoundRadix::Binary)||subdivisions<2||subdivisions>4096||
       polynomial.basis!=PolynomialBasis::Monomial||polynomial.decimalCoefficients.empty()||polynomial.decimalCoefficients.size()>257)
        throw std::invalid_argument("Invalid binary digit interval proof domain/polynomial/work limit");
    constexpr mpfr_prec_t precision=384;
    std::vector<mpq_class> exactCoefficients;
    std::vector<Interval> coefficients;
    for(const auto& text:polynomial.decimalCoefficients) {
        exactCoefficients.push_back(parseExactDecimal(text)); coefficients.emplace_back(precision);
        mpfr_set_q(coefficients.back().lo.get(),exactCoefficients.back().get_mpq_t(),MPFR_RNDD);
        mpfr_set_q(coefficients.back().hi.get(),exactCoefficients.back().get_mpq_t(),MPFR_RNDU);
    }
    const mpq_class exactRho(problem.rho);
    Real directMaximum(precision),centeredMaximum(precision),absolute(precision);
    mpfr_set_zero(directMaximum.get(),0);mpfr_set_zero(centeredMaximum.get(),0);
    for(std::int64_t I=-std::int64_t(problem.K);I<=std::int64_t(problem.K);++I) {
        const auto target=evalRoundIntegerDigits(I,problem.K,EvalRoundRadix::Binary)[digit];
        auto shifted=shiftExactCoefficients(exactCoefficients,I);
        shifted[0]-=target;shifted[0].canonicalize();
        std::vector<Interval> shiftedIntervals;
        shiftedIntervals.reserve(shifted.size());
        for(const auto& exact:shifted) {
            shiftedIntervals.emplace_back(precision);
            mpfr_set_q(shiftedIntervals.back().lo.get(),exact.get_mpq_t(),MPFR_RNDD);
            mpfr_set_q(shiftedIntervals.back().hi.get(),exact.get_mpq_t(),MPFR_RNDU);
        }
        for(std::size_t cell=0;cell<subdivisions;++cell) {
            // Exact dyadic rho and complete closed cells; no point/grid inference.
            mpq_class leftFraction(2*cell,subdivisions),rightFraction(2*(cell+1),subdivisions);
            leftFraction.canonicalize();rightFraction.canonicalize();
            const mpq_class yLo=exactRho*(leftFraction-1);
            const mpq_class yHi=exactRho*(rightFraction-1);
            const mpq_class lo=mpq_class(static_cast<long>(I))+yLo;
            const mpq_class hi=mpq_class(static_cast<long>(I))+yHi;
            Interval x(precision),value(precision),product(precision);
            mpfr_set_q(x.lo.get(),lo.get_mpq_t(),MPFR_RNDD);mpfr_set_q(x.hi.get(),hi.get_mpq_t(),MPFR_RNDU);
            mpfr_set_zero(value.lo.get(),0);mpfr_set_zero(value.hi.get(),0);
            for(std::size_t k=coefficients.size();k-->0;) {
                multiplyInterval(value,x,product);
                mpfr_add(value.lo.get(),product.lo.get(),coefficients[k].lo.get(),MPFR_RNDD);
                mpfr_add(value.hi.get(),product.hi.get(),coefficients[k].hi.get(),MPFR_RNDU);
            }
            mpfr_sub_si(value.lo.get(),value.lo.get(),target,MPFR_RNDD);
            mpfr_sub_si(value.hi.get(),value.hi.get(),target,MPFR_RNDU);
            for(auto endpoint:{value.lo.get(),value.hi.get()}) {
                mpfr_abs(absolute.get(),endpoint,MPFR_RNDU);mpfr_max(directMaximum.get(),directMaximum.get(),absolute.get(),MPFR_RNDU);
            }
            Interval y(precision),centeredValue(precision),centeredProduct(precision);
            mpfr_set_q(y.lo.get(),yLo.get_mpq_t(),MPFR_RNDD);mpfr_set_q(y.hi.get(),yHi.get_mpq_t(),MPFR_RNDU);
            mpfr_set_zero(centeredValue.lo.get(),0);mpfr_set_zero(centeredValue.hi.get(),0);
            for(std::size_t k=shiftedIntervals.size();k-->0;) {
                multiplyInterval(centeredValue,y,centeredProduct);
                mpfr_add(centeredValue.lo.get(),centeredProduct.lo.get(),shiftedIntervals[k].lo.get(),MPFR_RNDD);
                mpfr_add(centeredValue.hi.get(),centeredProduct.hi.get(),shiftedIntervals[k].hi.get(),MPFR_RNDU);
            }
            for(auto endpoint:{centeredValue.lo.get(),centeredValue.hi.get()}) {
                mpfr_abs(absolute.get(),endpoint,MPFR_RNDU);mpfr_max(centeredMaximum.get(),centeredMaximum.get(),absolute.get(),MPFR_RNDU);
            }
        }
    }
    EvalRoundDigitPolynomial out;out.digitIndex=digit;out.polynomial=polynomial;
    out.certifiedK=problem.K;out.certifiedRho=problem.rho;
    out.proof=EvalRoundPolynomialProof::OutwardInterval;out.intervalSubdivisions=subdivisions;
    out.intervalProofPrecisionBits=precision;
    const std::string directProvenance="384-bit outward MPFR direct-x Horner on every complete closed subcell of D_K,rho";
    const std::string centeredProvenance="exact rational Taylor shift x=I+y followed by 384-bit outward MPFR Horner on a complete closed partition of [-rho,+rho]";
    out.directXApproximationError={mpfr_get_d(directMaximum.get(),MPFR_RNDU),BootstrapBoundKind::Deterministic,directProvenance,{}};
    out.centeredApproximationError={mpfr_get_d(centeredMaximum.get(),MPFR_RNDU),BootstrapBoundKind::Deterministic,centeredProvenance,{}};
    const bool useCentered=mpfr_less_p(centeredMaximum.get(),directMaximum.get());
    out.selectedIntervalProofMethod=useCentered?EvalRoundIntervalProofMethod::CenteredShiftHorner:EvalRoundIntervalProofMethod::DirectXHorner;
    out.provenance=std::string("minimum of two independent whole-domain bounds for the same exact polynomial and bit_j(I+K): ")
        +(useCentered?centeredProvenance:directProvenance)+"; grid excluded";
    out.approximationError=useCentered?out.centeredApproximationError:out.directXApproximationError;
    out.approximationError.provenance=out.provenance;
    out.verified=mpfr_number_p(directMaximum.get())&&mpfr_number_p(centeredMaximum.get());
    return out;
}

} // namespace m2424::experimental
