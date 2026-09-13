#include "m2424/experimental/evalmod_analysis/evalround_reference.hpp"
#include <mpfr.h>
#include "m2424/experimental/evalmod_analysis/exact_decimal.hpp"
#include <cmath>
#include <stdexcept>

namespace m2424::experimental {
namespace {
void validPrecision(std::size_t p) {
    if (p < 64 || p > 16384) throw std::invalid_argument("EvalRound reference precision must be 64..16384 bits");
}
struct Real {
    mpfr_t value;
    explicit Real(std::size_t p) { validPrecision(p); mpfr_init2(value, static_cast<mpfr_prec_t>(p)); }
    ~Real() { mpfr_clear(value); }
    Real(const Real&) = delete;
    Real& operator=(const Real&) = delete;
};
mpq_class asRational(const Real& x) { mpq_class q; mpfr_get_q(q.get_mpq_t(), x.value); return q; }
mpq_class rounded(const mpq_class& q, std::size_t p) {
    Real x(p); mpfr_set_q(x.value, q.get_mpq_t(), MPFR_RNDN); return asRational(x);
}
EvalRoundExactComplex plus(const EvalRoundExactComplex& a, const EvalRoundExactComplex& b) {
    return {a.real + b.real, a.imag + b.imag};
}
EvalRoundExactComplex times(const EvalRoundExactComplex& a, const EvalRoundExactComplex& b) {
    return {a.real*b.real - a.imag*b.imag, a.real*b.imag + a.imag*b.real};
}
EvalRoundExactComplex scaled(const EvalRoundExactComplex& x, const mpq_class& factor) {
    return {x.real*factor, x.imag*factor};
}
}
ExactInteger evalRoundReferenceCenter(const EvalRoundProblem& p, const mpq_class& z) {
    if (!std::isfinite(p.rho) || p.rho < 0 || p.rho >= 0.5) throw std::invalid_argument("invalid EvalRound reference domain");
    mpq_class shifted = z + mpq_class(1, 2);
    ExactInteger center;
    mpz_fdiv_q(center.get_mpz_t(), shifted.get_num_mpz_t(), shifted.get_den_mpz_t());
    mpq_class rho = 0;
    if (p.rho != 0) {
        const auto exact = ExactScale::fromBinaryDouble(p.rho);
        rho = mpq_class(exact.numerator, exact.denominator);
    }
    const mpq_class difference = z - center;
    const ExactInteger K(std::to_string(p.K));
    if (center < -K || center > K || abs(difference) > rho)
        throw std::invalid_argument("reference input outside union of EvalRound intervals");
    return center;
}
std::vector<int> extractEvalRoundReferenceDigits(const EvalRoundProblem& p, EvalRoundRadix radix, const mpq_class& z) {
    const auto center = evalRoundReferenceCenter(p, z);
    return evalRoundIntegerDigits(std::stoll(center.get_str()), p.K, radix);
}
ExactInteger reconstructEvalRoundDigitsExact(const std::vector<int>& digits, std::uint32_t K, EvalRoundRadix radix) {
    // Validate the alphabet, length and range independently of the GMP sum below.
    (void)reconstructEvalRoundInteger(digits, K, radix);
    ExactInteger result = radix == EvalRoundRadix::Binary ? -ExactInteger(std::to_string(K)) : ExactInteger(0);
    ExactInteger weight = 1;
    for (int digit : digits) { result += weight * digit; weight *= static_cast<unsigned>(radix); }
    return result;
}
mpq_class evalRoundBinaryCleanerExact(const mpq_class& x) { return 3*x*x - 2*x*x*x; }
EvalRoundExactComplex evalRoundTernaryCleanerExact(const EvalRoundExactComplex& x) {
    const EvalRoundExactComplex y{x.real, -x.imag};
    return scaled(plus(plus(times(y,y), scaled(x,4)), scaled(times(times(x,x),y),-2)), mpq_class(1,3));
}
EvalRoundExactComplex evalRoundRootReference(int trit, std::size_t precisionBits) {
    validPrecision(precisionBits);
    if (trit < -1 || trit > 1) throw std::invalid_argument("invalid trit");
    if (!trit) return {1,0};
    Real imaginary(precisionBits);
    mpfr_sqrt_ui(imaginary.value, 3, MPFR_RNDN);
    mpfr_div_ui(imaginary.value, imaginary.value, 2, MPFR_RNDN);
    if (trit < 0) mpfr_neg(imaginary.value, imaginary.value, MPFR_RNDN);
    return {mpq_class(-1,2), asRational(imaginary)};
}
EvalRoundReferenceTrace evaluateEvalRoundReference(const EvalRoundProblem& p, const EvalRoundCandidate& candidate,
    const mpq_class& z, const std::vector<std::size_t>& cleaningCounts, std::size_t precisionBits) {
    validPrecision(precisionBits);
    EvalRoundReferenceTrace trace;
    trace.center = evalRoundReferenceCenter(p,z);
    trace.targetDigits = extractEvalRoundReferenceDigits(p,candidate.radix,z);
    trace.precisionBits = precisionBits;
    if (cleaningCounts.size() != trace.targetDigits.size()) throw std::invalid_argument("reference cleaning count mismatch");
    for (const auto rounds : cleaningCounts) if (rounds > 64) throw std::invalid_argument("reference cleaning limit exceeded");
    std::vector<EvalRoundExactComplex> extracted;
    switch (candidate.extraction.method) {
    case EvalRoundExtractionMethod::PiecewiseReference:
        for (int digit : trace.targetDigits) extracted.push_back(candidate.radix == EvalRoundRadix::Binary
            ? EvalRoundExactComplex{digit,0} : evalRoundRootReference(digit,precisionBits));
        break;
    case EvalRoundExtractionMethod::BinaryQuadraticK1:
        if (p.K != 1 || candidate.radix != EvalRoundRadix::Binary) throw std::invalid_argument("quadratic reference requires binary K=1");
        extracted = {{1-z*z,0},{(z*z+z)/2,0}};
        break;
    case EvalRoundExtractionMethod::TernaryPhaseReferenceK1: {
        if (p.K != 1 || candidate.radix != EvalRoundRadix::BalancedTernary) throw std::invalid_argument("phase reference requires ternary K=1");
        Real angle(precisionBits), input(precisionBits), real(precisionBits), imaginary(precisionBits);
        mpfr_const_pi(angle.value,MPFR_RNDN);
        mpfr_set_q(input.value,z.get_mpq_t(),MPFR_RNDN);
        mpfr_mul(angle.value,angle.value,input.value,MPFR_RNDN);
        mpfr_mul_ui(angle.value,angle.value,2,MPFR_RNDN);
        mpfr_div_ui(angle.value,angle.value,3,MPFR_RNDN);
        mpfr_sin_cos(imaginary.value,real.value,angle.value,MPFR_RNDN);
        extracted = {{asRational(real),asRational(imaginary)}};
        break;
    }
    case EvalRoundExtractionMethod::ExternalPolynomial:
    case EvalRoundExtractionMethod::DigitExtract:
        if(candidate.radix!=EvalRoundRadix::Binary||candidate.extraction.polynomials.size()!=trace.targetDigits.size())
            throw std::invalid_argument("Concrete binary polynomials required for reference evaluation");
        for(const auto& digit:candidate.extraction.polynomials) {
            if(digit.polynomial.basis!=PolynomialBasis::Monomial||digit.polynomial.decimalCoefficients.empty())
                throw std::invalid_argument("Known monomial coefficients required");
            mpq_class value=0;
            for(auto it=digit.polynomial.decimalCoefficients.rbegin();it!=digit.polynomial.decimalCoefficients.rend();++it)
                value=value*z+parseExactDecimal(*it);
            extracted.push_back({value,0});
        }
        break;
    default: throw std::invalid_argument("Unknown extraction method");
    }
    for (std::size_t j=0; j<extracted.size(); ++j) {
        auto value = extracted[j];
        value = {rounded(value.real,precisionBits),rounded(value.imag,precisionBits)};
        std::vector<EvalRoundExactComplex> stages{value};
        for (std::size_t r=0; r<cleaningCounts[j]; ++r) {
            value = candidate.radix == EvalRoundRadix::Binary ? EvalRoundExactComplex{evalRoundBinaryCleanerExact(value.real),0}
                : evalRoundTernaryCleanerExact(value);
            value = {rounded(value.real,precisionBits),rounded(value.imag,precisionBits)};
            stages.push_back(value);
        }
        trace.digitStages.push_back(std::move(stages));
    }
    Real gain(precisionBits);
    mpfr_sqrt_ui(gain.value,3,MPFR_RNDN);
    mpfr_ui_div(gain.value,2,gain.value,MPFR_RNDN);
    const mpq_class ternaryGain = asRational(gain);
    trace.reconstructed = candidate.radix == EvalRoundRadix::Binary ? -mpq_class(ExactInteger(std::to_string(p.K))) : mpq_class(0);
    ExactInteger weight = 1;
    for (const auto& stages : trace.digitStages) {
        const auto& value = stages.back();
        trace.reconstructed += weight * (candidate.radix == EvalRoundRadix::Binary ? value.real : ternaryGain*value.imag);
        weight *= static_cast<unsigned>(candidate.radix);
    }
    trace.reconstructed = rounded(trace.reconstructed,precisionBits);
    return trace;
}
} // namespace m2424::experimental
