#include "m2424/experimental/evalmod_analysis/domain_analysis.hpp"

#include <mpfr.h>
#include <gmpxx.h>

#include <cmath>
#include <cstdlib>
#include <limits>
#include <stdexcept>

namespace m2424::experimental {
namespace {

class Real {
public:
    explicit Real(mpfr_prec_t precision) { mpfr_init2(value_, precision); }
    ~Real() { mpfr_clear(value_); }
    mpfr_ptr get() { return value_; }
    mpfr_srcptr get() const { return value_; }
private:
    mpfr_t value_;
};

bool finiteNonNegative(double value) { return std::isfinite(value) && value >= 0.0; }

std::string decimalUp(mpfr_srcptr value, std::size_t precisionBits) {
    char* text = nullptr;
    const int digits = static_cast<int>(std::ceil(precisionBits * 0.30103)) + 3;
    mpfr_asprintf(&text, "%.*RUg", digits, value);
    std::string result = text ? text : "";
    mpfr_free_str(text);
    return result;
}

} // namespace

EvalModDomain estimateEvalModDomain(const EvalModCiphertextModel& model) {
    if (model.coefficientCount == 0 || !finiteNonNegative(model.deterministicIntegerOffset)
        || !finiteNonNegative(model.integerNoiseSubgaussianSigma)
        || !finiteNonNegative(model.normalizedMessageAbsBound)
        || !finiteNonNegative(model.normalizedEncodingErrorAbsBound)
        || !finiteNonNegative(model.normalizedCoeffToSlotErrorAbsBound)
        || !finiteNonNegative(model.relativePeriodMismatchAbsBound)
        || !finiteNonNegative(model.additiveNormalizationErrorAbsBound)
        || !std::isfinite(model.failureProbabilityLog2) || model.failureProbabilityLog2 >= 0.0
        || model.analysisPrecisionBits < 64
        || model.analysisPrecisionBits > static_cast<std::size_t>(std::numeric_limits<int>::max() / 2)
        || model.provenance.derivation.empty() || model.provenance.includedNoiseSources.empty()
        || model.provenance.assumptions.empty()
        || (model.tailModel == TailModel::Deterministic
            && model.integerNoiseSubgaussianSigma != 0.0)) {
        throw std::invalid_argument("invalid EvalMod ciphertext model");
    }

    const auto precision = static_cast<mpfr_prec_t>(model.analysisPrecisionBits);
    Real tail(precision), work(precision), logarithm(precision), sigma(precision), offset(precision);
    mpfr_set_d(offset.get(), model.deterministicIntegerOffset, MPFR_RNDU);
    if (model.tailModel == TailModel::Deterministic) {
        mpfr_set_zero(tail.get(), 0);
    } else {
        mpfr_set_ui(work.get(), model.coefficientCount, MPFR_RNDU);
        mpfr_mul_ui(work.get(), work.get(), 2, MPFR_RNDU);
        mpfr_log(logarithm.get(), work.get(), MPFR_RNDU);
        mpfr_const_log2(work.get(), MPFR_RNDU);
        mpfr_mul_d(work.get(), work.get(), -model.failureProbabilityLog2, MPFR_RNDU);
        mpfr_add(logarithm.get(), logarithm.get(), work.get(), MPFR_RNDU);
        mpfr_mul_ui(logarithm.get(), logarithm.get(), 2, MPFR_RNDU);
        mpfr_sqrt(logarithm.get(), logarithm.get(), MPFR_RNDU);
        mpfr_set_d(sigma.get(), model.integerNoiseSubgaussianSigma, MPFR_RNDU);
        mpfr_mul(tail.get(), sigma.get(), logarithm.get(), MPFR_RNDU);
    }
    mpfr_add(work.get(), offset.get(), tail.get(), MPFR_RNDU);
    mpfr_ceil(work.get(), work.get());
    mpz_class sizeMaximum;
    const std::size_t nativeMaximum = std::numeric_limits<std::size_t>::max();
    mpz_import(sizeMaximum.get_mpz_t(), 1, 1, sizeof(nativeMaximum), 0, 0, &nativeMaximum);
    Real sizeMaximumMpfr(precision);
    mpfr_set_z(sizeMaximumMpfr.get(), sizeMaximum.get_mpz_t(), MPFR_RNDN);
    if (mpfr_cmp(work.get(), sizeMaximumMpfr.get()) > 0) {
        throw std::overflow_error("EvalMod integer bound does not fit size_t");
    }
    mpz_class integerBoundExact;
    mpfr_get_z(integerBoundExact.get_mpz_t(), work.get(), MPFR_RNDU);
    std::size_t integerBound{};
    std::size_t imported{};
    mpz_export(&integerBound, &imported, 1, sizeof(integerBound), 0, 0,
               integerBoundExact.get_mpz_t());

    Real residual(precision), term(precision);
    mpfr_set_d(residual.get(), model.normalizedMessageAbsBound, MPFR_RNDU);
    mpfr_set_d(term.get(), model.normalizedEncodingErrorAbsBound, MPFR_RNDU);
    mpfr_add(residual.get(), residual.get(), term.get(), MPFR_RNDU);
    mpfr_set_d(term.get(), model.normalizedCoeffToSlotErrorAbsBound, MPFR_RNDU);
    mpfr_add(residual.get(), residual.get(), term.get(), MPFR_RNDU);
    mpfr_set_ui(term.get(), integerBound, MPFR_RNDU);
    mpfr_add(term.get(), term.get(), residual.get(), MPFR_RNDU);
    Real mismatch(precision);
    mpfr_set_d(mismatch.get(), model.relativePeriodMismatchAbsBound, MPFR_RNDU);
    mpfr_mul(term.get(), term.get(), mismatch.get(), MPFR_RNDU);
    mpfr_add(residual.get(), residual.get(), term.get(), MPFR_RNDU);
    mpfr_set_d(term.get(), model.additiveNormalizationErrorAbsBound, MPFR_RNDU);
    mpfr_add(residual.get(), residual.get(), term.get(), MPFR_RNDU);

    const double rho = mpfr_get_d(residual.get(), MPFR_RNDU);
    EvalModDomain result{integerBound, rho, decimalUp(residual.get(), model.analysisPrecisionBits),
                         model.failureProbabilityLog2};
    if (!isEvalModDomainValid(result)) {
        throw std::domain_error("EvalMod residual reaches a rounding discontinuity");
    }
    return result;
}

double EvalModDomain::discontinuityMargin() const {
    const double raw = 0.5 - normalizedResidualBound;
    return std::nextafter(raw, -std::numeric_limits<double>::infinity());
}

bool isEvalModDomainValid(const EvalModDomain& domain) {
    if (!std::isfinite(domain.normalizedResidualBound)
        || domain.normalizedResidualBound < 0.0 || domain.normalizedResidualBound >= 0.5
        || domain.normalizedResidualBoundDecimal.empty()
        || !std::isfinite(domain.failureProbabilityLog2) || domain.failureProbabilityLog2 >= 0.0) {
        return false;
    }
    Real exact(128);
    if (mpfr_set_str(exact.get(), domain.normalizedResidualBoundDecimal.c_str(), 10, MPFR_RNDN) != 0
        || !mpfr_number_p(exact.get()) || mpfr_sgn(exact.get()) < 0
        || mpfr_cmp_d(exact.get(), domain.normalizedResidualBound) > 0) {
        return false;
    }
    const double previous = std::nextafter(domain.normalizedResidualBound,
                                           -std::numeric_limits<double>::infinity());
    return mpfr_cmp_d(exact.get(), previous) > 0 && domain.discontinuityMargin() > 0.0;
}

} // namespace m2424::experimental
