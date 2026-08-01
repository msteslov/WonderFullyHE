#include "m2424/experimental/evalmod_analysis/exact_modular_oracle.hpp"

#include <mpfr.h>

#include <cstdlib>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <utility>

namespace m2424::experimental {
namespace {

ExactInteger extendedGcd(ExactInteger a, ExactInteger b, ExactInteger& x, ExactInteger& y) {
    if (b == 0) { x = 1; y = 0; return a; }
    ExactInteger x1, y1;
    const ExactInteger gcd = extendedGcd(b, a % b, x1, y1);
    x = y1;
    y = x1 - (a / b) * y1;
    return gcd;
}

ExactInteger inverseMod(const ExactInteger& value, const ExactInteger& modulus) {
    ExactInteger x, y;
    const ExactInteger gcd = extendedGcd(value % modulus, modulus, x, y);
    if (gcd != 1) throw std::invalid_argument("CRT moduli must be pairwise coprime");
    x %= modulus;
    return x < 0 ? x + modulus : x;
}

ExactInteger fromUint64(std::uint64_t value) {
    ExactInteger result;
    mpz_import(result.get_mpz_t(), 1, 1, sizeof(value), 0, 0, &value);
    return result;
}

} // namespace

ExactCoefficientOracleResult exactCoefficientOracle(const std::vector<std::uint64_t>& residues,
                                                     const std::vector<std::uint64_t>& moduli,
                                                     const ExactInteger& qSource,
                                                     const ExactScale& outputScale,
                                                     std::size_t precisionBits) {
    if (residues.empty() || residues.size() != moduli.size() || qSource <= 1
        || outputScale.numerator <= 0 || outputScale.denominator <= 0 || precisionBits < 64
        || precisionBits > static_cast<std::size_t>(std::numeric_limits<int>::max() / 2)
        || precisionBits > static_cast<std::size_t>(std::numeric_limits<mpfr_prec_t>::max())) {
        throw std::invalid_argument("invalid exact coefficient oracle input");
    }
    ExactInteger reconstructed = 0;
    ExactInteger product = 1;
    for (std::size_t i = 0; i < residues.size(); ++i) {
        if (moduli[i] < 2 || residues[i] >= moduli[i]) {
            throw std::invalid_argument("invalid RNS residue or modulus");
        }
        const ExactInteger modulus = fromUint64(moduli[i]);
        const ExactInteger delta = (fromUint64(residues[i]) - reconstructed) % modulus;
        const ExactInteger correction = ((delta < 0 ? delta + modulus : delta)
            * inverseMod(product % modulus, modulus)) % modulus;
        reconstructed += product * correction;
        product *= modulus;
    }
    if (qSource > product || product % qSource != 0) {
        throw std::invalid_argument("qSource must divide the full CRT modulus");
    }

    ExactInteger centered = reconstructed % qSource;
    if (centered < 0) centered += qSource;
    if (2 * centered >= qSource) centered -= qSource;

    mpfr_t expected, error;
    const auto precision = static_cast<mpfr_prec_t>(precisionBits);
    mpfr_inits2(precision, expected, error, static_cast<mpfr_ptr>(nullptr));
    mpq_class exactValue(centered * outputScale.denominator, outputScale.numerator);
    exactValue.canonicalize();
    const int inexact = mpfr_set_q(expected, exactValue.get_mpq_t(), MPFR_RNDN);
    char* digits = nullptr;
    mpfr_asprintf(&digits, "%.*Rg", static_cast<int>(std::ceil(precisionBits * 0.30103)) + 2, expected);
    std::string decimal = digits ? digits : "";
    mpfr_free_str(digits);
    if (inexact == 0) {
        mpfr_set_zero(error, 0);
    } else {
        mpfr_set_ui_2exp(error, 1, mpfr_get_exp(expected) - precision, MPFR_RNDU);
    }
    digits = nullptr;
    mpfr_asprintf(&digits, "%.Re", error);
    std::string errorDecimal = digits ? digits : "";
    mpfr_free_str(digits);
    mpfr_clears(expected, error, static_cast<mpfr_ptr>(nullptr));
    return {reconstructed, product, centered, outputScale.numerator, outputScale.denominator,
            std::move(decimal), std::move(errorDecimal), precisionBits};
}

ExactScale ExactScale::rational(ExactInteger numerator, ExactInteger denominator) {
    if (numerator <= 0 || denominator <= 0) throw std::invalid_argument("invalid exact scale");
    ExactInteger divisor;
    mpz_gcd(divisor.get_mpz_t(), numerator.get_mpz_t(), denominator.get_mpz_t());
    return {numerator / divisor, denominator / divisor};
}

ExactScale ExactScale::fromBinaryDouble(double value) {
    if (!std::isfinite(value) || value <= 0.0) throw std::invalid_argument("invalid binary scale");
    std::uint64_t bits{};
    static_assert(sizeof(bits) == sizeof(value), "unexpected double size");
    std::memcpy(&bits, &value, sizeof(bits));
    const std::uint64_t fraction = bits & ((std::uint64_t{1} << 52) - 1);
    const unsigned exponentBits = static_cast<unsigned>((bits >> 52) & 0x7ffU);
    ExactInteger mantissa;
    int exponent;
    if (exponentBits == 0) {
        mantissa = fromUint64(fraction);
        exponent = -1074;
    } else {
        mantissa = fromUint64((std::uint64_t{1} << 52) | fraction);
        exponent = static_cast<int>(exponentBits) - 1023 - 52;
    }
    ExactInteger numerator = mantissa;
    ExactInteger denominator = 1;
    if (exponent >= 0) numerator <<= exponent;
    else denominator <<= -exponent;
    return rational(std::move(numerator), std::move(denominator));
}

} // namespace m2424::experimental
