#include "m2424/exact_modular_oracle.hpp"

#include <stdexcept>

namespace m2424 {
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
                                                     const OracleFloat& outputScale) {
    if (residues.empty() || residues.size() != moduli.size() || qSource <= 1 || outputScale <= 0) {
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
    ExactInteger centered = reconstructed % qSource;
    if (centered < 0) centered += qSource;
    if (2 * centered >= qSource) centered -= qSource;
    return {reconstructed, centered, OracleFloat(centered) / outputScale};
}

} // namespace m2424
