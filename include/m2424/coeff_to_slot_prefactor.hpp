#pragma once

#include "m2424/bootstrap_contract.hpp"

namespace m2424 {

/// Positive exact dyadic / factored-integer scalar. Default is exactly one.
/// No qSource product or scalar division is performed in floating point.
class CoeffToSlotPrefactor {
public:
    CoeffToSlotPrefactor() = default;
    static CoeffToSlotPrefactor sourceNormalization(const BootstrapInputContext&);
    bool isIdentity() const noexcept;
    bool operator==(const CoeffToSlotPrefactor&) const noexcept;
    /// Round the exact product scalar * binary64 operand once, ties to even.
    /// Used only at the plaintext-constant boundary (or plaintext diagnostics).
    double multiplyRounded(double operand) const;
    std::uint64_t numeratorScaleBits() const noexcept { return scaleBits_; }
    const std::vector<std::uint64_t>& denominatorFactors() const noexcept { return primes_; }
private:
    std::uint64_t scaleBits_{0x3ff0000000000000ULL};
    std::vector<std::uint64_t> primes_;
    std::vector<std::uint64_t> denominator_{1};
};

} // namespace m2424
