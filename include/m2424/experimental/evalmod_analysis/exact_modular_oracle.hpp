#pragma once

#include <gmpxx.h>

#include <cstdint>
#include <cstddef>
#include <string>
#include <vector>

namespace m2424::experimental {

using ExactInteger = mpz_class;
struct ExactCoefficientOracleResult {
    ExactInteger reconstructed;
    ExactInteger crtModulus;
    ExactInteger centeredSourceCoefficient;
    std::string expectedValueDecimal;
    std::size_t precisionBits{};
};

/// Exact CRT in [0,Q), followed by exact centered reduction modulo qSource.
ExactCoefficientOracleResult exactCoefficientOracle(const std::vector<std::uint64_t>& residues,
                                                     const std::vector<std::uint64_t>& moduli,
                                                     const ExactInteger& qSource,
                                                     const std::string& outputScaleDecimal,
                                                     std::size_t precisionBits);

} // namespace m2424::experimental
