#pragma once

#include <gmpxx.h>

#include <cstdint>
#include <vector>

namespace m2424 {

using ExactInteger = mpz_class;
using OracleFloat = mpf_class;

struct ExactCoefficientOracleResult {
    ExactInteger reconstructed;
    ExactInteger centeredSourceCoefficient;
    OracleFloat expectedValue;
};

/// Exact CRT in [0,Q), followed by exact centered reduction modulo qSource.
ExactCoefficientOracleResult exactCoefficientOracle(const std::vector<std::uint64_t>& residues,
                                                     const std::vector<std::uint64_t>& moduli,
                                                     const ExactInteger& qSource,
                                                     const OracleFloat& outputScale);

} // namespace m2424
