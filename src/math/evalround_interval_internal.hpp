#pragma once

#include "m2424/decimal_polynomial.hpp"

#include <cstdint>
#include <gmpxx.h>
#include <vector>

namespace m2424::experimental::detail {

/// Exact Taylor shift: returns coefficients d_k such that
/// p(center+y) = sum_k d_k y^k. Decimal inputs are parsed as exact rationals.
std::vector<mpq_class> shiftPolynomialToIntegerCenterExact(
    const EvalModPolynomial& polynomial, std::int64_t center);

} // namespace m2424::experimental::detail
