#pragma once
#include <string>
#include <vector>
// Shared storage only; exact parsing remains in the optional analysis library.
namespace m2424::experimental {
enum class PolynomialBasis { Monomial, Chebyshev, Composite };

struct EvalModPolynomial {
    PolynomialBasis basis{PolynomialBasis::Monomial};
    std::vector<std::string> decimalCoefficients;
};

}
