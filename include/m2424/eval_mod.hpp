#pragma once

#include "m2424/seal_adapter.hpp"
#include "m2424/diagonal_transform.hpp"

#include <vector>

namespace m2424 {

class EvalModPolynomial {
public:
    static constexpr double approximation_bound = 0.0009765625;
    static constexpr double a1 = 1.0;
    static constexpr double a3 = -6.579736267393;
    static constexpr double a5 = 12.98787880453;
    static constexpr double a7 = -12.20811674381;

    double evaluate_plain(double u) const;
    std::vector<double> evaluate_plain(const std::vector<double>& input) const;
    Complex evaluate_plain(Complex u) const;
    ComplexVector evaluate_plain(const ComplexVector& input) const;
    double sine_reference(double u) const;
    std::vector<double> sine_reference(const std::vector<double>& input) const;
    Cipher evaluate(SealAdapter& adapter, const Cipher& input) const;
};

} // namespace m2424
