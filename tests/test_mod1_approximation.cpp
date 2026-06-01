#include "m2424/mod1_approximation.hpp"

#include <cmath>
#include <cstdio>
#include <exception>
#include <stdexcept>
#include <vector>

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

void require(bool condition, const char* message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

double sine_reference(double input) {
    return std::sin(2.0 * kPi * input) / (2.0 * kPi);
}

m2424::Complex sine_reference(m2424::Complex input) {
    return std::sin(2.0 * kPi * input) / (2.0 * kPi);
}

} // namespace

int main() {
    bool ok = true;
    try {
        const auto approximation = m2424::make_mod1_approximation(
            {m2424::BootstrapMod1Type::CosDiscrete, 31, 3, 8, 60.0});
        require(approximation.type == m2424::BootstrapMod1Type::CosDiscrete,
                "Mod1 approximation type mismatch");
        require(approximation.polynomial_degree == 31,
                "Mod1 approximation degree mismatch");
        require(approximation.double_angle == 3,
                "Mod1 approximation double-angle mismatch");
        require(approximation.input_bound >= 0.1,
                "Mod1 approximation input bound should cover plaintext validation matrix");
        require(approximation.coefficients.size() == 32,
                "Mod1 approximation coefficients should be inspectable");
        require(approximation.strategy == m2424::PolynomialEvaluationStrategy::PatersonStockmeyer,
                "degree-31 approximation must not use direct odd powers");
        require(!approximation.construction_note.empty(),
                "Mod1 approximation construction note should be explicit");

        const std::vector<double> inputs{0.0, 1e-5, -1e-5, 1e-3, -1e-3, 1e-2, -1e-2, 1e-1, -1e-1};
        double max_error = 0.0;
        for (double input : inputs) {
            const double actual = m2424::evaluate_polynomial_plain(approximation, input);
            const double expected = sine_reference(input);
            max_error = std::max(max_error, std::abs(actual - expected));
        }
        require(max_error < 1e-14, "degree-31 coefficient polynomial should match sine reference on test grid");

        const m2424::ComplexVector complex_inputs{{1e-3, 1e-4}, {-1e-2, 5e-4}, {1e-1, -1e-3}};
        const auto complex_outputs = m2424::evaluate_polynomial_plain(approximation, complex_inputs);
        double max_complex_error = 0.0;
        for (std::size_t i = 0; i < complex_inputs.size(); ++i) {
            max_complex_error = std::max(max_complex_error,
                                         std::abs(complex_outputs[i] - sine_reference(complex_inputs[i])));
        }
        require(max_complex_error < 1e-13,
                "degree-31 complex coefficient polynomial should match sine reference");

        const auto double_angle = m2424::evaluate_mod1_plain_with_double_angle(approximation, inputs);
        require(double_angle.size() == inputs.size(), "double-angle vector output size mismatch");

        bool invalid_degree_threw = false;
        try {
            (void)m2424::make_mod1_approximation({m2424::BootstrapMod1Type::CosDiscrete, 30, 0, 8, 60.0});
        } catch (const std::invalid_argument&) {
            invalid_degree_threw = true;
        }
        require(invalid_degree_threw, "invalid even CosDiscrete degree should throw");

        bool invalid_scale_threw = false;
        try {
            (void)m2424::make_mod1_approximation({m2424::BootstrapMod1Type::CosDiscrete, 31, 0, 8, -1.0});
        } catch (const std::invalid_argument&) {
            invalid_scale_threw = true;
        }
        require(invalid_scale_threw, "invalid Mod1 scale should throw");

        std::printf("[test_mod1_approximation] PASS max_error=%.3e max_complex_error=%.3e coeffs=%zu strategy=%s note=%s\n",
                    max_error,
                    max_complex_error,
                    approximation.coefficients.size(),
                    m2424::to_string(approximation.strategy),
                    approximation.construction_note.c_str());
    } catch (const std::exception& error) {
        ok = false;
        std::printf("[test_mod1_approximation] FAIL: %s\n", error.what());
    }

    return ok ? 0 : 1;
}
