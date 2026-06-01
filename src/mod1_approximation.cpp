#include "m2424/mod1_approximation.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace m2424 {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

void validate_common(const BootstrapMod1Model& model) {
    if (model.degree == 0) {
        throw std::invalid_argument("Mod1 approximation degree must be positive");
    }
    if (!std::isfinite(model.evalmod_log_scale) || model.evalmod_log_scale <= 0.0) {
        throw std::invalid_argument("Mod1 approximation evalmod_log_scale must be positive and finite");
    }
}

void validate_input_bound(double input_bound) {
    if (!std::isfinite(input_bound) || input_bound <= 0.0) {
        throw std::invalid_argument("Mod1 approximation input_bound must be positive and finite");
    }
}

void validate_input(double input, double input_bound) {
    if (!std::isfinite(input)) {
        throw std::invalid_argument("Mod1 approximation input must be finite");
    }
    if (std::abs(input) > input_bound) {
        throw std::invalid_argument("Mod1 approximation input is outside approximation interval");
    }
}

void validate_input(Complex input, double input_bound) {
    if (!std::isfinite(input.real()) || !std::isfinite(input.imag())) {
        throw std::invalid_argument("Mod1 approximation complex input must be finite");
    }
    if (std::abs(input) > input_bound) {
        throw std::invalid_argument("Mod1 approximation complex input is outside approximation interval");
    }
}

std::vector<double> sine_coefficients_for_degree(std::size_t degree) {
    switch (degree) {
    case 3:
        return {0.0, EvalModPolynomial::a1, 0.0, EvalModPolynomial::a3};
    case 5:
        return {0.0, EvalModPolynomial::a1, 0.0, EvalModPolynomial::a3, 0.0, EvalModPolynomial::a5};
    case 7:
        return {0.0,
                EvalModPolynomial::a1,
                0.0,
                EvalModPolynomial::a3,
                0.0,
                EvalModPolynomial::a5,
                0.0,
                EvalModPolynomial::a7};
    default:
        if (degree < 3 || (degree % 2) == 0) {
            throw std::invalid_argument("sine polynomial degree must be odd and at least 3");
        }
        std::vector<double> coefficients(degree + 1, 0.0);
        for (std::size_t power = 1; power <= degree; power += 2) {
            const std::size_t k = (power - 1) / 2;
            const double sign = (k % 2 == 0) ? 1.0 : -1.0;
            double factorial = 1.0;
            for (std::size_t i = 2; i <= power; ++i) {
                factorial *= static_cast<double>(i);
            }
            coefficients[power] = sign * std::pow(2.0 * kPi, static_cast<double>(power - 1)) / factorial;
        }
        return coefficients;
    }
}

std::size_t direct_odd_power_depth(std::size_t degree) {
    switch (degree) {
    case 3:
        return 3;
    case 5:
        return 4;
    case 7:
        return 5;
    default:
        return static_cast<std::size_t>(std::ceil(std::log2(static_cast<double>(degree)))) + 2;
    }
}

template <class T>
T evaluate_coefficients(const Mod1Approximation& approximation, T input) {
    validate_input(input, approximation.input_bound);
    T result{};
    T power{1.0};
    for (double coefficient : approximation.coefficients) {
        result += coefficient * power;
        power *= input;
    }
    return result;
}

template <>
double evaluate_coefficients<double>(const Mod1Approximation& approximation, double input) {
    validate_input(input, approximation.input_bound);
    double result = 0.0;
    double power = 1.0;
    for (double coefficient : approximation.coefficients) {
        result += coefficient * power;
        power *= input;
    }
    return result;
}

template <class T>
T apply_double_angle(const Mod1Approximation& approximation, T input) {
    const double divisor = std::exp2(static_cast<double>(approximation.double_angle));
    T value = evaluate_coefficients(approximation, input / divisor);
    for (std::size_t i = 0; i < approximation.double_angle; ++i) {
        value = 2.0 * value - 4.0 * kPi * kPi * value * value * value;
    }
    return value;
}

} // namespace

Mod1Approximation make_mod1_approximation(const BootstrapMod1Model& model) {
    validate_common(model);

    Mod1Approximation approximation;
    approximation.type = model.type;
    approximation.polynomial_degree = model.degree;
    approximation.double_angle = model.double_angle;
    approximation.input_bound = EvalModPolynomial::approximation_bound;
    approximation.evalmod_log_scale = model.evalmod_log_scale;

    if (model.type == BootstrapMod1Type::LegacySineP3) {
        if (model.degree != 3) {
            throw std::invalid_argument("LegacySineP3 approximation supports only degree 3");
        }
        approximation.coefficients = sine_coefficients_for_degree(3);
        approximation.estimated_depth = model.double_angle == 0 ? 3 : 6;
        approximation.strategy = PolynomialEvaluationStrategy::DirectOddPowers;
        approximation.construction_note = "legacy sine P3 approximation; not CosDiscrete";
        return approximation;
    }

    if (model.type == BootstrapMod1Type::CosDiscrete) {
        if (model.degree < 3 || (model.degree % 2) == 0) {
            throw std::invalid_argument("CosDiscrete placeholder degree must be odd and at least 3");
        }
        approximation.input_bound = 0.125;
        approximation.coefficients = sine_coefficients_for_degree(model.degree);
        approximation.estimated_depth = direct_odd_power_depth(model.degree) + model.double_angle;
        approximation.strategy = model.degree <= 7
            ? PolynomialEvaluationStrategy::DirectOddPowers
            : PolynomialEvaluationStrategy::PatersonStockmeyer;
        approximation.construction_note =
            "small-degree sine-polynomial placeholder for CosDiscrete; real CosDiscrete is not implemented";
        return approximation;
    }

    throw std::invalid_argument("unknown Mod1 approximation type");
}

const char* to_string(PolynomialEvaluationStrategy strategy) noexcept {
    switch (strategy) {
    case PolynomialEvaluationStrategy::DirectOddPowers:
        return "DirectOddPowers";
    case PolynomialEvaluationStrategy::PatersonStockmeyer:
        return "PatersonStockmeyer";
    case PolynomialEvaluationStrategy::BabyStepGiantStep:
        return "BabyStepGiantStep";
    }
    return "unknown";
}

double evaluate_polynomial_plain(const Mod1Approximation& approximation, double input) {
    validate_input_bound(approximation.input_bound);
    return evaluate_coefficients(approximation, input);
}

Complex evaluate_polynomial_plain(const Mod1Approximation& approximation, Complex input) {
    validate_input_bound(approximation.input_bound);
    return evaluate_coefficients(approximation, input);
}

std::vector<double> evaluate_polynomial_plain(const Mod1Approximation& approximation,
                                              const std::vector<double>& input) {
    if (input.empty()) {
        throw std::invalid_argument("Mod1 approximation input vector must not be empty");
    }
    std::vector<double> result;
    result.reserve(input.size());
    for (double value : input) {
        result.push_back(evaluate_polynomial_plain(approximation, value));
    }
    return result;
}

ComplexVector evaluate_polynomial_plain(const Mod1Approximation& approximation,
                                        const ComplexVector& input) {
    if (input.empty()) {
        throw std::invalid_argument("Mod1 approximation complex input vector must not be empty");
    }
    ComplexVector result;
    result.reserve(input.size());
    for (Complex value : input) {
        result.push_back(evaluate_polynomial_plain(approximation, value));
    }
    return result;
}

double evaluate_mod1_plain_with_double_angle(const Mod1Approximation& approximation, double input) {
    validate_input_bound(approximation.input_bound);
    return apply_double_angle(approximation, input);
}

Complex evaluate_mod1_plain_with_double_angle(const Mod1Approximation& approximation, Complex input) {
    validate_input_bound(approximation.input_bound);
    return apply_double_angle(approximation, input);
}

std::vector<double> evaluate_mod1_plain_with_double_angle(const Mod1Approximation& approximation,
                                                          const std::vector<double>& input) {
    if (input.empty()) {
        throw std::invalid_argument("Mod1 approximation input vector must not be empty");
    }
    std::vector<double> result;
    result.reserve(input.size());
    for (double value : input) {
        result.push_back(evaluate_mod1_plain_with_double_angle(approximation, value));
    }
    return result;
}

ComplexVector evaluate_mod1_plain_with_double_angle(const Mod1Approximation& approximation,
                                                    const ComplexVector& input) {
    if (input.empty()) {
        throw std::invalid_argument("Mod1 approximation complex input vector must not be empty");
    }
    ComplexVector result;
    result.reserve(input.size());
    for (Complex value : input) {
        result.push_back(evaluate_mod1_plain_with_double_angle(approximation, value));
    }
    return result;
}

} // namespace m2424
