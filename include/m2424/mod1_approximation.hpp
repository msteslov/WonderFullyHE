#pragma once

#include "m2424/bootstrap_plan.hpp"
#include "m2424/eval_mod.hpp"

#include <string>
#include <vector>

namespace m2424 {

enum class PolynomialEvaluationStrategy {
    DirectOddPowers,
    PatersonStockmeyer,
    BabyStepGiantStep
};

struct Mod1Approximation {
    BootstrapMod1Type type{BootstrapMod1Type::LegacySineP3};
    std::size_t polynomial_degree{3};
    std::size_t double_angle{};
    double input_bound{EvalModPolynomial::approximation_bound};
    double target_error{};
    double evalmod_log_scale{40.0};
    std::vector<double> coefficients;
    std::size_t estimated_depth{};
    PolynomialEvaluationStrategy strategy{PolynomialEvaluationStrategy::DirectOddPowers};
    std::string construction_note;
};

struct BootstrapWideMod1Plan {
    double period_log2{};
    double input_bound{};
    double compressed_bound{0.125};
    std::size_t double_angle_steps{};
    std::size_t polynomial_degree{15};
    double target_error{1e-9};
    double evalmod_scale_log2{59.0};
    BootstrapMod1Type type{BootstrapMod1Type::CosDiscrete};
};

const char* to_string(PolynomialEvaluationStrategy strategy) noexcept;
Mod1Approximation make_mod1_approximation(const BootstrapMod1Model& model);
Mod1Approximation make_wide_mod1_approximation(const BootstrapWideMod1Plan& plan);

double evaluate_polynomial_plain(const Mod1Approximation& approximation, double input);
Complex evaluate_polynomial_plain(const Mod1Approximation& approximation, Complex input);
std::vector<double> evaluate_polynomial_plain(const Mod1Approximation& approximation,
                                              const std::vector<double>& input);
ComplexVector evaluate_polynomial_plain(const Mod1Approximation& approximation,
                                        const ComplexVector& input);
double evaluate_mod1_plain_with_double_angle(const Mod1Approximation& approximation, double input);
Complex evaluate_mod1_plain_with_double_angle(const Mod1Approximation& approximation, Complex input);
std::vector<double> evaluate_mod1_plain_with_double_angle(const Mod1Approximation& approximation,
                                                          const std::vector<double>& input);
ComplexVector evaluate_mod1_plain_with_double_angle(const Mod1Approximation& approximation,
                                                    const ComplexVector& input);
double evaluate_wide_mod1_plain(const Mod1Approximation& approximation, double input);
Complex evaluate_wide_mod1_plain(const Mod1Approximation& approximation, Complex input);

} // namespace m2424
