#pragma once

#include "m2424/bootstrap_plan.hpp"
#include "m2424/eval_mod.hpp"

#include <string>
#include <vector>

namespace m2424 {

struct Mod1Approximation {
    BootstrapMod1Type type{BootstrapMod1Type::LegacySineP3};
    std::size_t polynomial_degree{3};
    std::size_t double_angle{};
    double input_bound{EvalModPolynomial::approximation_bound};
    double target_error{};
    double evalmod_log_scale{40.0};
    std::vector<double> coefficients;
    std::size_t estimated_depth{};
    std::string construction_note;
};

Mod1Approximation make_mod1_approximation(const BootstrapMod1Model& model);

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

} // namespace m2424
