#include "m2424/mod1_circuit.hpp"

#include <cmath>
#include <stdexcept>

namespace m2424 {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

void validate_model(const BootstrapMod1Model& model) {
    if (model.degree == 0) {
        throw std::invalid_argument("Mod1 circuit degree must be positive");
    }
    if (model.type == BootstrapMod1Type::LegacySineP3 && model.degree != 3) {
        throw std::invalid_argument("LegacySineP3 Mod1 circuit supports only degree 3");
    }
    if (model.type == BootstrapMod1Type::CosDiscrete && model.degree < 3) {
        throw std::invalid_argument("CosDiscrete Mod1 circuit degree must be at least 3");
    }
    if (model.type == BootstrapMod1Type::CosDiscrete && (model.degree % 2) == 0) {
        throw std::invalid_argument("CosDiscrete Mod1 circuit degree must be odd");
    }
    if (!std::isfinite(model.evalmod_log_scale) || model.evalmod_log_scale <= 0.0) {
        throw std::invalid_argument("Mod1 circuit evalmod_log_scale must be positive and finite");
    }
}

EvalModDegree legacy_degree(const BootstrapMod1Model& model) {
    return model.double_angle == 0 ? EvalModDegree::P3 : EvalModDegree::P3DoubleAngle;
}

double sine_mod1_reference(double input) {
    if (!std::isfinite(input)) {
        throw std::invalid_argument("Mod1 input must be finite");
    }
    return std::sin(2.0 * kPi * input) / (2.0 * kPi);
}

Complex sine_mod1_reference(Complex input) {
    if (!std::isfinite(input.real()) || !std::isfinite(input.imag())) {
        throw std::invalid_argument("Mod1 complex input must be finite");
    }
    return std::sin(2.0 * kPi * input) / (2.0 * kPi);
}

} // namespace

Mod1Circuit::Mod1Circuit(BootstrapMod1Model model) : model_(model) {
    validate_model(model_);
}

const BootstrapMod1Model& Mod1Circuit::model() const noexcept {
    return model_;
}

std::size_t Mod1Circuit::estimated_levels() const {
    return estimated_bootstrap_mod1_levels(model_);
}

bool Mod1Circuit::encrypted_evaluation_available() const noexcept {
    return model_.type == BootstrapMod1Type::LegacySineP3
        || (model_.type == BootstrapMod1Type::CosDiscrete && model_.degree == 3 && model_.double_angle <= 1);
}

double Mod1Circuit::evaluate_plain(double input) const {
    if (model_.type == BootstrapMod1Type::LegacySineP3) {
        return EvalModPolynomial{}.evaluate_plain(input, legacy_degree(model_));
    }
    return sine_mod1_reference(input);
}

Complex Mod1Circuit::evaluate_plain(Complex input) const {
    if (model_.type == BootstrapMod1Type::LegacySineP3) {
        return EvalModPolynomial{}.evaluate_plain(input, legacy_degree(model_));
    }
    return sine_mod1_reference(input);
}

std::vector<double> Mod1Circuit::evaluate_plain(const std::vector<double>& input) const {
    if (input.empty()) {
        throw std::invalid_argument("Mod1 input vector must not be empty");
    }
    std::vector<double> result;
    result.reserve(input.size());
    for (double value : input) {
        result.push_back(evaluate_plain(value));
    }
    return result;
}

ComplexVector Mod1Circuit::evaluate_plain(const ComplexVector& input) const {
    if (input.empty()) {
        throw std::invalid_argument("Mod1 complex input vector must not be empty");
    }
    ComplexVector result;
    result.reserve(input.size());
    for (Complex value : input) {
        result.push_back(evaluate_plain(value));
    }
    return result;
}

Cipher Mod1Circuit::evaluate(SealAdapter& adapter, const Cipher& input) const {
    if (!encrypted_evaluation_available()) {
        throw std::runtime_error("encrypted CosDiscrete Mod1 circuit is not implemented for this model");
    }
    if (model_.type == BootstrapMod1Type::LegacySineP3) {
        return EvalModPolynomial{}.evaluate(adapter, input, legacy_degree(model_));
    }
    const EvalModDegree degree = model_.double_angle == 0 ? EvalModDegree::P3 : EvalModDegree::P3DoubleAngle;
    return EvalModPolynomial{}.evaluate(adapter, input, degree);
}

} // namespace m2424
