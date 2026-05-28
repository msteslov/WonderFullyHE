#include "m2424/eval_mod.hpp"

#include <cmath>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace m2424 {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

void validate_input(double u) {
    if (!std::isfinite(u)) {
        throw std::invalid_argument("EvalMod input must be finite");
    }
    if (std::fabs(u) > EvalModPolynomial::approximation_bound) {
        throw std::invalid_argument("EvalMod input is outside approximation interval");
    }
}

void validate_input(Complex u) {
    if (!std::isfinite(u.real()) || !std::isfinite(u.imag())) {
        throw std::invalid_argument("EvalMod input must be finite");
    }
    if (std::abs(u) > EvalModPolynomial::approximation_bound) {
        throw std::invalid_argument("EvalMod input is outside approximation interval");
    }
}

Cipher multiply_same_level(SealAdapter& adapter, const Cipher& lhs, const Cipher& rhs) {
    return adapter.mul_relin_rescale(lhs, rhs);
}

Cipher weighted_term(SealAdapter& adapter, const Cipher& value, double coefficient) {
    return adapter.mul_plain_rescale(value, adapter.encode_scalar_like(coefficient, value));
}

Cipher add_terms(SealAdapter& adapter, Cipher lhs, Cipher rhs) {
    try {
        lhs = adapter.match_level_and_scale(lhs, rhs);
    } catch (const std::exception& e) {
        const auto lhs_info = adapter.info(lhs);
        const auto rhs_info = adapter.info(rhs);
        std::ostringstream out;
        out << e.what()
            << "; lhs_chain=" << lhs_info.chain_index
            << "; rhs_chain=" << rhs_info.chain_index
            << "; lhs_scale_log2=" << std::log2(lhs_info.scale)
            << "; rhs_scale_log2=" << std::log2(rhs_info.scale)
            << "; lhs_coeff_modulus_log2=" << lhs_info.coeff_modulus_log2
            << "; rhs_coeff_modulus_log2=" << rhs_info.coeff_modulus_log2;
        throw std::runtime_error(out.str());
    }
    return adapter.add(lhs, rhs);
}

} // namespace

double EvalModPolynomial::evaluate_plain(double u) const {
    return evaluate_plain(u, EvalModDegree::P7);
}

double EvalModPolynomial::evaluate_plain(double u, EvalModDegree degree) const {
    validate_input(u);
    const double u2 = u * u;
    const double u3 = u * u2;
    double result = a1 * u + a3 * u3;
    if (degree == EvalModDegree::P3) {
        return result;
    }
    const double u5 = u3 * u2;
    result += a5 * u5;
    if (degree == EvalModDegree::P5) {
        return result;
    }
    const double u7 = u5 * u2;
    return result + a7 * u7;
}

std::vector<double> EvalModPolynomial::evaluate_plain(const std::vector<double>& input) const {
    if (input.empty()) {
        throw std::invalid_argument("EvalMod input vector must not be empty");
    }
    std::vector<double> result;
    result.reserve(input.size());
    for (double value : input) {
        result.push_back(evaluate_plain(value));
    }
    return result;
}

Complex EvalModPolynomial::evaluate_plain(Complex u) const {
    return evaluate_plain(u, EvalModDegree::P7);
}

Complex EvalModPolynomial::evaluate_plain(Complex u, EvalModDegree degree) const {
    validate_input(u);
    const Complex u2 = u * u;
    const Complex u3 = u * u2;
    Complex result = a1 * u + a3 * u3;
    if (degree == EvalModDegree::P3) {
        return result;
    }
    const Complex u5 = u3 * u2;
    result += a5 * u5;
    if (degree == EvalModDegree::P5) {
        return result;
    }
    const Complex u7 = u5 * u2;
    return result + a7 * u7;
}

ComplexVector EvalModPolynomial::evaluate_plain(const ComplexVector& input) const {
    if (input.empty()) {
        throw std::invalid_argument("EvalMod input vector must not be empty");
    }
    ComplexVector result;
    result.reserve(input.size());
    for (Complex value : input) {
        result.push_back(evaluate_plain(value));
    }
    return result;
}

double EvalModPolynomial::sine_reference(double u) const {
    validate_input(u);
    return std::sin(2.0 * kPi * u) / (2.0 * kPi);
}

std::vector<double> EvalModPolynomial::sine_reference(const std::vector<double>& input) const {
    if (input.empty()) {
        throw std::invalid_argument("EvalMod input vector must not be empty");
    }
    std::vector<double> result;
    result.reserve(input.size());
    for (double value : input) {
        result.push_back(sine_reference(value));
    }
    return result;
}

Cipher EvalModPolynomial::evaluate(SealAdapter& adapter, const Cipher& input) const {
    return evaluate(adapter, input, EvalModDegree::P7);
}

Cipher EvalModPolynomial::evaluate(SealAdapter& adapter, const Cipher& input, EvalModDegree degree) const {
    const Cipher u2 = multiply_same_level(adapter, input, input);
    const Cipher u3 = multiply_same_level(adapter, adapter.mod_switch_to(input, u2), u2);
    if (degree == EvalModDegree::P3) {
        Cipher result = weighted_term(adapter, adapter.mod_switch_to(input, u3), a1);
        result = add_terms(adapter, std::move(result), weighted_term(adapter, u3, a3));
        return result;
    }

    const Cipher u5 = multiply_same_level(adapter, u3, adapter.mod_switch_to(u2, u3));
    if (degree == EvalModDegree::P5) {
        Cipher result = weighted_term(adapter, adapter.mod_switch_to(input, u5), a1);
        result = add_terms(adapter, std::move(result), weighted_term(adapter, adapter.mod_switch_to(u3, u5), a3));
        result = add_terms(adapter, std::move(result), weighted_term(adapter, u5, a5));
        return result;
    }

    const Cipher u7 = multiply_same_level(adapter, u5, adapter.mod_switch_to(u2, u5));
    Cipher result = weighted_term(adapter, adapter.mod_switch_to(input, u7), a1);
    result = add_terms(adapter, std::move(result), weighted_term(adapter, adapter.mod_switch_to(u3, u7), a3));
    result = add_terms(adapter, std::move(result), weighted_term(adapter, adapter.mod_switch_to(u5, u7), a5));
    result = add_terms(adapter, std::move(result), weighted_term(adapter, u7, a7));
    return result;
}

const char* to_string(EvalModDegree degree) noexcept {
    switch (degree) {
    case EvalModDegree::P3:
        return "P3";
    case EvalModDegree::P5:
        return "P5";
    case EvalModDegree::P7:
        return "P7";
    }
    return "unknown";
}

} // namespace m2424
