#include "m2424/eval_mod.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
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
    return multiplyRelinearizeAndRescale(adapter, lhs, rhs);
}

Cipher weighted_term(SealAdapter& adapter, const Cipher& value, double coefficient) {
    return multiplyPlainAndRescale(adapter, value, adapter.encodeScalarFor(coefficient, value));
}

Cipher squash_scale(SealAdapter& adapter, Cipher value, double max_scale_log2) {
    while (std::log2(adapter.info(value).scale) > max_scale_log2 + 0.5) {
        if (adapter.info(value).chainIndex == 0) {
            throw std::runtime_error("not enough levels to squash EvalMod scale");
        }
        value = adapter.rescaleToNext(value);
    }
    return value;
}

template <class Fn>
Cipher evalmod_step(const char* name, Fn&& fn) {
    try {
        return fn();
    } catch (const std::exception& e) {
        std::ostringstream out;
        out << "EvalMod " << name << " failed: " << e.what();
        throw std::runtime_error(out.str());
    }
}

double next_rescale_drop_log2(SealAdapter& adapter, const Cipher& value) {
    const auto info = adapter.info(value);
    const auto bits = adapter.coeffModulusBits();
    if (info.coeffModulusSize == 0 || info.coeffModulusSize > bits.size()) {
        throw std::runtime_error("cannot infer next rescale modulus size");
    }
    return static_cast<double>(bits[info.coeffModulusSize - 1]);
}

Cipher weighted_term_to_scale(SealAdapter& adapter,
                              const Cipher& value,
                              double coefficient,
                              double target_scale_log2) {
    const auto info = adapter.info(value);
    const double plain_scale_log2 =
        target_scale_log2 + next_rescale_drop_log2(adapter, value) - std::log2(info.scale);
    if (!std::isfinite(plain_scale_log2) || plain_scale_log2 <= 0.0) {
        throw std::runtime_error("cannot align EvalMod term scale");
    }
    return multiplyPlainAndRescale(adapter, 
        value,
        adapter.encodeScalarAtScaleFor(coefficient, std::exp2(plain_scale_log2), value));
}

Cipher add_terms(SealAdapter& adapter, Cipher lhs, Cipher rhs) {
    try {
        lhs = adapter.alignForAddition(lhs, rhs);
    } catch (const std::exception& e) {
        const auto lhs_info = adapter.info(lhs);
        const auto rhs_info = adapter.info(rhs);
        std::ostringstream out;
        out << e.what()
            << "; lhs_chain=" << lhs_info.chainIndex
            << "; rhs_chain=" << rhs_info.chainIndex
            << "; lhs_scale_log2=" << std::log2(lhs_info.scale)
            << "; rhs_scale_log2=" << std::log2(rhs_info.scale)
            << "; lhs_coeff_modulus_log2=" << lhs_info.coeffModulusLog2
            << "; rhs_coeff_modulus_log2=" << rhs_info.coeffModulusLog2;
        throw std::runtime_error(out.str());
    }
    return adapter.add(lhs, rhs);
}

bool evalmod_diagnostics_enabled() {
    const char* value = std::getenv("M2424_DIAG_EVALMOD");
    return value != nullptr && value[0] != '\0' && value[0] != '0';
}

void print_evalmod_info(SealAdapter& adapter, const char* label, const Cipher& value) {
    if (!evalmod_diagnostics_enabled()) {
        return;
    }
    const auto info = adapter.info(value);
    std::printf("[eval_mod] %s chain=%zu scale_log2=%.6f coeffModulusLog2=%.0f\n",
                label,
                info.chainIndex,
                std::log2(info.scale),
                info.coeffModulusLog2);
}

} // namespace

double EvalModPolynomial::evaluate_plain(double u) const {
    return evaluate_plain(u, EvalModDegree::P7);
}

double EvalModPolynomial::evaluate_plain(double u, EvalModDegree degree) const {
    validate_input(u);
    if (degree == EvalModDegree::P3DoubleAngle) {
        const double half = 0.5 * u;
        const double half2 = half * half;
        const double half3 = half * half2;
        const double p3 = a1 * half + a3 * half3;
        return 2.0 * p3 - 4.0 * kPi * kPi * p3 * p3 * p3;
    }
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
    if (degree == EvalModDegree::P3DoubleAngle) {
        const Complex half = 0.5 * u;
        const Complex half2 = half * half;
        const Complex half3 = half * half2;
        const Complex p3 = a1 * half + a3 * half3;
        return 2.0 * p3 - 4.0 * kPi * kPi * p3 * p3 * p3;
    }
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
    if (degree == EvalModDegree::P3DoubleAngle) {
        const Cipher half = evalmod_step("P3DoubleAngle half-scale", [&] {
            return weighted_term(adapter, input, 0.5);
        });
        const Cipher p3 = evalmod_step("P3DoubleAngle inner P3", [&] {
            return evaluate(adapter, half, EvalModDegree::P3);
        });
        const Cipher p3_scaled = evalmod_step("P3DoubleAngle scale squash", [&] {
            return squash_scale(adapter, p3, 60.0);
        });
        const Cipher p3_squared = evalmod_step("P3DoubleAngle square", [&] {
            return multiply_same_level(adapter, p3_scaled, p3_scaled);
        });
        const Cipher p3_cubed = evalmod_step("P3DoubleAngle cube", [&] {
            return multiply_same_level(adapter, adapter.modSwitchTo(p3_scaled, p3_squared), p3_squared);
        });
        Cipher cubic = evalmod_step("P3DoubleAngle cubic weight", [&] {
            return weighted_term(adapter, p3_cubed, -4.0 * kPi * kPi);
        });
        const double target_scale_log2 = std::log2(adapter.info(cubic).scale);
        Cipher linear = evalmod_step("P3DoubleAngle linear weight", [&] {
            return weighted_term_to_scale(
                adapter,
                adapter.modSwitchTo(p3_scaled, p3_cubed),
                2.0,
                target_scale_log2);
        });
        return evalmod_step("P3DoubleAngle add", [&] {
            return add_terms(adapter, std::move(linear), std::move(cubic));
        });
    }

    const Cipher u2 = multiply_same_level(adapter, input, input);
    print_evalmod_info(adapter, "P3 input", input);
    print_evalmod_info(adapter, "P3 u2", u2);
    const Cipher u3 = multiply_same_level(adapter, adapter.modSwitchTo(input, u2), u2);
    print_evalmod_info(adapter, "P3 u3", u3);
    if (degree == EvalModDegree::P3) {
        const double target_scale_log2 = std::log2(adapter.info(input).scale);
        Cipher cubic = weighted_term_to_scale(adapter, u3, a3, target_scale_log2);
        print_evalmod_info(adapter, "P3 cubic_weighted", cubic);
        Cipher linear = weighted_term_to_scale(adapter, adapter.modSwitchTo(input, u3), a1, target_scale_log2);
        print_evalmod_info(adapter, "P3 linear_weighted", linear);
        Cipher result = add_terms(adapter, std::move(linear), std::move(cubic));
        print_evalmod_info(adapter, "P3 add_result", result);
        return result;
    }

    const Cipher u5 = multiply_same_level(adapter, u3, adapter.modSwitchTo(u2, u3));
    if (degree == EvalModDegree::P5) {
        Cipher result = weighted_term(adapter, adapter.modSwitchTo(input, u5), a1);
        result = add_terms(adapter, std::move(result), weighted_term(adapter, adapter.modSwitchTo(u3, u5), a3));
        result = add_terms(adapter, std::move(result), weighted_term(adapter, u5, a5));
        return result;
    }

    const Cipher u7 = multiply_same_level(adapter, u5, adapter.modSwitchTo(u2, u5));
    Cipher result = weighted_term(adapter, adapter.modSwitchTo(input, u7), a1);
    result = add_terms(adapter, std::move(result), weighted_term(adapter, adapter.modSwitchTo(u3, u7), a3));
    result = add_terms(adapter, std::move(result), weighted_term(adapter, adapter.modSwitchTo(u5, u7), a5));
    result = add_terms(adapter, std::move(result), weighted_term(adapter, u7, a7));
    return result;
}

const char* to_string(EvalModDegree degree) noexcept {
    switch (degree) {
    case EvalModDegree::P3:
        return "P3";
    case EvalModDegree::P3DoubleAngle:
        return "P3DoubleAngle";
    case EvalModDegree::P5:
        return "P5";
    case EvalModDegree::P7:
        return "P7";
    }
    return "unknown";
}

} // namespace m2424
