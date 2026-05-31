#include "bootstrap_prototype_detail.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace m2424::bootstrap_prototype_detail {

double max_complex_error(const ComplexVector& expected, const ComplexVector& actual) {
    if (expected.size() != actual.size()) {
        throw std::invalid_argument("vectors must have equal size");
    }
    const std::size_t n = expected.size();
    double result = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        result = std::max(result, std::abs(expected[i] - actual[i]));
    }
    return result;
}

ComplexVector head(ComplexVector values, std::size_t n) {
    values.resize(std::min(values.size(), n));
    return values;
}

ComplexVector scaled(const ComplexVector& values, double factor) {
    ComplexVector result;
    result.reserve(values.size());
    for (const auto& value : values) {
        result.push_back(value * factor);
    }
    return result;
}

ComplexVector evaluate_plain(const EvalModPolynomial& eval_mod,
                             const ComplexVector& input,
                             EvalModDegree degree) {
    ComplexVector result;
    result.reserve(input.size());
    for (const auto& value : input) {
        result.push_back(eval_mod.evaluate_plain(value, degree));
    }
    return result;
}

double max_abs_value(const ComplexVector& values) {
    double result = 0.0;
    for (const auto& value : values) {
        result = std::max(result, std::abs(value));
    }
    return result;
}

double normalization_factor_for(const ComplexVector& values) {
    const double max_abs = max_abs_value(values);
    if (max_abs == 0.0) {
        return 1.0;
    }
    constexpr double target = EvalModPolynomial::approximation_bound * 0.5;
    return max_abs > target ? target / max_abs : 1.0;
}

BootstrapPrototypeStage make_stage(const std::string& name,
                                   const CipherInfo& before,
                                   const CipherInfo& after,
                                   double max_error,
                                   double tolerance,
                                   double duration_ms,
                                   bool checked) {
    return BootstrapPrototypeStage{
        name,
        checked ? (max_error <= tolerance ? "PASS" : "FAIL") : "RUN",
        before.chain_index,
        after.chain_index,
        before.coeff_modulus_size,
        after.coeff_modulus_size,
        before.coeff_modulus_log2,
        after.coeff_modulus_log2,
        before.scale,
        after.scale,
        max_error,
        duration_ms
    };
}

BootstrapPrototypeStage make_harness_stage(const CipherInfo& after) {
    return BootstrapPrototypeStage{
        "mod_raise_harness",
        "PASS",
        after.chain_index,
        after.chain_index,
        after.coeff_modulus_size,
        after.coeff_modulus_size,
        after.coeff_modulus_log2,
        after.coeff_modulus_log2,
        after.scale,
        after.scale,
        0.0,
        0.0
    };
}

Cipher apply_output_scale_repair(SealAdapter& adapter, const Cipher& input, double target_scale_log2) {
    auto current = input;
    auto info = adapter.info(current);
    double current_scale_log2 = std::log2(info.scale);
    if (current_scale_log2 <= target_scale_log2 + 0.5) {
        return current;
    }
    if (info.chain_index == 0) {
        throw std::runtime_error("cannot repair bootstrap output scale without remaining levels");
    }
    const auto bits = adapter.coeff_modulus_bits();
    while (true) {
        if (info.coeff_modulus_size == 0 || info.coeff_modulus_size > bits.size()) {
            throw std::runtime_error("cannot infer bootstrap output rescale modulus size");
        }
        const double next_drop_log2 = static_cast<double>(bits[info.coeff_modulus_size - 1]);
        if (current_scale_log2 <= target_scale_log2 + next_drop_log2 + 0.5) {
            break;
        }
        if (info.chain_index == 0) {
            throw std::runtime_error("cannot repair bootstrap output scale without remaining levels");
        }
        current = adapter.rescale_to_next(current);
        info = adapter.info(current);
        current_scale_log2 = std::log2(info.scale);
        if (current_scale_log2 <= target_scale_log2 + 0.5) {
            return current;
        }
    }
    if (info.coeff_modulus_size == 0 || info.coeff_modulus_size > bits.size()) {
        throw std::runtime_error("cannot infer bootstrap output rescale modulus size");
    }
    const double next_drop_log2 = static_cast<double>(bits[info.coeff_modulus_size - 1]);
    const double plain_scale_log2 = target_scale_log2 + next_drop_log2 - current_scale_log2;
    if (!std::isfinite(plain_scale_log2) || plain_scale_log2 <= 0.0) {
        std::ostringstream out;
        out << "cannot repair bootstrap output scale"
            << "; current_scale_log2=" << current_scale_log2
            << "; target_scale_log2=" << target_scale_log2
            << "; next_drop_log2=" << next_drop_log2
            << "; required_plain_scale_log2=" << plain_scale_log2
            << "; chain_index=" << info.chain_index;
        throw std::runtime_error(out.str());
    }
    return adapter.mul_plain_rescale(
        current,
        adapter.encode_scalar_at_scale_like(1.0, std::exp2(plain_scale_log2), current));
}

void mark_stage_structural(BootstrapPrototypeStage& stage) {
    stage.status = "STRUCTURAL";
    stage.max_abs_error = 0.0;
}

void mark_stage_diagnostic(BootstrapPrototypeStage& stage) {
    stage.status = "DIAG";
}

double finite_exp2_or_zero(double exponent) {
    if (exponent < -1074.0) {
        return 0.0;
    }
    if (exponent > 1023.0) {
        return std::numeric_limits<double>::infinity();
    }
    return std::exp2(exponent);
}

} // namespace m2424::bootstrap_prototype_detail
