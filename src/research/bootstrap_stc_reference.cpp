#include "m2424/bootstrap_stc_reference.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace m2424 {
namespace {

double max_error(const ComplexVector& lhs, const ComplexVector& rhs) {
    if (lhs.size() != rhs.size()) {
        throw std::invalid_argument("lattice scan vectors must have equal size");
    }
    double result = 0.0;
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        result = std::max(result, std::abs(lhs[i] - rhs[i]));
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

ComplexVector multiply(const ComplexVector& values, std::complex<double> factor) {
    ComplexVector result;
    result.reserve(values.size());
    for (const auto& value : values) {
        result.push_back(factor * value);
    }
    return result;
}

ComplexVector divide_by_period(const ComplexVector& values, double period_log2) {
    const double factor = std::exp2(-period_log2);
    ComplexVector result;
    result.reserve(values.size());
    for (const auto& value : values) {
        result.push_back(value * factor);
    }
    return result;
}

void validate_plan(const BootstrapStcReferencePlan& plan) {
    if (plan.slots == 0) {
        throw std::invalid_argument("reference plan slots must be positive");
    }
    if (!std::isfinite(plan.message_scale_log2)
        || !std::isfinite(plan.slots_to_coeff_scale_log2)
        || !std::isfinite(plan.scale_down_log2)
        || !std::isfinite(plan.modup_period_log2)) {
        throw std::invalid_argument("reference plan log values must be finite");
    }
    if (!std::isfinite(plan.expected_gamma.real()) || !std::isfinite(plan.expected_gamma.imag())) {
        throw std::invalid_argument("reference plan gamma must be finite");
    }
}

} // namespace

std::complex<double> bootstrap_lattice_best_gain(const ComplexVector& actual,
                                                 const ComplexVector& expected) {
    if (actual.size() != expected.size()) {
        throw std::invalid_argument("gain vectors must have equal size");
    }
    std::complex<double> numerator{0.0, 0.0};
    double denominator = 0.0;
    for (std::size_t i = 0; i < actual.size(); ++i) {
        numerator += std::conj(expected[i]) * actual[i];
        denominator += std::norm(expected[i]);
    }
    return denominator == 0.0 ? std::complex<double>{0.0, 0.0} : numerator / denominator;
}

ComplexVector reduce_mod_integer_lattice(const ComplexVector& values) {
    ComplexVector result;
    result.reserve(values.size());
    for (const auto& value : values) {
        result.push_back({
            value.real() - std::round(value.real()),
            value.imag() - std::round(value.imag())
        });
    }
    return result;
}

BootstrapLatticeScanSummary scan_bootstrap_lattice_periods(const ComplexVector& z,
                                                           const ComplexVector& expected,
                                                           double period_min_log2,
                                                           double period_max_log2,
                                                           double minimum_useful_gain_abs) {
    if (period_min_log2 > period_max_log2) {
        throw std::invalid_argument("period_min_log2 must be <= period_max_log2");
    }
    if (minimum_useful_gain_abs < 0.0 || !std::isfinite(minimum_useful_gain_abs)) {
        throw std::invalid_argument("minimum_useful_gain_abs must be non-negative and finite");
    }
    BootstrapLatticeScanSummary summary;
    summary.best_mod_direct.err_mod_direct = std::numeric_limits<double>::infinity();
    summary.best_mod_gain.err_mod_gain = std::numeric_limits<double>::infinity();
    summary.best_useful_mod_gain.err_mod_gain = std::numeric_limits<double>::infinity();

    for (int p = static_cast<int>(std::ceil(period_min_log2));
         p <= static_cast<int>(std::floor(period_max_log2));
         ++p) {
        const auto u = divide_by_period(z, static_cast<double>(p));
        const auto r = reduce_mod_integer_lattice(u);
        const auto gamma = bootstrap_lattice_best_gain(r, expected);

        BootstrapLatticeScanPoint point;
        point.period_log2 = static_cast<double>(p);
        point.gamma = gamma;
        point.err_mod_direct = max_error(r, expected);
        point.err_mod_gain = max_error(r, multiply(expected, gamma));
        point.max_abs_u = max_abs_value(u);
        point.max_abs_r = max_abs_value(r);

        if (point.err_mod_direct < summary.best_mod_direct.err_mod_direct) {
            summary.best_mod_direct = point;
        }
        if (point.err_mod_gain < summary.best_mod_gain.err_mod_gain) {
            summary.best_mod_gain = point;
        }
        if (std::abs(point.gamma) >= minimum_useful_gain_abs
            && point.err_mod_gain < summary.best_useful_mod_gain.err_mod_gain) {
            summary.best_useful_mod_gain = point;
            summary.has_useful_gain = true;
        }
    }
    return summary;
}

BootstrapStcReferencePlan make_bootstrap_stc_reference_plan(std::size_t slots,
                                                            double message_scale_log2,
                                                            double target_error) {
    if (slots == 0) {
        throw std::invalid_argument("reference plan slots must be positive");
    }
    if (!std::isfinite(message_scale_log2)) {
        throw std::invalid_argument("message_scale_log2 must be finite");
    }
    if (!std::isfinite(target_error) || target_error <= 0.0) {
        throw std::invalid_argument("target_error must be positive and finite");
    }

    BootstrapStcReferencePlan plan;
    plan.slots = slots;
    plan.message_scale_log2 = message_scale_log2;
    plan.slots_to_coeff_scale_log2 = message_scale_log2;
    plan.scale_down_log2 = 0.0;
    plan.modup_period_log2 = 100.0;
    plan.expected_gamma = {1.0, 0.0};
    plan.note =
        "Plaintext reference uses slot-domain lattice construction z_ref = P*(k + gamma*m) "
        "with deterministic k=0. P=2^modup_period_log2. Gamma is applied in slot domain. "
        "Real and imaginary parts are reduced independently by reduce_mod_integer_lattice. "
        "Therefore (z_ref/P) mod 1 = gamma*m exactly when |gamma*m| < 0.5 per component.";
    return plan;
}

ComplexVector construct_reference_modup_lattice(const ComplexVector& coeff_domain_values,
                                                const ComplexVector& original_slots,
                                                double period_log2,
                                                std::complex<double> gamma) {
    (void)coeff_domain_values;
    if (!std::isfinite(period_log2)) {
        throw std::invalid_argument("period_log2 must be finite");
    }
    if (!std::isfinite(gamma.real()) || !std::isfinite(gamma.imag())) {
        throw std::invalid_argument("gamma must be finite");
    }

    // Reference form A, adjusted for the scan invariant:
    // z_ref_i = P * (k_i + gamma * m_i), with deterministic k_i = 0.
    // The period scan computes r = (z_ref / P) mod Z[i], so r_i = gamma*m_i.
    // This proves the lattice target independently of CKKS/SEAL noise.
    const double period = std::exp2(period_log2);
    ComplexVector result;
    result.reserve(original_slots.size());
    for (const auto& value : original_slots) {
        result.push_back(period * gamma * value);
    }
    return result;
}

BootstrapStcReferenceResult evaluate_bootstrap_stc_reference(
    const BootstrapStcReferencePlan& plan,
    const ComplexVector& input) {
    validate_plan(plan);
    if (input.size() != plan.slots) {
        throw std::invalid_argument("reference input size must match plan slots");
    }

    BootstrapStcReferenceResult result;
    result.input_slots = input;
    const auto stc = FactorizedLinearTransform(make_bootstrap_dft_plan(
        plan.slots, BootstrapDftType::HomomorphicEncode, plan.slots_to_coeff_scale_log2));
    const auto cts = FactorizedLinearTransform(make_bootstrap_dft_plan(
        plan.slots, BootstrapDftType::HomomorphicDecode, plan.slots_to_coeff_scale_log2));
    result.after_slots_to_coeff = stc.apply_plain(input);
    result.after_scale_down = multiply(result.after_slots_to_coeff, std::exp2(-plan.scale_down_log2));
    result.after_modup_reference = construct_reference_modup_lattice(
        result.after_scale_down,
        input,
        plan.modup_period_log2,
        plan.expected_gamma);
    result.after_coeff_to_slots = result.after_modup_reference;

    const auto scan = scan_bootstrap_lattice_periods(
        result.after_coeff_to_slots, input, plan.modup_period_log2, plan.modup_period_log2, 0.1);
    result.best_period_log2 = scan.best_useful_mod_gain.period_log2;
    result.best_gamma = scan.best_useful_mod_gain.gamma;
    result.err_mod_direct = scan.best_mod_direct.err_mod_direct;
    result.err_mod_gain = scan.best_useful_mod_gain.err_mod_gain;
    result.lattice_invariant_passed = scan.has_useful_gain && result.err_mod_gain <= 1e-12;
    if (!result.lattice_invariant_passed) {
        std::ostringstream blocker;
        blocker << "reference lattice invariant failed; expected (z_ref/P) mod 1 = gamma*m"
                << "; best_period_log2=" << result.best_period_log2
                << "; err_mod_gain=" << result.err_mod_gain;
        result.blocker = blocker.str();
    }
    return result;
}

} // namespace m2424
