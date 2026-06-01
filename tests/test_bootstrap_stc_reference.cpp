#include "m2424/bootstrap_stc_reference.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>

namespace {

m2424::ComplexVector make_input(bool complex_mode) {
    constexpr std::size_t slots = 4;
    constexpr double amplitude = 1e-5;
    m2424::ComplexVector input;
    input.reserve(slots);
    for (std::size_t i = 0; i < slots; ++i) {
        const double x = static_cast<double>(i + 1);
        const double real_sign = (i % 2 == 0) ? 1.0 : -1.0;
        const double imag_sign = (i % 3 == 0) ? -1.0 : 1.0;
        input.push_back({
            real_sign * amplitude * (0.25 + 0.1 * std::sin(x)),
            complex_mode ? imag_sign * amplitude * (0.15 + 0.05 * std::cos(0.5 * x)) : 0.0
        });
    }
    return input;
}

double max_error(const m2424::ComplexVector& lhs, const m2424::ComplexVector& rhs) {
    double result = 0.0;
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        result = std::max(result, std::abs(lhs[i] - rhs[i]));
    }
    return result;
}

} // namespace

int main() {
    bool ok = true;
    const auto real_input = make_input(false);
    const auto complex_input = make_input(true);
    const auto plan = m2424::make_bootstrap_stc_reference_plan(4, 59.0, 1e-9);

    const auto stc = m2424::FactorizedLinearTransform(m2424::make_bootstrap_dft_plan(
        4, m2424::BootstrapDftType::HomomorphicEncode, 59.0));
    const auto cts = m2424::FactorizedLinearTransform(m2424::make_bootstrap_dft_plan(
        4, m2424::BootstrapDftType::HomomorphicDecode, 59.0));
    const double roundtrip_error = max_error(real_input, cts.apply_plain(stc.apply_plain(real_input)));
    ok = ok && roundtrip_error <= 1e-12;

    const auto real_result = m2424::evaluate_bootstrap_stc_reference(plan, real_input);
    ok = ok && real_result.lattice_invariant_passed;
    ok = ok && real_result.err_mod_gain <= 1e-12;
    ok = ok && std::abs(real_result.best_gamma) >= 0.1;

    const auto complex_result = m2424::evaluate_bootstrap_stc_reference(plan, complex_input);
    ok = ok && complex_result.lattice_invariant_passed;
    ok = ok && complex_result.err_mod_gain <= 1e-12;
    ok = ok && std::abs(complex_result.best_gamma) >= 0.1;

    if (!ok) {
        std::cerr
            << "STC reference test failed"
            << "; roundtrip_error=" << roundtrip_error
            << "; real_err_mod_gain=" << real_result.err_mod_gain
            << "; real_gamma_abs=" << std::abs(real_result.best_gamma)
            << "; real_blocker=" << real_result.blocker
            << "; complex_err_mod_gain=" << complex_result.err_mod_gain
            << "; complex_gamma_abs=" << std::abs(complex_result.best_gamma)
            << "; complex_blocker=" << complex_result.blocker
            << "\n";
        return 1;
    }
    std::cout
        << "STC reference test passed"
        << "; roundtrip_error=" << roundtrip_error
        << "; real_err_mod_gain=" << real_result.err_mod_gain
        << "; real_gamma_abs=" << std::abs(real_result.best_gamma)
        << "; complex_err_mod_gain=" << complex_result.err_mod_gain
        << "; complex_gamma_abs=" << std::abs(complex_result.best_gamma)
        << "\n";
    return 0;
}
