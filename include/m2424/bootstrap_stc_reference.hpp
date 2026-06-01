#pragma once

#include "m2424/bootstrap_dft.hpp"

#include <complex>
#include <cstddef>
#include <string>

namespace m2424 {

struct BootstrapLatticeScanPoint {
    double period_log2{};
    std::complex<double> gamma{};
    double err_mod_direct{};
    double err_mod_gain{};
    double max_abs_u{};
    double max_abs_r{};
};

struct BootstrapLatticeScanSummary {
    BootstrapLatticeScanPoint best_mod_direct;
    BootstrapLatticeScanPoint best_mod_gain;
    BootstrapLatticeScanPoint best_useful_mod_gain;
    bool has_useful_gain{};
};

struct BootstrapStcReferencePlan {
    std::size_t slots{};
    double message_scale_log2{};
    double slots_to_coeff_scale_log2{};
    double scale_down_log2{};
    double modup_period_log2{};
    std::complex<double> expected_gamma{};
    std::string note;
};

struct BootstrapStcReferenceResult {
    ComplexVector input_slots;
    ComplexVector after_slots_to_coeff;
    ComplexVector after_scale_down;
    ComplexVector after_modup_reference;
    ComplexVector after_coeff_to_slots;
    double best_period_log2{};
    std::complex<double> best_gamma{};
    double err_mod_direct{};
    double err_mod_gain{};
    bool lattice_invariant_passed{};
    std::string blocker;
};

std::complex<double> bootstrap_lattice_best_gain(const ComplexVector& actual,
                                                 const ComplexVector& expected);

ComplexVector reduce_mod_integer_lattice(const ComplexVector& values);

BootstrapLatticeScanSummary scan_bootstrap_lattice_periods(const ComplexVector& z,
                                                           const ComplexVector& expected,
                                                           double period_min_log2,
                                                           double period_max_log2,
                                                           double minimum_useful_gain_abs);

BootstrapStcReferencePlan make_bootstrap_stc_reference_plan(std::size_t slots,
                                                            double message_scale_log2,
                                                            double target_error);

ComplexVector construct_reference_modup_lattice(const ComplexVector& coeff_domain_values,
                                                const ComplexVector& original_slots,
                                                double period_log2,
                                                std::complex<double> gamma);

BootstrapStcReferenceResult evaluate_bootstrap_stc_reference(
    const BootstrapStcReferencePlan& plan,
    const ComplexVector& input);

} // namespace m2424
