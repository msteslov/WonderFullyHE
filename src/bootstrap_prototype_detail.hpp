#pragma once

#include "m2424/bootstrap_prototype.hpp"
#include "m2424/eval_mod.hpp"

#include <chrono>

namespace m2424::bootstrap_prototype_detail {

double max_complex_error(const ComplexVector& expected, const ComplexVector& actual);
ComplexVector head(ComplexVector values, std::size_t n);
ComplexVector scaled(const ComplexVector& values, double factor);
ComplexVector evaluate_plain(const EvalModPolynomial& eval_mod,
                             const ComplexVector& input,
                             EvalModDegree degree);
double max_abs_value(const ComplexVector& values);
double normalization_factor_for(const ComplexVector& values);

BootstrapPrototypeStage make_stage(const std::string& name,
                                   const CipherInfo& before,
                                   const CipherInfo& after,
                                   double max_error,
                                   double tolerance,
                                   double duration_ms,
                                   bool checked = true);
BootstrapPrototypeStage make_harness_stage(const CipherInfo& after);

Cipher apply_output_scale_repair(SealAdapter& adapter, const Cipher& input, double target_scale_log2);

void mark_stage_structural(BootstrapPrototypeStage& stage);
void mark_stage_diagnostic(BootstrapPrototypeStage& stage);

double finite_exp2_or_zero(double exponent);

template <typename Fn>
double elapsed_ms(Fn&& fn) {
    const auto started = std::chrono::steady_clock::now();
    fn();
    const auto finished = std::chrono::steady_clock::now();
    return std::chrono::duration<double, std::milli>(finished - started).count();
}

} // namespace m2424::bootstrap_prototype_detail
