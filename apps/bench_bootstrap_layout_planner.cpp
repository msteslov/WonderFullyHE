#include "m2424/bootstrap_plan.hpp"

#include <cmath>
#include <cstdio>

namespace {

void print_layout(const char* label, const m2424::BootstrapLayoutPlanningRequest& request) {
    const auto plan = m2424::plan_bootstrap_layout(request);
    std::printf("%s,%s,%s,%zu,%zu,%.0f,%.0f,%.0f,%zu,%zu,%zu,%zu,%zu,%zu,%d,%d,%.6e,%.6e,%.6e,%.6e,%s,%s\n",
                label,
                m2424::to_string(plan.transform_backend),
                m2424::to_string(plan.security_level),
                plan.slots,
                plan.poly_modulus_degree,
                request.max_abs_after_coeff_to_slot_log2,
                plan.period_log2,
                request.transform_output_scale_log2,
                plan.coeff_to_slot_levels,
                plan.normalization_levels,
                plan.scale_squash_levels,
                plan.evalmod_levels,
                plan.slot_to_coeff_levels,
                plan.total_levels,
                plan.total_coeff_modulus_bits,
                plan.security_budget_bits,
                plan.scale_after_normalization_log2,
                plan.scale_after_squash_log2,
                plan.remaining_modulus_before_evalmod_log2,
                plan.first_evalmod_product_scale_log2,
                m2424::to_string(plan.status),
                plan.blocker.c_str());
}

} // namespace

int main() {
    std::printf("case,backend,security,slots,poly_modulus_degree,max_abs_after_coeff_to_slot_log2,period_log2,transform_output_scale_log2,coeff_to_slot_levels,normalization_levels,scale_squash_levels,evalmod_levels,slot_to_coeff_levels,total_levels,total_coeff_modulus_bits,security_budget_bits,scale_after_normalization_log2,scale_after_squash_log2,remaining_modulus_before_evalmod_log2,first_evalmod_product_scale_log2,status,blocker\n");

    print_layout("dense_current_like", {
        16,
        16384,
        m2424::SecurityLevel::TC128,
        m2424::BootstrapTransformBackend::DenseDiagonal,
        252.0,
        40.0,
        160.0,
        60.0,
        2.0,
        1,
        3,
        1,
        1,
        60,
        40,
        60
    });

    print_layout("fft_coeff_to_slot_current_inverse", {
        16,
        32768,
        m2424::SecurityLevel::TC128,
        m2424::BootstrapTransformBackend::FftLike,
        132.0,
        40.0,
        160.0,
        60.0,
        2.0
    });

    print_layout("fft_target_inverse_factored_estimate", {
        16,
        32768,
        m2424::SecurityLevel::TC128,
        m2424::BootstrapTransformBackend::FftLike,
        132.0,
        40.0,
        160.0,
        60.0,
        2.0,
        0,
        3,
        8,
        1,
        60,
        40,
        60
    });

    return 0;
}
