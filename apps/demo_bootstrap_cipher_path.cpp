#include "m2424/bootstrap.hpp"
#include "m2424/profiles.hpp"
#include "m2424/seal_adapter.hpp"

#include <cmath>
#include <cstdio>
#include <vector>

int main() {
    constexpr std::size_t slots = 16;
    constexpr double tolerance = 2e-5;

    auto rotation_steps = m2424::Bootstrapper::scalable_refresh_rotation_steps(slots);
    auto adapter = m2424::SealAdapter::create(m2424::profiles::boot_deep_ckks());
    adapter.keygen(rotation_steps, true);

    std::vector<double> input;
    input.reserve(slots);
    for (std::size_t i = 0; i < slots; ++i) {
        input.push_back(1e-5 * std::sin(static_cast<double>(i) / 3.0));
    }

    auto encrypted = adapter.encrypt(adapter.encode(input));
    auto lowered = adapter.mul_plain_rescale(encrypted, adapter.encode_scalar_like(1.0, encrypted));
    while (adapter.info(lowered).chain_index > 2) {
        lowered = adapter.mul_plain_rescale(lowered, adapter.encode_scalar_like(1.0, lowered));
    }

    m2424::Bootstrapper bootstrapper(adapter);
    m2424::ComplexVector expected;
    expected.reserve(input.size());
    for (double value : input) {
        expected.push_back({value, 0.0});
    }
    m2424::CkksOperationBudget continuation_budget;
    continuation_budget.ciphertext_muls = adapter.info(lowered).chain_index + 1;
    const auto guarded = bootstrapper.refresh_slots_to_coeffs_first_checked_guarded(
        lowered,
        expected,
        continuation_budget,
        1e-9,
        slots,
        tolerance);
    const auto& report = guarded.refresh;

    std::printf("bootstrap_cipher_path\n");
    std::printf("profile,slots,tolerance,normalization_factor,rotation_keys,planning_status,refresh_executed,continuation_levels,restore_level\n");
    std::printf("boot_deep_ckks,%zu,%.6e,%.6e,%zu,%s,%s,%zu,%s\n",
                report.slots,
                report.tolerance,
                report.normalization_factor,
                rotation_steps.size(),
                m2424::to_string(guarded.planning.status),
                guarded.refresh_executed ? "true" : "false",
                report.continuation_levels,
                report.restore_level_criterion ? "true" : "false");
    std::printf("stage,status,chain_before,chain_after,scale_before,scale_after,duration_ms\n");
    for (const auto& stage : report.stages) {
        std::printf("%s,%s,%zu,%zu,%.6e,%.6e,%.6f\n",
                    stage.name.c_str(),
                    stage.status.c_str(),
                    stage.chain_before,
                    stage.chain_after,
                    stage.scale_before,
                    stage.scale_after,
                    stage.duration_ms);
    }

    return guarded.refresh_executed && report.preserve_value_criterion ? 0 : 1;
}
