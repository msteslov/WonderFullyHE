#include "m2424/bootstrap_prototype.hpp"
#include "m2424/eval_mod.hpp"
#include "m2424/profiles.hpp"
#include "m2424/seal_adapter.hpp"

#include <cmath>
#include <cstdio>
#include <vector>

int main() {
    constexpr std::size_t slots = 16;
    constexpr double tolerance = 2e-5;

    std::vector<m2424::Complex> input;
    input.reserve(slots);
    for (std::size_t i = 0; i < slots; ++i) {
        const double x = static_cast<double>(i);
        input.push_back({
            1e-5 * std::sin(x / 3.0),
            1e-5 * std::cos(x / 5.0)
        });
    }

    auto planning_adapter = m2424::SealAdapter::create(m2424::profiles::boot_ckks());
    m2424::BootstrapPrototype planning(planning_adapter, slots, tolerance);
    auto rotation_steps = planning.rotation_steps();

    auto adapter = m2424::SealAdapter::create(m2424::profiles::boot_ckks());
    adapter.keygen(rotation_steps, true);
    m2424::BootstrapPrototype prototype(adapter, slots, tolerance);
    auto report = prototype.refresh_harness(input);

    std::printf("bootstrap_prototype\n");
    std::printf("profile,slots,tolerance,normalization_factor,rotation_keys\n");
    std::printf("boot_ckks,%zu,%.6e,%.6e,%zu\n",
                report.slots,
                report.tolerance,
                report.normalization_factor,
                rotation_steps.size());
    std::printf("stage,status,chain_before,chain_after,scale_before,scale_after,max_abs_error\n");
    for (const auto& stage : report.stages) {
        std::printf("%s,%s,%zu,%zu,%.6e,%.6e,%.6e\n",
                    stage.name.c_str(),
                    stage.status.c_str(),
                    stage.chain_before,
                    stage.chain_after,
                    stage.scale_before,
                    stage.scale_after,
                    stage.max_abs_error);
    }
    std::printf("criterion,status\n");
    std::printf("Dec(c_prime)_approx_Dec(c),%s\n", report.preserve_value_criterion ? "PASS" : "FAIL");
    std::printf("level_after_available,%s\n", report.restore_level_criterion ? "PASS" : "FAIL");

    return report.preserve_value_criterion && report.restore_level_criterion ? 0 : 1;
}
