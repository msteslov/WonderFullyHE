#include "m2424/bootstrap_prototype.hpp"
#include "m2424/profiles.hpp"
#include "m2424/seal_adapter.hpp"

#include <cmath>
#include <cstdio>
#include <vector>

int main() {
    constexpr std::size_t slots = 16;
    constexpr double tolerance = 2e-5;

    auto rotation_steps = m2424::BootstrapPrototype::required_rotation_steps(slots);
    auto adapter = m2424::SealAdapter::create(m2424::profiles::boot_ckks());
    adapter.keygen(rotation_steps, true);

    std::vector<double> input;
    input.reserve(slots);
    for (std::size_t i = 0; i < slots; ++i) {
        input.push_back(1e-5 * std::sin(static_cast<double>(i) / 3.0));
    }

    auto encrypted = adapter.encrypt(adapter.encode(input));
    auto lowered = adapter.mul_plain_rescale(encrypted, adapter.encode_scalar_like(1.0, encrypted));

    m2424::BootstrapPrototype prototype(adapter, slots, tolerance);
    auto report = prototype.refresh_cipher_fast(lowered);

    std::printf("bootstrap_cipher_path\n");
    std::printf("profile,slots,tolerance,normalization_factor,rotation_keys\n");
    std::printf("boot_ckks,%zu,%.6e,%.6e,%zu\n",
                report.slots,
                report.tolerance,
                report.normalization_factor,
                rotation_steps.size());
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

    return report.restore_level_criterion ? 0 : 1;
}
