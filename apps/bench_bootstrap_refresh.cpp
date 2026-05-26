#include "m2424/bootstrap.hpp"
#include "m2424/profiles.hpp"
#include "m2424/seal_adapter.hpp"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <vector>

namespace {

template <typename Fn>
double elapsed_ms(Fn&& fn) {
    const auto started = std::chrono::steady_clock::now();
    fn();
    const auto finished = std::chrono::steady_clock::now();
    return std::chrono::duration<double, std::milli>(finished - started).count();
}

std::vector<double> make_input(std::size_t slots, double amplitude) {
    std::vector<double> input;
    input.reserve(slots);
    for (std::size_t i = 0; i < slots; ++i) {
        const double x = static_cast<double>(i);
        input.push_back(amplitude * std::sin(x / 4.0));
    }
    return input;
}

} // namespace

int main() {
    constexpr std::size_t slots = 16;
    constexpr double tolerance = 2e-5;
    constexpr double amplitude = 1e-5;

    std::vector<int> rotation_steps;
    const double plan_ms = elapsed_ms([&] {
        rotation_steps = m2424::Bootstrapper::refresh_rotation_steps(slots);
    });

    auto adapter = m2424::SealAdapter::create(m2424::profiles::boot_ckks());
    const double keygen_ms = elapsed_ms([&] {
        adapter.keygen(rotation_steps, true);
    });

    const auto input = make_input(slots, amplitude);
    auto encrypted = adapter.encrypt(adapter.encode(input));
    auto lowered = adapter.mul_plain_rescale(encrypted, adapter.encode_scalar_like(1.0, encrypted));

    m2424::Bootstrapper bootstrapper(adapter);
    m2424::BootstrapPrototypeReport report;
    const double refresh_ms = elapsed_ms([&] {
        report = bootstrapper.refresh(lowered, slots, tolerance);
    });

    const auto before = adapter.info(lowered);
    const auto after = adapter.info(report.result);

    std::printf("summary,profile,slots,tolerance,rotation_keys,plan_ms,keygen_ms,refresh_ms,chain_before,chain_after,scale_before,scale_after,normalization_factor,status\n");
    std::printf("refresh,boot_ckks,%zu,%.6e,%zu,%.6f,%.6f,%.6f,%zu,%zu,%.6e,%.6e,%.6e,%s\n",
                report.slots,
                report.tolerance,
                rotation_steps.size(),
                plan_ms,
                keygen_ms,
                refresh_ms,
                before.chain_index,
                after.chain_index,
                before.scale,
                after.scale,
                report.normalization_factor,
                report.restore_level_criterion ? "PASS" : "FAIL");

    std::printf("stage,status,chain_before,chain_after,scale_before,scale_after,max_abs_error,duration_ms\n");
    for (const auto& stage : report.stages) {
        std::printf("%s,%s,%zu,%zu,%.6e,%.6e,%.6e,%.6f\n",
                    stage.name.c_str(),
                    stage.status.c_str(),
                    stage.chain_before,
                    stage.chain_after,
                    stage.scale_before,
                    stage.scale_after,
                    stage.max_abs_error,
                    stage.duration_ms);
    }

    return report.restore_level_criterion ? 0 : 1;
}
