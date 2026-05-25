#include "m2424/bootstrap_prototype.hpp"
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

std::vector<m2424::Complex> make_input(std::size_t slots) {
    std::vector<m2424::Complex> input;
    input.reserve(slots);
    for (std::size_t i = 0; i < slots; ++i) {
        const double x = static_cast<double>(i);
        input.push_back({
            1e-5 * std::sin(x / 3.0),
            1e-5 * std::cos(x / 5.0)
        });
    }
    return input;
}

const m2424::BootstrapPrototypeStage* find_stage(const m2424::BootstrapPrototypeReport& report,
                                                 const char* name) {
    for (const auto& stage : report.stages) {
        if (stage.name == name) {
            return &stage;
        }
    }
    return nullptr;
}

} // namespace

int main() {
    constexpr std::size_t slots = 16;
    constexpr double tolerance = 2e-5;

    std::vector<int> rotation_steps;
    const double plan_ms = elapsed_ms([&] {
        rotation_steps = m2424::BootstrapPrototype::required_rotation_steps(slots);
    });

    m2424::SealAdapter adapter = m2424::SealAdapter::create(m2424::profiles::boot_ckks());
    const double keygen_ms = elapsed_ms([&] {
        adapter.keygen(rotation_steps, true);
    });

    auto input = make_input(slots);
    m2424::BootstrapPrototype prototype(adapter, slots, tolerance);

    m2424::BootstrapPrototypeReport report;
    const double refresh_ms = elapsed_ms([&] {
        report = prototype.refresh_harness(input);
    });

    const auto* coeff = find_stage(report, "coeff_to_slot");
    const auto* eval = find_stage(report, "eval_mod");
    const auto* slot = find_stage(report, "slot_to_coeff");
    const auto* final = find_stage(report, "refresh_result");

    std::printf("profile,slots,tolerance,rotation_keys,plan_ms,keygen_ms,refresh_ms,coeff_to_slot_ms,eval_mod_ms,slot_to_coeff_ms,final_max_abs_error,chain_after,status\n");
    std::printf("boot_ckks,%zu,%.6e,%zu,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6e,%zu,%s\n",
                report.slots,
                report.tolerance,
                rotation_steps.size(),
                plan_ms,
                keygen_ms,
                refresh_ms,
                coeff ? coeff->duration_ms : 0.0,
                eval ? eval->duration_ms : 0.0,
                slot ? slot->duration_ms : 0.0,
                final ? final->max_abs_error : 0.0,
                final ? final->chain_after : 0,
                report.preserve_value_criterion && report.restore_level_criterion ? "PASS" : "FAIL");

    return report.preserve_value_criterion && report.restore_level_criterion ? 0 : 1;
}
