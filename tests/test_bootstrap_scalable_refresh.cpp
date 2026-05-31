#include "m2424/bootstrap.hpp"
#include "m2424/profiles.hpp"
#include "m2424/seal_adapter.hpp"

#include <cmath>
#include <cstdio>
#include <vector>

namespace {

m2424::ComplexVector make_expected(std::size_t slots) {
    m2424::ComplexVector expected;
    expected.reserve(slots);
    for (std::size_t i = 0; i < slots; ++i) {
        expected.push_back({1e-5 * std::sin(static_cast<double>(i) / 4.0), 0.0});
    }
    return expected;
}

} // namespace

int main() {
    constexpr std::size_t slots = 4;
    constexpr double tolerance = 1e-3;

    const auto rotation_steps = m2424::Bootstrapper::scalable_refresh_rotation_steps(slots);
    auto adapter = m2424::SealAdapter::create(m2424::profiles::boot_deep_ckks());
    adapter.keygen(rotation_steps, true);

    const auto expected = make_expected(slots);
    auto current = adapter.encrypt(adapter.encode_complex(expected));
    while (adapter.info(current).chain_index > 2) {
        current = adapter.mul_plain_rescale(current, adapter.encode_scalar_like(1.0, current));
    }

    m2424::Bootstrapper bootstrapper(adapter);
    const auto plan = bootstrapper.plan(slots);
    const auto report = bootstrapper.refresh_slots_to_coeffs_first_checked(
        current, expected, slots, tolerance);

    bool saw_evalmod = false;
    bool evalmod_passed = false;
    for (const auto& stage : report.stages) {
        if (stage.name == "eval_mod") {
            saw_evalmod = true;
            evalmod_passed = stage.status == "DIAG" || stage.status == "PASS";
        }
    }

    const bool ok =
        !rotation_steps.empty()
        && plan.circuit_order == m2424::BootstrapCircuitOrder::SlotsToCoeffsFirst
        && plan.transform_backend == m2424::BootstrapTransformBackend::FftLike
        && report.circuit_order == m2424::BootstrapCircuitOrder::SlotsToCoeffsFirst
        && report.transform_backend == m2424::BootstrapTransformBackend::FftLike
        && report.evalmod_degree == m2424::EvalModDegree::P3
        && report.inside_evalmod_interval
        && report.preserve_value_criterion
        && saw_evalmod
        && evalmod_passed;

    std::printf("[test_bootstrap_scalable_refresh] %s rotations=%zu stages=%zu\n",
                ok ? "PASS" : "FAIL",
                rotation_steps.size(),
                report.stages.size());
    return ok ? 0 : 1;
}
