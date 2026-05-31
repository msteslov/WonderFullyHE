#include "m2424/bootstrap.hpp"
#include "m2424/profiles.hpp"
#include "m2424/seal_adapter.hpp"

#include <algorithm>
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

m2424::ComplexVector make_expected(std::size_t slots) {
    m2424::ComplexVector expected;
    expected.reserve(slots);
    for (std::size_t i = 0; i < slots; ++i) {
        expected.push_back({1e-5 * std::sin(static_cast<double>(i) / 4.0), 0.0});
    }
    return expected;
}

double max_error(m2424::SealAdapter& adapter,
                 const m2424::Cipher& cipher,
                 const m2424::ComplexVector& expected) {
    auto actual = adapter.decode_complex(adapter.decrypt(cipher));
    actual.resize(expected.size());
    double result = 0.0;
    for (std::size_t i = 0; i < expected.size(); ++i) {
        result = std::max(result, std::abs(actual[i] - expected[i]));
    }
    return result;
}

} // namespace

int main() {
    constexpr std::size_t slots = 16;
    constexpr double tolerance = 1e-3;
    constexpr double target_error = 1e-9;

    std::vector<int> rotation_steps;
    const double plan_ms = elapsed_ms([&] {
        rotation_steps = m2424::Bootstrapper::scalable_refresh_rotation_steps(slots);
    });

    auto adapter = m2424::SealAdapter::create(m2424::profiles::boot_deep_ckks());
    const double keygen_ms = elapsed_ms([&] {
        adapter.keygen(rotation_steps, true);
    });

    const auto expected = make_expected(slots);
    auto current = adapter.encrypt(adapter.encode_complex(expected));
    while (adapter.info(current).chain_index > 2) {
        current = adapter.mul_plain_rescale(current, adapter.encode_scalar_like(1.0, current));
    }
    const auto before = adapter.info(current);

    m2424::CkksOperationBudget continuation_budget;
    continuation_budget.ciphertext_muls = before.chain_index + 1;

    m2424::Bootstrapper bootstrapper(adapter);
    m2424::BootstrapGuardedRefreshResult guarded;
    double refresh_ms = 0.0;
    refresh_ms = elapsed_ms([&] {
        guarded = bootstrapper.refresh_slots_to_coeffs_first_checked_guarded(
            current,
            expected,
            continuation_budget,
            target_error,
            slots,
            tolerance);
    });

    std::printf("planner,status,required_levels,available_levels,needs_refresh,executed,blocker\n");
    std::printf("refresh_gate,%s,%zu,%zu,%s,%s,%s\n",
                m2424::to_string(guarded.planning.status),
                guarded.planning.required_compute_levels,
                guarded.planning.available_compute_levels,
                guarded.planning.needs_refresh ? "true" : "false",
                guarded.refresh_executed ? "true" : "false",
                guarded.blocker.c_str());

    const auto after = guarded.refresh_executed ? adapter.info(guarded.refresh.result) : before;
    const double error = guarded.refresh_executed
        ? max_error(adapter, guarded.refresh.result, expected)
        : max_error(adapter, current, expected);

    std::printf("summary,profile,slots,tolerance,rotation_keys,plan_ms,keygen_ms,refresh_ms,chain_before,chain_after,scale_before,scale_after,continuation_levels,max_abs_error,status\n");
    std::printf("refresh,boot_deep_ckks,%zu,%.6e,%zu,%.6f,%.6f,%.6f,%zu,%zu,%.6e,%.6e,%zu,%.6e,%s\n",
                slots,
                tolerance,
                rotation_steps.size(),
                plan_ms,
                keygen_ms,
                refresh_ms,
                before.chain_index,
                after.chain_index,
                before.scale,
                after.scale,
                guarded.refresh_executed ? guarded.refresh.continuation_levels : before.chain_index,
                error,
                guarded.refresh_executed && guarded.refresh.preserve_value_criterion ? "PASS" : "SKIPPED");
    return guarded.refresh_executed && guarded.refresh.preserve_value_criterion ? 0 : 1;
}
