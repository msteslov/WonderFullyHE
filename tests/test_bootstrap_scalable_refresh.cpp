#include "m2424/bootstrap.hpp"
#include "m2424/profiles.hpp"
#include "m2424/seal_adapter.hpp"

#include <cmath>
#include <cstdio>
#include <exception>
#include <stdexcept>
#include <vector>

namespace {

m2424::ComplexVector make_expected(std::size_t slots, double amplitude) {
    m2424::ComplexVector expected;
    expected.reserve(slots);
    for (std::size_t i = 0; i < slots; ++i) {
        const double x = static_cast<double>(i + 1);
        expected.push_back({amplitude * std::sin(x / 4.0), amplitude * 0.25 * std::cos(x / 3.0)});
    }
    return expected;
}

void require(bool condition, const char* message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

void run_case(std::size_t slots, double amplitude, double tolerance, m2424::EvalModDegree evalmod_degree) {
    const auto rotation_steps = m2424::Bootstrapper::scalable_refresh_rotation_steps(slots);
    require(!rotation_steps.empty(), "refresh must require rotations");

    auto adapter = m2424::SealAdapter::create(m2424::profiles::boot_deep_ckks());
    adapter.keygen(rotation_steps, true);

    const auto expected = make_expected(slots, amplitude);
    auto current = adapter.encrypt(adapter.encode_complex(expected));
    while (adapter.info(current).chain_index > 2) {
        current = adapter.mul_plain_rescale(current, adapter.encode_scalar_like(1.0, current));
    }

    m2424::Bootstrapper bootstrapper(adapter);
    const auto plan = bootstrapper.plan(slots);
    const auto report = bootstrapper.refresh_slots_to_coeffs_first_checked(
        current, expected, slots, tolerance, evalmod_degree);

    m2424::CkksOperationBudget small_budget;
    small_budget.plaintext_mul_rescales = 1;
    const auto skipped_guard = bootstrapper.refresh_slots_to_coeffs_first_checked_guarded(
        current, expected, small_budget, 1e-9, slots, tolerance, evalmod_degree);

    m2424::CkksOperationBudget refresh_budget;
    refresh_budget.ciphertext_muls = adapter.info(current).chain_index + 1;
    const auto executed_guard = bootstrapper.refresh_slots_to_coeffs_first_checked_guarded(
        current, expected, refresh_budget, 1e-9, slots, tolerance, evalmod_degree);

    bool saw_evalmod = false;
    bool evalmod_passed = false;
    bool saw_result = false;
    for (const auto& stage : report.stages) {
        if (stage.name == "eval_mod") {
            saw_evalmod = true;
            evalmod_passed = stage.status == "DIAG" || stage.status == "PASS";
        }
        if (stage.name == "refresh_result") {
            saw_result = true;
        }
        require(stage.status != "FAIL", "checked refresh should not contain failing stages");
    }

    require(plan.circuit_order == m2424::BootstrapCircuitOrder::SlotsToCoeffsFirst,
            "plan should use SlotsToCoeffsFirst");
    require(plan.transform_backend == m2424::BootstrapTransformBackend::FftLike,
            "plan should use FftLike backend");
    require(report.circuit_order == m2424::BootstrapCircuitOrder::SlotsToCoeffsFirst,
            "report should use SlotsToCoeffsFirst");
    require(report.transform_backend == m2424::BootstrapTransformBackend::FftLike,
            "report should use FftLike backend");
    require(report.evalmod_degree == evalmod_degree,
            "scalable refresh should preserve requested EvalMod degree");
    require(report.inside_evalmod_interval, "normalized values should fit EvalMod interval");
    require(report.preserve_value_criterion, "refresh should preserve checked value");
    require(report.continuation_levels >= 5, "refresh should leave continuation levels");
    require(!skipped_guard.refresh_executed, "guard should skip when compute fits");
    require(skipped_guard.planning.status == m2424::BootstrapRefreshPlanningStatus::ComputeFitsWithoutRefresh,
            "skip guard planning status mismatch");
    require(executed_guard.refresh_executed, "guard should execute when levels are insufficient");
    require(executed_guard.planning.status == m2424::BootstrapRefreshPlanningStatus::RefreshRequired,
            "execute guard planning status mismatch");
    require(executed_guard.refresh.evalmod_degree == evalmod_degree,
            "checked guarded refresh should preserve requested EvalMod degree");
    require(executed_guard.refresh.preserve_value_criterion,
            "checked guarded refresh should preserve value");
    require(executed_guard.blocker == "none", "checked guarded refresh should report no blocker after success");
    require(saw_evalmod && evalmod_passed, "EvalMod stage should be present and pass diagnostically");
    require(saw_result, "refresh_result stage should be present");

    std::printf("[test_bootstrap_scalable_refresh] case PASS slots=%zu amplitude=%.3e evalmod=%s rotations=%zu stages=%zu\n",
                slots,
                amplitude,
                m2424::to_string(evalmod_degree),
                rotation_steps.size(),
                report.stages.size());
}

} // namespace

int main() {
    bool ok = true;
    try {
        run_case(4, 1e-5, 1e-3, m2424::EvalModDegree::P3);
        run_case(8, 1e-5, 1e-3, m2424::EvalModDegree::P3);
        run_case(16, 1e-5, 1e-3, m2424::EvalModDegree::P3);
        run_case(4, 5e-5, 1e-3, m2424::EvalModDegree::P3);
        run_case(4, 1e-5, 1e-3, m2424::EvalModDegree::P3DoubleAngle);
    } catch (const std::exception& error) {
        ok = false;
        std::printf("[test_bootstrap_scalable_refresh] FAIL: %s\n", error.what());
    }
    return ok ? 0 : 1;
}
