#include "m2424/bootstrap.hpp"
#include "m2424/checked_evaluator.hpp"
#include "m2424/operation_budget.hpp"
#include "m2424/profiles.hpp"

#include <cmath>
#include <cstdio>
#include <exception>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace {

void require(bool condition, const char* message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

void test_static_estimators() {
    const auto sum_budget = m2424::estimate_sum_slots_budget(8);
    require(sum_budget.rotations == 3, "sum_slots rotations mismatch");
    require(sum_budget.additions == 3, "sum_slots additions mismatch");

    m2424::LinearTransform transform({
        {0, {0.5}},
        {1, {0.25}},
        {3, {0.125}},
        {4, {0.0}}
    });
    const auto linear_budget = m2424::estimate_linear_transform_budget(transform);
    require(linear_budget.linear_transforms == 1, "linear transform count mismatch");
    require(linear_budget.rotations == 2, "linear transform rotations mismatch");
    require(linear_budget.plaintext_mul_rescales == 3, "linear transform plaintext mul count mismatch");
    require(linear_budget.additions == 2, "linear transform additions mismatch");
}

void test_checked_evaluator_tracking() {
    const std::size_t payload_size = 8;
    auto adapter = m2424::SealAdapter::create(m2424::profiles::balanced_ckks());
    adapter.keygen(m2424::sum_slots_rotation_steps(payload_size), true);

    std::vector<double> input;
    input.reserve(payload_size);
    for (std::size_t i = 0; i < payload_size; ++i) {
        input.push_back(0.15 + 0.01 * static_cast<double>(i));
    }
    auto cipher = adapter.encrypt(adapter.encode(input));
    m2424::CheckedEvaluator checked(adapter, payload_size, 1e-5);

    auto shift_plain = adapter.encode_like(std::vector<double>(payload_size, 0.01), cipher);
    std::vector<double> add_plain_expected;
    add_plain_expected.reserve(payload_size);
    for (double value : input) {
        add_plain_expected.push_back(value + 0.01);
    }
    auto add_plain_result = checked.add_plain(cipher, shift_plain, add_plain_expected);

    std::vector<double> add_expected;
    add_expected.reserve(payload_size);
    for (double value : add_plain_expected) {
        add_expected.push_back(value + value);
    }
    auto add_result = checked.add(add_plain_result.cipher, add_plain_result.cipher, add_expected);

    auto scalar_plain = adapter.encode_scalar_like(1.25, add_result.cipher);
    std::vector<double> mul_plain_expected;
    mul_plain_expected.reserve(payload_size);
    for (double value : add_expected) {
        mul_plain_expected.push_back(1.25 * value);
    }
    auto mul_plain_result = checked.mul_plain_rescale(add_result.cipher, scalar_plain, mul_plain_expected);

    std::vector<double> mul_expected;
    mul_expected.reserve(payload_size);
    for (double value : mul_plain_expected) {
        mul_expected.push_back(value * value);
    }
    auto mul_result = checked.mul(mul_plain_result.cipher, mul_plain_result.cipher, mul_expected);

    std::vector<double> rotate_expected;
    rotate_expected.reserve(payload_size);
    for (std::size_t i = 0; i < payload_size; ++i) {
        rotate_expected.push_back(mul_expected[(i + 1) % payload_size]);
    }
    auto rotate_result = checked.rotate(mul_result.cipher, 1, rotate_expected);

    const double sum_value = std::accumulate(rotate_expected.begin(), rotate_expected.end(), 0.0);
    auto sum_result = checked.sum_slots(rotate_result.cipher, payload_size, std::vector<double>(payload_size, sum_value));
    (void)sum_result;

    const auto& budget = checked.operation_budget();
    require(budget.additions == 4, "tracked additions mismatch");
    require(budget.plaintext_additions == 1, "tracked plaintext additions mismatch");
    require(budget.ciphertext_muls == 1, "tracked ciphertext mul mismatch");
    require(budget.plaintext_mul_rescales == 1, "tracked plaintext mul mismatch");
    require(budget.rotations == 4, "tracked rotations mismatch");

    m2424::CkksPlanningRequest request;
    request.target_error = 1e-9;
    request.multiplicative_depth = 1;
    request.slots = 4096;
    request.use_operation_budget = true;
    request.operation_budget = budget;
    const auto plan = m2424::plan_ckks_parameters(request);
    require(plan.selected_work_levels == 2, "tracked budget level plan mismatch");
    require(plan.passes_target_error, "tracked budget should pass target");

    const auto refresh_plan = checked.plan_refresh_for_tracked_budget(
        adapter.info(sum_result.cipher),
        1e-9,
        4096);
    require(refresh_plan.status == m2424::BootstrapRefreshPlanningStatus::RefreshRequired,
            "checked evaluator refresh planning mismatch");

    m2424::Bootstrapper bootstrapper(adapter);
    const auto bootstrapper_plan = bootstrapper.plan_refresh_for_budget(sum_result.cipher,
                                                                        budget,
                                                                        1e-9,
                                                                        4096);
    require(bootstrapper_plan.status == refresh_plan.status,
            "bootstrapper refresh planning should match checked evaluator");
}

} // namespace

int main() {
    bool ok = true;
    try {
        test_static_estimators();
        test_checked_evaluator_tracking();
    } catch (const std::exception& error) {
        ok = false;
        std::printf("[test_operation_budget] FAIL: %s\n", error.what());
    }
    if (ok) {
        std::printf("[test_operation_budget] PASS\n");
    }
    return ok ? 0 : 1;
}
