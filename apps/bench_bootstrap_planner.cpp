#include "m2424/checked_evaluator.hpp"
#include "m2424/bootstrap_plan.hpp"
#include "m2424/profiles.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <numeric>
#include <vector>

namespace {

std::vector<double> make_input(std::size_t n) {
    std::vector<double> values;
    values.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        values.push_back(0.12 + 0.02 * std::sin(static_cast<double>(i) / 3.0));
    }
    return values;
}

void print_plan(const char* label,
                const m2424::CipherInfo& info,
                const m2424::CkksOperationBudget& budget,
                std::size_t slots) {
    const auto plan = m2424::plan_bootstrap_refresh({
        info,
        budget,
        1e-9,
        slots,
        128,
        m2424::ParameterOptimizeFor::Speed,
        0
    });
    std::printf("%s,%zu,%zu,%zu,%zu,%zu,%s,%zu,%zu,%s,%d,%.6e\n",
                label,
                info.chain_index,
                budget.additions,
                budget.ciphertext_muls,
                budget.rotations,
                budget.plaintext_mul_rescales,
                m2424::to_string(plan.status),
                plan.required_compute_levels,
                plan.available_compute_levels,
                plan.needs_refresh ? "true" : "false",
                plan.parameter_plan.selected_work_bits,
                plan.parameter_plan.estimated_abs_error_bound);
}

} // namespace

int main() {
    const std::size_t payload_size = 8;
    auto adapter = m2424::SealAdapter::create(m2424::profiles::balanced_ckks());
    adapter.keygen(m2424::sum_slots_rotation_steps(payload_size), true);

    const auto input = make_input(payload_size);
    auto cipher = adapter.encrypt(adapter.encode(input));
    m2424::CheckedEvaluator checked(adapter, payload_size, 1e-5);

    std::vector<double> add_expected;
    add_expected.reserve(payload_size);
    for (double value : input) {
        add_expected.push_back(value + value);
    }
    auto add_result = checked.add(cipher, cipher, add_expected);

    std::vector<double> mul_expected;
    mul_expected.reserve(payload_size);
    for (double value : add_expected) {
        mul_expected.push_back(value * value);
    }
    auto mul_result = checked.mul(add_result.cipher, add_result.cipher, mul_expected);

    std::vector<double> rotate_expected;
    rotate_expected.reserve(payload_size);
    for (std::size_t i = 0; i < payload_size; ++i) {
        rotate_expected.push_back(mul_expected[(i + 1) % payload_size]);
    }
    auto rotate_result = checked.rotate(mul_result.cipher, 1, rotate_expected);

    std::printf("case,current_chain,additions,ciphertext_muls,rotations,plaintext_mul_rescales,status,required_levels,available_levels,needs_refresh,planned_work_bits,estimated_abs_error_bound\n");
    const auto tracked_plan = checked.plan_refresh_for_tracked_budget(rotate_result.info, 1e-9, 4096);
    std::printf("%s,%zu,%zu,%zu,%zu,%zu,%s,%zu,%zu,%s,%d,%.6e\n",
                "repeat_tracked_block",
                rotate_result.info.chain_index,
                checked.operation_budget().additions,
                checked.operation_budget().ciphertext_muls,
                checked.operation_budget().rotations,
                checked.operation_budget().plaintext_mul_rescales,
                m2424::to_string(tracked_plan.status),
                tracked_plan.required_compute_levels,
                tracked_plan.available_compute_levels,
                tracked_plan.needs_refresh ? "true" : "false",
                tracked_plan.parameter_plan.selected_work_bits,
                tracked_plan.parameter_plan.estimated_abs_error_bound);

    m2424::CkksOperationBudget deeper_budget = checked.operation_budget();
    deeper_budget.ciphertext_muls += 3;
    print_plan("deeper_next_block", rotate_result.info, deeper_budget, 4096);

    return 0;
}
