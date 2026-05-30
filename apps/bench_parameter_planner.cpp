#include "m2424/accuracy.hpp"
#include "m2424/parameter_planner.hpp"
#include "m2424/seal_adapter.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>
#include <vector>

namespace {

std::vector<double> make_input(std::size_t n) {
    std::vector<double> values;
    values.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        const double x = static_cast<double>(i);
        values.push_back(0.18 + 0.11 * std::sin(x / 7.0) + 0.04 * std::cos(x / 13.0));
    }
    return values;
}

std::vector<double> head(const std::vector<double>& values, std::size_t n) {
    return std::vector<double>(values.begin(), values.begin() + std::min(values.size(), n));
}

const char* pass_fail(bool pass) {
    return pass ? "PASS" : "FAIL";
}

void run_case(double target_error, std::size_t depth, std::size_t slots, m2424::ParameterOptimizeFor optimize_for) {
    const std::size_t payload_size = 64;
    try {
        m2424::CkksPlanningRequest request;
        request.target_error = target_error;
        request.multiplicative_depth = depth;
        request.slots = slots;
        request.security_bits = 128;
        request.optimize_for = optimize_for;

        const auto plan = m2424::plan_ckks_parameters(request);
        auto adapter = m2424::SealAdapter::create(plan.profile);
        adapter.keygen(true, false);

        auto reference = make_input(payload_size);
        auto current = adapter.encrypt(adapter.encode(reference));
        for (std::size_t i = 0; i < depth; ++i) {
            current = adapter.mul_relin_rescale(current, current);
            for (double& value : reference) {
                value *= value;
            }
        }

        const auto decoded = head(adapter.decode(adapter.decrypt(current)), payload_size);
        const auto accuracy = m2424::compare(reference, decoded, target_error);
        const auto info = adapter.info(current);

        std::printf("%.1e,%zu,%zu,%s,%zu,%d,%d,%d,%.6e,%s,%.6e,%.6e,%s,%zu,%.6e,%s\n",
                    target_error,
                    depth,
                    slots,
                    m2424::to_string(optimize_for),
                    plan.profile.poly_modulus_degree,
                    plan.selected_work_bits,
                    plan.selected_scale_log2,
                    plan.estimated_precision_bits,
                    plan.estimated_abs_error_bound,
                    pass_fail(plan.passes_target_error),
                    accuracy.max_abs_error,
                    accuracy.mean_abs_error,
                    pass_fail(accuracy.ok),
                    info.chain_index,
                    info.scale,
                    "ok");
    } catch (const std::exception& error) {
        std::printf("%.1e,%zu,%zu,%s,0,0,0,0,0,FAIL,0,0,FAIL,0,0,failed:%s\n",
                    target_error,
                    depth,
                    slots,
                    m2424::to_string(optimize_for),
                    error.what());
    }
}

void print_budget_plan(double target_error,
                       const char* label,
                       const m2424::CkksOperationBudget& budget,
                       std::size_t slots,
                       m2424::ParameterOptimizeFor optimize_for) {
    try {
        m2424::CkksPlanningRequest request;
        request.target_error = target_error;
        request.multiplicative_depth = 1;
        request.slots = slots;
        request.security_bits = 128;
        request.optimize_for = optimize_for;
        request.use_operation_budget = true;
        request.operation_budget = budget;

        const auto plan = m2424::plan_ckks_parameters(request);
        std::printf("budget,%s,%.1e,%zu,%s,%zu,%d,%d,%d,%.6e,%s,%zu,%s\n",
                    label,
                    target_error,
                    slots,
                    m2424::to_string(optimize_for),
                    plan.profile.poly_modulus_degree,
                    plan.selected_work_bits,
                    plan.selected_scale_log2,
                    plan.estimated_precision_bits,
                    plan.estimated_abs_error_bound,
                    pass_fail(plan.passes_target_error),
                    plan.selected_work_levels,
                    "ok");
    } catch (const std::exception& error) {
        std::printf("budget,%s,%.1e,%zu,%s,0,0,0,0,0,FAIL,0,failed:%s\n",
                    label,
                    target_error,
                    slots,
                    m2424::to_string(optimize_for),
                    error.what());
    }
}

} // namespace

int main() {
    std::printf("target_error,depth,slots,optimize_for,poly_modulus_degree,work_bits,scale_log2,estimated_precision_bits,estimated_abs_error_bound,planned_pass,max_abs_error,mean_abs_error,measured_pass,final_chain_index,final_scale,status\n");
    for (std::size_t depth : {1UL, 2UL, 3UL}) {
        run_case(1e-9, depth, 4096, m2424::ParameterOptimizeFor::Speed);
    }
    run_case(1e-9, 2, 4096, m2424::ParameterOptimizeFor::Conservative);

    std::printf("kind,label,target_error,slots,optimize_for,poly_modulus_degree,work_bits,scale_log2,estimated_precision_bits,estimated_abs_error_bound,planned_pass,work_levels,status\n");
    m2424::CkksOperationBudget linear_budget;
    linear_budget.additions = 3;
    linear_budget.ciphertext_muls = 1;
    linear_budget.plaintext_mul_rescales = 1;
    linear_budget.mod_switches = 2;
    linear_budget.rotations = 4;
    linear_budget.linear_transforms = 1;
    print_budget_plan(1e-9, "mul_plus_linear_transform", linear_budget, 4096, m2424::ParameterOptimizeFor::Speed);

    m2424::CkksOperationBudget evalmod_budget;
    evalmod_budget.ciphertext_muls = 1;
    evalmod_budget.evalmod_p3 = 1;
    print_budget_plan(1e-9, "mul_plus_evalmod_p3", evalmod_budget, 4096, m2424::ParameterOptimizeFor::Speed);
    return 0;
}
