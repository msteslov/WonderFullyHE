#pragma once

#include "m2424/accuracy.hpp"
#include "m2424/bootstrap_plan.hpp"
#include "m2424/linear_transform.hpp"
#include "m2424/operation_budget.hpp"
#include "m2424/seal_adapter.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace m2424 {

struct CheckedResult {
    std::string operation;
    Cipher cipher;
    CipherInfo info;
    AccuracyReport accuracy;
    bool ok{};
};

class CheckedEvaluator {
public:
    CheckedEvaluator(SealAdapter& adapter, std::size_t payload_size, double tolerance);

    CheckedResult add(const Cipher& lhs, const Cipher& rhs, const std::vector<double>& expected);
    CheckedResult sub(const Cipher& lhs, const Cipher& rhs, const std::vector<double>& expected);
    CheckedResult add_plain(const Cipher& lhs, const Plain& rhs, const std::vector<double>& expected);
    CheckedResult sub_plain(const Cipher& lhs, const Plain& rhs, const std::vector<double>& expected);
    CheckedResult mul(const Cipher& lhs, const Cipher& rhs, const std::vector<double>& expected);
    CheckedResult mul_plain_rescale(const Cipher& lhs, const Plain& rhs, const std::vector<double>& expected);
    CheckedResult rescale_to_next(const Cipher& input, const std::vector<double>& expected);
    CheckedResult rotate(const Cipher& input, int steps, const std::vector<double>& expected);
    CheckedResult sum_slots(const Cipher& input, std::size_t slot_count, const std::vector<double>& expected);
    CheckedResult linear_transform(const Cipher& input, const LinearTransform& transform,
                                   const std::vector<double>& expected);
    const CkksOperationBudget& operation_budget() const noexcept;
    void reset_operation_budget() noexcept;
    BootstrapRefreshPlanningResult plan_refresh_for_tracked_budget(
        const CipherInfo& current_info,
        double target_error,
        std::size_t slots,
        int security_bits = 128,
        ParameterOptimizeFor optimize_for = ParameterOptimizeFor::Speed,
        std::size_t min_chain_remaining_after_compute = 0) const;

private:
    CheckedResult finalize(std::string operation, Cipher cipher, const std::vector<double>& expected);

    SealAdapter& adapter_;
    std::size_t payload_size_{};
    double tolerance_{};
    OperationBudgetBuilder budget_builder_;
};

} // namespace m2424
