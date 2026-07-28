#pragma once

#include "m2424/linear_transform.hpp"
#include "m2424/parameter_planner.hpp"

#include <cstddef>

namespace m2424 {

class OperationBudgetBuilder {
public:
    void reset() noexcept;
    const CkksOperationBudget& budget() const noexcept;

    void record_add(std::size_t count = 1) noexcept;
    void record_plaintext_add(std::size_t count = 1) noexcept;
    void record_ciphertext_mul(std::size_t count = 1) noexcept;
    void record_plaintext_mul_rescale(std::size_t count = 1) noexcept;
    void record_rescale_to_next(std::size_t count = 1) noexcept;
    void record_mod_switch(std::size_t count = 1) noexcept;
    void record_rotation(std::size_t count = 1) noexcept;
    void record_linear_transform(const LinearTransform& transform);
    void record_sum_slots(std::size_t slotCount);
    void record_evalmod_p3(std::size_t count = 1) noexcept;
    void record_bootstrap_refresh(std::size_t count = 1) noexcept;

private:
    CkksOperationBudget budget_{};
};

CkksOperationBudget estimate_linear_transform_budget(const LinearTransform& transform);
CkksOperationBudget estimate_sum_slots_budget(std::size_t slotCount);
CkksOperationBudget merge_operation_budgets(CkksOperationBudget lhs,
                                            const CkksOperationBudget& rhs) noexcept;

} // namespace m2424
