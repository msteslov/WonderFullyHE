#include "m2424/operation_budget.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace m2424 {
namespace {

bool has_nonzero_coefficient(const std::vector<double>& coefficients) {
    return std::any_of(coefficients.begin(), coefficients.end(), [](double value) {
        return std::fabs(value) > 0.0;
    });
}

std::size_t nonzero_term_count(const LinearTransform& transform) {
    std::size_t count = 0;
    for (const auto& term : transform.terms()) {
        if (has_nonzero_coefficient(term.coefficients)) {
            ++count;
        }
    }
    return count;
}

} // namespace

void OperationBudgetBuilder::reset() noexcept {
    budget_ = CkksOperationBudget{};
}

const CkksOperationBudget& OperationBudgetBuilder::budget() const noexcept {
    return budget_;
}

void OperationBudgetBuilder::record_add(std::size_t count) noexcept {
    budget_.additions += count;
}

void OperationBudgetBuilder::record_plaintext_add(std::size_t count) noexcept {
    budget_.plaintext_additions += count;
}

void OperationBudgetBuilder::record_ciphertext_mul(std::size_t count) noexcept {
    budget_.ciphertext_muls += count;
}

void OperationBudgetBuilder::record_plaintext_mul_rescale(std::size_t count) noexcept {
    budget_.plaintext_mul_rescales += count;
}

void OperationBudgetBuilder::record_rescale_to_next(std::size_t count) noexcept {
    budget_.rescaleToNext += count;
}

void OperationBudgetBuilder::record_mod_switch(std::size_t count) noexcept {
    budget_.mod_switches += count;
}

void OperationBudgetBuilder::record_rotation(std::size_t count) noexcept {
    budget_.rotations += count;
}

void OperationBudgetBuilder::record_linear_transform(const LinearTransform& transform) {
    budget_ = merge_operation_budgets(budget_, estimate_linear_transform_budget(transform));
}

void OperationBudgetBuilder::record_sum_slots(std::size_t slotCount) {
    budget_ = merge_operation_budgets(budget_, estimate_sum_slots_budget(slotCount));
}

void OperationBudgetBuilder::record_evalmod_p3(std::size_t count) noexcept {
    budget_.evalmod_p3 += count;
}

void OperationBudgetBuilder::record_bootstrap_refresh(std::size_t count) noexcept {
    budget_.bootstrap_refreshes += count;
}

CkksOperationBudget estimate_linear_transform_budget(const LinearTransform& transform) {
    const std::size_t terms = nonzero_term_count(transform);
    if (terms == 0) {
        return CkksOperationBudget{};
    }
    CkksOperationBudget budget;
    budget.linear_transforms = 1;
    budget.rotations = transform.rotationSteps().size();
    budget.plaintext_mul_rescales = terms;
    budget.additions = terms > 0 ? terms - 1 : 0;
    return budget;
}

CkksOperationBudget estimate_sum_slots_budget(std::size_t slotCount) {
    if (slotCount == 0) {
        throw std::invalid_argument("slotCount must be positive");
    }
    const std::size_t steps = sum_slots_rotation_steps(slotCount).size();
    CkksOperationBudget budget;
    budget.rotations = steps;
    budget.additions = steps;
    return budget;
}

CkksOperationBudget merge_operation_budgets(CkksOperationBudget lhs,
                                            const CkksOperationBudget& rhs) noexcept {
    lhs.additions += rhs.additions;
    lhs.plaintext_additions += rhs.plaintext_additions;
    lhs.ciphertext_muls += rhs.ciphertext_muls;
    lhs.plaintext_mul_rescales += rhs.plaintext_mul_rescales;
    lhs.rescaleToNext += rhs.rescaleToNext;
    lhs.mod_switches += rhs.mod_switches;
    lhs.rotations += rhs.rotations;
    lhs.linear_transforms += rhs.linear_transforms;
    lhs.evalmod_p3 += rhs.evalmod_p3;
    lhs.bootstrap_refreshes += rhs.bootstrap_refreshes;
    return lhs;
}

} // namespace m2424
