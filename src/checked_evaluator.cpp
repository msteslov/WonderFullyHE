#include "m2424/checked_evaluator.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace m2424 {

CheckedEvaluator::CheckedEvaluator(SealAdapter& adapter, std::size_t payload_size, double tolerance)
    : adapter_(adapter), payload_size_(payload_size), tolerance_(tolerance) {
    if (payload_size_ == 0) {
        throw std::invalid_argument("payload_size must be positive");
    }
    if (!std::isfinite(tolerance_) || tolerance_ < 0.0) {
        throw std::invalid_argument("tolerance must be a non-negative finite value");
    }
}

CheckedResult CheckedEvaluator::add(const Cipher& lhs, const Cipher& rhs, const std::vector<double>& expected) {
    budget_builder_.record_add();
    return finalize("add", adapter_.add(lhs, rhs), expected);
}

CheckedResult CheckedEvaluator::sub(const Cipher& lhs, const Cipher& rhs, const std::vector<double>& expected) {
    budget_builder_.record_add();
    return finalize("sub", adapter_.sub(lhs, rhs), expected);
}

CheckedResult CheckedEvaluator::add_plain(const Cipher& lhs, const Plain& rhs, const std::vector<double>& expected) {
    budget_builder_.record_plaintext_add();
    return finalize("add_plain", adapter_.add_plain(lhs, rhs), expected);
}

CheckedResult CheckedEvaluator::sub_plain(const Cipher& lhs, const Plain& rhs, const std::vector<double>& expected) {
    budget_builder_.record_plaintext_add();
    return finalize("sub_plain", adapter_.sub_plain(lhs, rhs), expected);
}

CheckedResult CheckedEvaluator::mul(const Cipher& lhs, const Cipher& rhs, const std::vector<double>& expected) {
    budget_builder_.record_ciphertext_mul();
    return finalize("mul", adapter_.mul_relin_rescale(lhs, rhs), expected);
}

CheckedResult CheckedEvaluator::mul_plain_rescale(const Cipher& lhs,
                                                  const Plain& rhs,
                                                  const std::vector<double>& expected) {
    budget_builder_.record_plaintext_mul_rescale();
    return finalize("mul_plain_rescale", adapter_.mul_plain_rescale(lhs, rhs), expected);
}

CheckedResult CheckedEvaluator::rescale_to_next(const Cipher& input, const std::vector<double>& expected) {
    budget_builder_.record_rescale_to_next();
    return finalize("rescale_to_next", adapter_.rescale_to_next(input), expected);
}

CheckedResult CheckedEvaluator::rotate(const Cipher& input, int steps, const std::vector<double>& expected) {
    budget_builder_.record_rotation();
    return finalize("rotate", adapter_.rotate(input, steps), expected);
}

CheckedResult CheckedEvaluator::sum_slots(const Cipher& input, std::size_t slot_count,
                                          const std::vector<double>& expected) {
    budget_builder_.record_sum_slots(slot_count);
    return finalize("sum_slots", m2424::sum_slots(adapter_, input, slot_count), expected);
}

CheckedResult CheckedEvaluator::linear_transform(const Cipher& input, const LinearTransform& transform,
                                                 const std::vector<double>& expected) {
    budget_builder_.record_linear_transform(transform);
    return finalize("linear_transform", transform.apply(adapter_, input), expected);
}

const CkksOperationBudget& CheckedEvaluator::operation_budget() const noexcept {
    return budget_builder_.budget();
}

void CheckedEvaluator::reset_operation_budget() noexcept {
    budget_builder_.reset();
}

BootstrapRefreshPlanningResult CheckedEvaluator::plan_refresh_for_tracked_budget(
    const CipherInfo& current_info,
    double target_error,
    std::size_t slots,
    int security_bits,
    ParameterOptimizeFor optimize_for,
    std::size_t min_chain_remaining_after_compute) const {
    return plan_bootstrap_refresh({
        current_info,
        budget_builder_.budget(),
        target_error,
        slots,
        security_bits,
        optimize_for,
        min_chain_remaining_after_compute
    });
}

CheckedResult CheckedEvaluator::finalize(std::string operation, Cipher cipher, const std::vector<double>& expected) {
    if (expected.size() != payload_size_) {
        throw std::invalid_argument("expected vector size must match payload_size");
    }
    auto decoded = adapter_.decode(adapter_.decrypt(cipher));
    if (decoded.size() < payload_size_) {
        throw std::runtime_error("decoded vector is shorter than payload_size");
    }
    decoded.resize(payload_size_);
    auto accuracy = compare(expected, decoded, tolerance_);
    auto info = adapter_.info(cipher);
    return CheckedResult{
        std::move(operation),
        std::move(cipher),
        info,
        accuracy,
        accuracy.ok
    };
}

} // namespace m2424
