#pragma once

#include "m2424/bootstrap_plan.hpp"
#include "m2424/bootstrap_prototype.hpp"
#include "m2424/seal_adapter.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace m2424 {

enum class BootstrapStageStatus {
    Ready,
    PrimitiveReady,
    SpecificationReady,
    Blocked
};

struct BootstrapCipherMetrics {
    bool available{};
    double scale{};
    std::size_t chain_index{};
    std::size_t coeff_modulus_size{};
    std::size_t ciphertext_size{};
    std::size_t serialized_bytes{};
};

struct BootstrapStage {
    std::string name;
    BootstrapStageStatus status{};
    BootstrapCipherMetrics before;
    BootstrapCipherMetrics after;
    std::string note;
};

struct BootstrapReport {
    BootstrapCipherMetrics input;
    BootstrapCipherMetrics depth_boundary;
    std::size_t successful_multiplications{};
    std::size_t next_exponent{};
    std::string stop_reason;
    bool preserve_value_criterion{};
    bool restore_level_criterion{};
    std::vector<BootstrapStage> stages;
};

struct BootstrapGuardedRefreshResult {
    BootstrapRefreshPlanningResult planning;
    bool refresh_executed{};
    std::string blocker;
    BootstrapPrototypeReport refresh;
};

class Bootstrapper {
public:
    explicit Bootstrapper(SealAdapter& adapter);

    static std::vector<int> refresh_rotation_steps(std::size_t slots);
    static std::vector<int> scalable_refresh_rotation_steps(std::size_t slots);

    BootstrapReport analyze_depth(const std::vector<double>& input, std::size_t max_steps);
    BootstrapPrototypeReport refresh(const Cipher& input, std::size_t slots, double tolerance);
    BootstrapPrototypeReport refresh(const Cipher& input,
                                     std::size_t slots,
                                     double tolerance,
                                     double normalization_factor);
    BootstrapPrototypeReport refresh_checked(const Cipher& input,
                                             const ComplexVector& expected,
                                             std::size_t slots,
                                             double tolerance);
    BootstrapPrototypeReport refresh_checked(const Cipher& input,
                                             const ComplexVector& expected,
                                             std::size_t slots,
                                             double tolerance,
                                             double normalization_factor);
    BootstrapPrototypeReport refresh_slots_to_coeffs_first(const Cipher& input,
                                                           std::size_t slots,
                                                           double tolerance);
    BootstrapPrototypeReport refresh_slots_to_coeffs_first_checked(const Cipher& input,
                                                                   const ComplexVector& expected,
                                                                   std::size_t slots,
                                                                   double tolerance);
    BootstrapGuardedRefreshResult refresh_slots_to_coeffs_first_guarded(
        const Cipher& input,
        const CkksOperationBudget& operation_budget,
        double target_error,
        std::size_t slots,
        double tolerance,
        int security_bits = 128,
        ParameterOptimizeFor optimize_for = ParameterOptimizeFor::Speed,
        std::size_t min_chain_remaining_after_compute = 0);
    BootstrapGuardedRefreshResult refresh_slots_to_coeffs_first_checked_guarded(
        const Cipher& input,
        const ComplexVector& expected,
        const CkksOperationBudget& operation_budget,
        double target_error,
        std::size_t slots,
        double tolerance,
        int security_bits = 128,
        ParameterOptimizeFor optimize_for = ParameterOptimizeFor::Speed,
        std::size_t min_chain_remaining_after_compute = 0);
    BootstrapRefreshPlanningResult plan_refresh_for_budget(
        const Cipher& input,
        const CkksOperationBudget& operation_budget,
        double target_error,
        std::size_t slots,
        int security_bits = 128,
        ParameterOptimizeFor optimize_for = ParameterOptimizeFor::Speed,
        std::size_t min_chain_remaining_after_compute = 0) const;
    BootstrapPipelinePlan plan(std::size_t slots) const;
    const std::vector<BootstrapStage>& pipeline() const noexcept;

private:
    SealAdapter* adapter_{};
    std::vector<BootstrapStage> stages_;
};

const char* to_string(BootstrapStageStatus status) noexcept;

} // namespace m2424
