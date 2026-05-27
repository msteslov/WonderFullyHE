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

class Bootstrapper {
public:
    explicit Bootstrapper(SealAdapter& adapter);

    static std::vector<int> refresh_rotation_steps(std::size_t slots);

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
    BootstrapPipelinePlan plan(std::size_t slots) const;
    const std::vector<BootstrapStage>& pipeline() const noexcept;

private:
    SealAdapter* adapter_{};
    std::vector<BootstrapStage> stages_;
};

const char* to_string(BootstrapStageStatus status) noexcept;

} // namespace m2424
