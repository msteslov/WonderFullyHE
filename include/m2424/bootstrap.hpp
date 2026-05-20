#pragma once

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

    BootstrapReport analyze_depth(const std::vector<double>& input, std::size_t max_steps);
    const std::vector<BootstrapStage>& pipeline() const noexcept;

private:
    SealAdapter* adapter_{};
    std::vector<BootstrapStage> stages_;
};

const char* to_string(BootstrapStageStatus status) noexcept;

} // namespace m2424
