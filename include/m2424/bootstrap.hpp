#pragma once

#include "m2424/seal_adapter.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace m2424 {

enum class BootstrapStageStatus {
    Done,
    Prototype,
    ModelOnly,
    Planned,
    Failed
};

struct BootstrapStage {
    std::string name;
    BootstrapStageStatus status{};
    std::string note;
};

struct BootstrapReport {
    std::size_t successful_multiplications{};
    std::size_t next_exponent{};
    std::string stop_reason;
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
