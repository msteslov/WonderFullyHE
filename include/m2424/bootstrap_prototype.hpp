#pragma once

#include "m2424/diagonal_transform.hpp"
#include "m2424/seal_adapter.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace m2424 {

struct BootstrapPrototypeStage {
    std::string name;
    std::string status;
    std::size_t chain_before{};
    std::size_t chain_after{};
    double scale_before{};
    double scale_after{};
    double max_abs_error{};
    double duration_ms{};
};

struct BootstrapPrototypeReport {
    std::size_t slots{};
    double tolerance{};
    double normalization_factor{};
    bool preserve_value_criterion{};
    bool restore_level_criterion{};
    std::vector<BootstrapPrototypeStage> stages;
    Cipher result;
};

class BootstrapPrototype {
public:
    BootstrapPrototype(SealAdapter& adapter, std::size_t slots, double tolerance);

    static std::vector<int> required_rotation_steps(std::size_t slots);

    std::vector<int> rotation_steps() const;
    BootstrapPrototypeReport refresh_harness(const ComplexVector& input) const;

private:
    SealAdapter& adapter_;
    std::size_t slots_{};
    double tolerance_{};
    DiagonalLinearTransform coeff_to_slot_;
    DiagonalLinearTransform slot_to_coeff_;
};

} // namespace m2424
