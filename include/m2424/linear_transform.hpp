#pragma once

#include "m2424/seal_adapter.hpp"

#include <vector>

namespace m2424 {

struct LinearTerm {
    int rotation{};
    std::vector<double> coefficients;
};

class LinearTransform {
public:
    explicit LinearTransform(std::vector<LinearTerm> terms);

    const std::vector<LinearTerm>& terms() const noexcept;
    std::vector<int> rotation_steps() const;
    Cipher apply(SealAdapter& adapter, const Cipher& input) const;

private:
    std::vector<LinearTerm> terms_;
};

Cipher sum_slots(SealAdapter& adapter, const Cipher& input, std::size_t slot_count);
std::vector<int> power_of_two_rotation_steps(std::size_t slot_count);
std::vector<int> sum_slots_rotation_steps(std::size_t slot_count);

} // namespace m2424
