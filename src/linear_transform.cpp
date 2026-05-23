#include "m2424/linear_transform.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace m2424 {

LinearTransform::LinearTransform(std::vector<LinearTerm> terms) : terms_(std::move(terms)) {
    if (terms_.empty()) {
        throw std::invalid_argument("linear transform must contain at least one term");
    }
    for (const auto& term : terms_) {
        if (term.coefficients.empty()) {
            throw std::invalid_argument("linear transform term coefficients must not be empty");
        }
        for (double coefficient : term.coefficients) {
            if (!std::isfinite(coefficient)) {
                throw std::invalid_argument("linear transform coefficients must be finite");
            }
        }
    }
}

const std::vector<LinearTerm>& LinearTransform::terms() const noexcept {
    return terms_;
}

std::vector<int> LinearTransform::rotation_steps() const {
    std::vector<int> steps;
    for (const auto& term : terms_) {
        if (term.rotation != 0) {
            steps.push_back(term.rotation);
        }
    }
    std::sort(steps.begin(), steps.end());
    steps.erase(std::unique(steps.begin(), steps.end()), steps.end());
    return steps;
}

Cipher LinearTransform::apply(SealAdapter& adapter, const Cipher& input) const {
    bool has_result = false;
    Cipher result;

    for (const auto& term : terms_) {
        Cipher rotated = term.rotation == 0 ? input : adapter.rotate(input, term.rotation);
        Plain encoded = term.coefficients.size() == 1
            ? adapter.encode_scalar_like(term.coefficients.front(), rotated)
            : adapter.encode_like(term.coefficients, rotated);
        Cipher weighted = adapter.mul_plain_rescale(rotated, encoded);

        if (!has_result) {
            result = std::move(weighted);
            has_result = true;
        } else {
            result = adapter.add(result, weighted);
        }
    }

    return result;
}

std::vector<int> power_of_two_rotation_steps(std::size_t slot_count) {
    if (slot_count == 0) {
        throw std::invalid_argument("slot_count must be positive");
    }
    std::vector<int> steps;
    for (std::size_t step = 1; step < slot_count; step <<= 1) {
        steps.push_back(static_cast<int>(step));
    }
    return steps;
}

std::vector<int> sum_slots_rotation_steps(std::size_t slot_count) {
    return power_of_two_rotation_steps(slot_count);
}

Cipher sum_slots(SealAdapter& adapter, const Cipher& input, std::size_t slot_count) {
    if (slot_count == 0) {
        throw std::invalid_argument("slot_count must be positive");
    }
    Cipher result = input;
    for (int step : power_of_two_rotation_steps(slot_count)) {
        Cipher rotated = adapter.rotate(result, step);
        result = adapter.add(result, rotated);
    }
    return result;
}

} // namespace m2424
