#include "m2424/linear_transform.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <unordered_map>
#include <utility>

namespace m2424 {
namespace {

constexpr double kZeroCoefficientTolerance = 0.0;

bool has_nonzero_coefficient(const std::vector<double>& coefficients) {
    return std::any_of(coefficients.begin(), coefficients.end(), [](double coefficient) {
        return std::fabs(coefficient) > kZeroCoefficientTolerance;
    });
}

Cipher pairwise_add(SealAdapter& adapter, std::vector<Cipher> terms) {
    if (terms.empty()) {
        throw std::invalid_argument("linear transform has no non-zero terms");
    }
    while (terms.size() > 1) {
        std::vector<Cipher> next;
        next.reserve((terms.size() + 1) / 2);
        for (std::size_t i = 0; i < terms.size(); i += 2) {
            if (i + 1 < terms.size()) {
                next.push_back(adapter.add(terms[i], terms[i + 1]));
            } else {
                next.push_back(std::move(terms[i]));
            }
        }
        terms = std::move(next);
    }
    return std::move(terms.front());
}

} // namespace

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

std::vector<int> LinearTransform::rotationSteps() const {
    std::vector<int> steps;
    for (const auto& term : terms_) {
        if (term.rotation != 0 && has_nonzero_coefficient(term.coefficients)) {
            steps.push_back(term.rotation);
        }
    }
    std::sort(steps.begin(), steps.end());
    steps.erase(std::unique(steps.begin(), steps.end()), steps.end());
    return steps;
}

Cipher LinearTransform::apply(SealAdapter& adapter, const Cipher& input) const {
    std::unordered_map<int, Cipher> rotated_cache;
    std::vector<Cipher> weighted_terms;
    weighted_terms.reserve(terms_.size());

    for (const auto& term : terms_) {
        if (!has_nonzero_coefficient(term.coefficients)) {
            continue;
        }

        const Cipher* rotated = nullptr;
        if (term.rotation == 0) {
            rotated = &input;
        } else {
            auto cached = rotated_cache.find(term.rotation);
            if (cached == rotated_cache.end()) {
                cached = rotated_cache.emplace(term.rotation, adapter.rotate(input, term.rotation)).first;
            }
            rotated = &cached->second;
        }

        Plain encoded = term.coefficients.size() == 1
            ? adapter.encodeScalarFor(term.coefficients.front(), *rotated)
            : adapter.encodeFor(term.coefficients, *rotated);
        weighted_terms.push_back(multiplyPlainAndRescale(adapter, *rotated, encoded));
    }

    return pairwise_add(adapter, std::move(weighted_terms));
}

std::vector<int> power_of_two_rotation_steps(std::size_t slotCount) {
    if (slotCount == 0) {
        throw std::invalid_argument("slotCount must be positive");
    }
    std::vector<int> steps;
    for (std::size_t step = 1; step < slotCount; step <<= 1) {
        steps.push_back(static_cast<int>(step));
    }
    return steps;
}

std::vector<int> sum_slots_rotation_steps(std::size_t slotCount) {
    return power_of_two_rotation_steps(slotCount);
}

Cipher sum_slots(SealAdapter& adapter, const Cipher& input, std::size_t slotCount) {
    if (slotCount == 0) {
        throw std::invalid_argument("slotCount must be positive");
    }
    Cipher result = input;
    for (int step : power_of_two_rotation_steps(slotCount)) {
        Cipher rotated = adapter.rotate(result, step);
        result = adapter.add(result, rotated);
    }
    return result;
}

} // namespace m2424
