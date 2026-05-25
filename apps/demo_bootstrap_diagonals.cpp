#include "m2424/accuracy.hpp"
#include "m2424/diagonal_transform.hpp"
#include "m2424/profiles.hpp"
#include "m2424/seal_adapter.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

namespace {

double max_complex_error(const m2424::ComplexVector& expected, const m2424::ComplexVector& actual) {
    const std::size_t n = std::min(expected.size(), actual.size());
    double result = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        result = std::max(result, std::abs(expected[i] - actual[i]));
    }
    return result;
}

m2424::ComplexVector head(m2424::ComplexVector values, std::size_t n) {
    values.resize(std::min(values.size(), n));
    return values;
}

} // namespace

int main() {
    const std::size_t slots = 16;
    const double tolerance = 1e-5;

    auto coeff_to_slot_matrix = m2424::canonical_embedding_matrix(slots);
    auto slot_to_coeff_matrix = m2424::invert_matrix(coeff_to_slot_matrix);
    auto coeff_to_slot = m2424::DiagonalLinearTransform::from_matrix(coeff_to_slot_matrix);
    auto slot_to_coeff = m2424::DiagonalLinearTransform::from_matrix(slot_to_coeff_matrix);

    std::vector<int> rotation_steps = coeff_to_slot.rotation_steps();
    auto inverse_steps = slot_to_coeff.rotation_steps();
    rotation_steps.insert(rotation_steps.end(), inverse_steps.begin(), inverse_steps.end());
    std::sort(rotation_steps.begin(), rotation_steps.end());
    rotation_steps.erase(std::unique(rotation_steps.begin(), rotation_steps.end()), rotation_steps.end());

    m2424::ComplexVector input;
    input.reserve(slots);
    for (std::size_t i = 0; i < slots; ++i) {
        const double x = static_cast<double>(i);
        input.push_back({0.01 * std::sin(x / 3.0), 0.01 * std::cos(x / 5.0)});
    }

    const auto expected_slots = coeff_to_slot.apply_plain(input);
    const auto expected_roundtrip = slot_to_coeff.apply_plain(expected_slots);
    const double cpu_roundtrip_error = max_complex_error(input, expected_roundtrip);

    auto adapter = m2424::SealAdapter::create(m2424::profiles::high_precision_ckks());
    adapter.keygen(rotation_steps, true);

    m2424::ComplexVector packed(adapter.slot_count());
    for (std::size_t i = 0; i < packed.size(); ++i) {
        packed[i] = input[i % slots];
    }

    auto encrypted = adapter.encrypt(adapter.encode_complex(packed));
    auto encrypted_slots = coeff_to_slot.apply(adapter, encrypted);
    auto decoded_slots = head(adapter.decode_complex(adapter.decrypt(encrypted_slots)), slots);
    const double coeff_to_slot_error = max_complex_error(expected_slots, decoded_slots);

    m2424::ComplexVector packed_slots(adapter.slot_count());
    for (std::size_t i = 0; i < packed_slots.size(); ++i) {
        packed_slots[i] = expected_slots[i % slots];
    }
    auto encrypted_expected_slots = adapter.encrypt(adapter.encode_complex(packed_slots));
    auto encrypted_coeffs = slot_to_coeff.apply(adapter, encrypted_expected_slots);
    auto decoded_coeffs = head(adapter.decode_complex(adapter.decrypt(encrypted_coeffs)), slots);
    const double slot_to_coeff_error = max_complex_error(input, decoded_coeffs);

    std::printf("metric,value,status\n");
    std::printf("slots,%zu,PASS\n", slots);
    std::printf("coeff_to_slot_terms,%zu,PASS\n", coeff_to_slot.terms().size());
    std::printf("slot_to_coeff_terms,%zu,PASS\n", slot_to_coeff.terms().size());
    std::printf("rotation_keys,%zu,PASS\n", rotation_steps.size());
    std::printf("first_rotation,%d,PASS\n", rotation_steps.empty() ? 0 : rotation_steps.front());
    std::printf("last_rotation,%d,PASS\n", rotation_steps.empty() ? 0 : rotation_steps.back());
    std::printf("cpu_roundtrip_max_error,%.6e,%s\n", cpu_roundtrip_error, cpu_roundtrip_error <= 1e-10 ? "PASS" : "FAIL");
    std::printf("encrypted_coeff_to_slot_max_error,%.6e,%s\n", coeff_to_slot_error, coeff_to_slot_error <= tolerance ? "PASS" : "FAIL");
    std::printf("encrypted_slot_to_coeff_max_error,%.6e,%s\n", slot_to_coeff_error, slot_to_coeff_error <= tolerance ? "PASS" : "FAIL");

    return cpu_roundtrip_error <= 1e-10
        && coeff_to_slot_error <= tolerance
        && slot_to_coeff_error <= tolerance ? 0 : 1;
}
