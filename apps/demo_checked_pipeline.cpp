#include "m2424/abft.hpp"
#include "m2424/checked_evaluator.hpp"
#include "m2424/linear_transform.hpp"
#include "m2424/profiles.hpp"
#include "m2424/seal_adapter.hpp"

#include <cmath>
#include <cstdio>
#include <numeric>
#include <string>
#include <vector>

namespace {

void print_checked(const m2424::CheckedResult& result, const char* abft_status) {
    std::printf("%s,%.6e,%.6e,%.6e,%zu,%zu,%zu,%s,%s\n",
                result.operation.c_str(),
                result.accuracy.max_abs_error,
                result.accuracy.mean_abs_error,
                result.info.scale,
                result.info.chain_index,
                result.info.coeff_modulus_size,
                result.info.ciphertext_size,
                abft_status,
                result.ok ? "PASS" : "FAIL");
}

std::vector<double> rotate_left_zero_tail(const std::vector<double>& input, int steps) {
    std::vector<double> result(input.size(), 0.0);
    const auto offset = static_cast<std::size_t>(steps);
    for (std::size_t i = 0; i < input.size(); ++i) {
        if (i + offset < input.size()) {
            result[i] = input[i + offset];
        }
    }
    return result;
}

} // namespace

int main() {
    constexpr std::size_t payload_size = 16;
    constexpr double tolerance = 1e-6;

    auto adapter = m2424::SealAdapter::create(m2424::profiles::high_precision_ckks());
    adapter.keygen(m2424::sum_slots_rotation_steps(payload_size), true);
    m2424::CheckedEvaluator checked(adapter, payload_size, tolerance);

    std::vector<double> input;
    input.reserve(payload_size);
    for (std::size_t i = 0; i < payload_size; ++i) {
        input.push_back(0.25 + 0.05 * std::sin(static_cast<double>(i) / 5.0));
    }

    auto encrypted = adapter.encrypt(adapter.encode(input));

    std::printf("operation,max_abs_error,mean_abs_error,scale,chain_index,coeff_modulus_size,ciphertext_size,abft_status,status\n");

    std::vector<double> add_ref;
    add_ref.reserve(payload_size);
    for (double value : input) {
        add_ref.push_back(value + value);
    }
    auto add_result = checked.add(encrypted, encrypted, add_ref);
    auto add_decoded = adapter.decode(adapter.decrypt(add_result.cipher));
    auto add_abft = m2424::abft::verify_checksum_value(add_decoded, payload_size,
                                                       m2424::abft::checksum(add_ref), tolerance);
    print_checked(add_result, add_abft.ok ? "PASS" : "FAIL");

    std::vector<double> mul_ref;
    mul_ref.reserve(payload_size);
    for (double value : add_ref) {
        mul_ref.push_back(value * value);
    }
    auto mul_result = checked.mul(add_result.cipher, add_result.cipher, mul_ref);
    auto mul_decoded = adapter.decode(adapter.decrypt(mul_result.cipher));
    auto mul_abft = m2424::abft::verify_checksum_value(mul_decoded, payload_size,
                                                       m2424::abft::checksum(mul_ref), tolerance);
    print_checked(mul_result, mul_abft.ok ? "PASS" : "FAIL");

    const auto rotate_ref = rotate_left_zero_tail(mul_ref, 1);
    auto rotate_result = checked.rotate(mul_result.cipher, 1, rotate_ref);
    auto rotate_decoded = adapter.decode(adapter.decrypt(rotate_result.cipher));
    auto rotate_abft = m2424::abft::verify_checksum_value(rotate_decoded, payload_size,
                                                          m2424::abft::checksum(rotate_ref), tolerance);
    print_checked(rotate_result, rotate_abft.ok ? "PASS" : "FAIL");

    const double sum_ref = std::accumulate(rotate_ref.begin(), rotate_ref.end(), 0.0);
    std::vector<double> sum_expected = rotate_ref;
    for (std::size_t step = 1; step < payload_size; step <<= 1) {
        auto rotated = rotate_left_zero_tail(sum_expected, static_cast<int>(step));
        for (std::size_t i = 0; i < sum_expected.size(); ++i) {
            sum_expected[i] += rotated[i];
        }
    }
    auto sum_result = checked.sum_slots(rotate_result.cipher, payload_size, sum_expected);
    auto sum_decoded = adapter.decode(adapter.decrypt(sum_result.cipher));
    const bool sum_abft_ok = std::fabs(sum_decoded.front() - sum_ref) <= tolerance;
    print_checked(sum_result, sum_abft_ok ? "PASS" : "FAIL");

    return add_result.ok && add_abft.ok
        && mul_result.ok && mul_abft.ok
        && rotate_result.ok && rotate_abft.ok
        && sum_result.ok && sum_abft_ok ? 0 : 1;
}
