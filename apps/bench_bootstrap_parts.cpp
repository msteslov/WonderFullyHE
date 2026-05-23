#include "m2424/accuracy.hpp"
#include "m2424/linear_transform.hpp"
#include "m2424/polynomial.hpp"
#include "m2424/seal_adapter.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <numeric>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

template <class Fn>
double elapsed_ms(Fn&& fn) {
    const auto start = Clock::now();
    fn();
    const auto finish = Clock::now();
    return std::chrono::duration<double, std::milli>(finish - start).count();
}

std::vector<double> head(const std::vector<double>& values, std::size_t n) {
    return std::vector<double>(values.begin(), values.begin() + std::min(values.size(), n));
}

void print_row(const char* block,
               double time_ms,
               const m2424::SealAdapter& adapter,
               const m2424::Cipher& cipher,
               const std::vector<double>& expected,
               const std::vector<double>& actual) {
    const auto err = m2424::compare(expected, actual, 1e-3);
    const auto info = adapter.info(cipher);
    std::printf("%s,%.6f,%zu,%zu,%.6e,%.6e,%zu\n",
                block,
                time_ms,
                info.chain_index,
                info.coeff_modulus_size,
                err.max_abs_error,
                err.mean_abs_error,
                adapter.serialized_size(cipher));
}

} // namespace

int main() {
    const std::size_t poly_modulus_degree = 16384;
    const std::size_t payload_size = 64;
    const std::size_t sum_size = 64;
    m2424::CkksProfile profile{poly_modulus_degree, {60, 40, 40, 40, 40, 60}, std::pow(2.0, 40), poly_modulus_degree / 2};

    auto adapter = m2424::SealAdapter::create(profile);
    auto rotation_steps = m2424::sum_slots_rotation_steps(sum_size);
    rotation_steps.push_back(1);
    rotation_steps.push_back(2);
    std::sort(rotation_steps.begin(), rotation_steps.end());
    rotation_steps.erase(std::unique(rotation_steps.begin(), rotation_steps.end()), rotation_steps.end());
    adapter.keygen(rotation_steps, true);

    std::vector<double> input;
    input.reserve(payload_size);
    for (std::size_t i = 0; i < payload_size; ++i) {
        input.push_back(0.1 + 0.2 * std::sin(static_cast<double>(i) / 9.0));
    }

    auto encrypted = adapter.encrypt(adapter.encode(input));

    std::printf("block,time_ms,chain_index,coeff_modulus_size,max_abs_error,mean_abs_error,serialized_bytes\n");

    m2424::Cipher mul_plain_result;
    const auto mul_plain_ms = elapsed_ms([&] {
        auto coeff = adapter.encode_scalar_like(1.5, encrypted);
        mul_plain_result = adapter.mul_plain_rescale(encrypted, coeff);
    });
    std::vector<double> mul_plain_ref;
    mul_plain_ref.reserve(payload_size);
    for (double x : input) {
        mul_plain_ref.push_back(1.5 * x);
    }
    print_row("mul_plain_rescale", mul_plain_ms, adapter, mul_plain_result, mul_plain_ref,
              head(adapter.decode(adapter.decrypt(mul_plain_result)), payload_size));

    m2424::LinearTransform transform({
        {0, {0.5}},
        {1, {0.25}},
        {2, {-0.125}}
    });
    m2424::Cipher linear_result;
    const auto linear_ms = elapsed_ms([&] {
        linear_result = transform.apply(adapter, encrypted);
    });
    std::vector<double> linear_ref;
    linear_ref.reserve(payload_size);
    for (std::size_t i = 0; i < payload_size; ++i) {
        const double r1 = i + 1 < payload_size ? input[i + 1] : 0.0;
        const double r2 = i + 2 < payload_size ? input[i + 2] : 0.0;
        linear_ref.push_back(0.5 * input[i]
            + 0.25 * r1
            - 0.125 * r2);
    }
    print_row("linear_transform", linear_ms, adapter, linear_result, linear_ref,
              head(adapter.decode(adapter.decrypt(linear_result)), payload_size));

    m2424::Cipher sum_result;
    const auto sum_ms = elapsed_ms([&] {
        sum_result = m2424::sum_slots(adapter, encrypted, sum_size);
    });
    const double total = std::accumulate(input.begin(), input.end(), 0.0);
    const auto sum_actual = head(adapter.decode(adapter.decrypt(sum_result)), 1);
    print_row("sum_slots", sum_ms, adapter, sum_result, std::vector<double>{total}, sum_actual);

    m2424::PolynomialEvaluator polynomial({
        {1, 0.75},
        {3, -0.125}
    });
    m2424::Cipher polynomial_result;
    const auto polynomial_ms = elapsed_ms([&] {
        polynomial_result = polynomial.evaluate(adapter, encrypted);
    });
    std::vector<double> polynomial_ref;
    polynomial_ref.reserve(payload_size);
    for (double x : input) {
        polynomial_ref.push_back(0.75 * x - 0.125 * x * x * x);
    }
    print_row("polynomial_eval", polynomial_ms, adapter, polynomial_result, polynomial_ref,
              head(adapter.decode(adapter.decrypt(polynomial_result)), payload_size));

    return 0;
}
