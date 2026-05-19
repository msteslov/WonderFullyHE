#include "m2424/seal_adapter.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <exception>
#include <functional>
#include <string>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

double elapsed_ms(const std::function<void()>& fn) {
    const auto start = Clock::now();
    fn();
    const auto finish = Clock::now();
    return std::chrono::duration<double, std::milli>(finish - start).count();
}

double max_abs_error(const std::vector<double>& a, const std::vector<double>& b) {
    const std::size_t n = std::min(a.size(), b.size());
    double result = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        result = std::max(result, std::fabs(a[i] - b[i]));
    }
    return result;
}

double mean_abs_error(const std::vector<double>& a, const std::vector<double>& b) {
    const std::size_t n = std::min(a.size(), b.size());
    if (n == 0) return 0.0;
    double total = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        total += std::fabs(a[i] - b[i]);
    }
    return total / static_cast<double>(n);
}

std::vector<double> payload_head(const std::vector<double>& values, std::size_t payload_size) {
    return std::vector<double>(values.begin(), values.begin() + std::min(values.size(), payload_size));
}

void print_row(std::size_t step, std::size_t exponent, double time_ms, double max_error,
               double mean_error, std::size_t bytes, const char* status) {
    std::printf("%zu,%zu,%.6f,%.6e,%.6e,%zu,%s\n",
                step, exponent, time_ms, max_error, mean_error, bytes, status);
}

} // namespace

int main() {
    const std::size_t poly_degree = 16384;
    const std::size_t payload_size = 32;
    const std::size_t max_steps = 8;

    m2424::CkksProfile profile{poly_degree, {60, 40, 40, 40, 40, 60}, std::pow(2.0, 40), poly_degree / 2};
    auto adapter = m2424::SealAdapter::create(profile);
    adapter.keygen(true, false);

    std::vector<double> reference;
    reference.reserve(payload_size);
    for (std::size_t i = 0; i < payload_size; ++i) {
        reference.push_back(0.25 + 0.25 * std::sin(static_cast<double>(i) / 7.0));
    }

    auto current = adapter.encrypt(adapter.encode(reference));
    std::size_t exponent = 1;

    std::printf("step,exponent,time_ms,max_abs_error,mean_abs_error,serialized_bytes,status\n");
    auto initial = payload_head(adapter.decode(adapter.decrypt(current)), payload_size);
    print_row(0, exponent, 0.0, max_abs_error(reference, initial), mean_abs_error(reference, initial),
              adapter.serialized_size(current), "ok");

    for (std::size_t step = 1; step <= max_steps; ++step) {
        try {
            double time_ms = elapsed_ms([&] {
                current = adapter.mul_relin_rescale(current, current);
            });
            exponent *= 2;
            for (double& value : reference) {
                value *= value;
            }

            auto decoded = payload_head(adapter.decode(adapter.decrypt(current)), payload_size);
            print_row(step, exponent, time_ms, max_abs_error(reference, decoded), mean_abs_error(reference, decoded),
                      adapter.serialized_size(current), "ok");
        } catch (const std::exception& error) {
            std::printf("%zu,%zu,0.000000,0.000000e+00,0.000000e+00,0,failed:%s\n",
                        step, exponent * 2, error.what());
            return 0;
        }
    }

    return 0;
}
