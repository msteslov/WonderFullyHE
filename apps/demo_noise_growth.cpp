#include "m2424/accuracy.hpp"
#include "m2424/profiles.hpp"
#include "m2424/seal_adapter.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <exception>
#include <functional>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

double elapsed_ms(const std::function<void()>& fn) {
    const auto start = Clock::now();
    fn();
    const auto finish = Clock::now();
    return std::chrono::duration<double, std::milli>(finish - start).count();
}

std::vector<double> payload_head(const std::vector<double>& values, std::size_t payload_size) {
    return std::vector<double>(values.begin(), values.begin() + std::min(values.size(), payload_size));
}

void print_row(std::size_t step, std::size_t exponent, double time_ms, const m2424::AccuracyReport& accuracy,
               const m2424::CipherInfo& info, std::size_t bytes, const char* status) {
    std::printf("%zu,%zu,%.6f,%.6e,%.6e,%.6e,%zu,%zu,%zu,%s\n",
                step, exponent, time_ms, accuracy.max_abs_error, accuracy.mean_abs_error,
                info.scale, info.chainIndex, info.coeffModulusSize, bytes, status);
}

} // namespace

int main() {
    const std::size_t payload_size = 32;
    const std::size_t max_steps = 8;

    const auto profile = m2424::profiles::depth_ckks();
    auto adapter = m2424::SealAdapter::create(profile);
    adapter.generateKeys(true, false);

    std::vector<double> reference;
    reference.reserve(payload_size);
    for (std::size_t i = 0; i < payload_size; ++i) {
        reference.push_back(0.25 + 0.25 * std::sin(static_cast<double>(i) / 7.0));
    }

    auto current = adapter.encrypt(adapter.encode(reference));
    std::size_t exponent = 1;

    std::printf("step,exponent,time_ms,max_abs_error,mean_abs_error,scale,chainIndex,coeffModulusSize,serialized_bytes,status\n");
    auto initial = payload_head(adapter.decode(adapter.decrypt(current)), payload_size);
    print_row(0, exponent, 0.0, m2424::compare(reference, initial, 1e-5), adapter.info(current),
              adapter.serializedSize(current), "ok");

    for (std::size_t step = 1; step <= max_steps; ++step) {
        try {
            double time_ms = elapsed_ms([&] {
                current = multiplyRelinearizeAndRescale(adapter, current, current);
            });
            exponent *= 2;
            for (double& value : reference) {
                value *= value;
            }

            auto decoded = payload_head(adapter.decode(adapter.decrypt(current)), payload_size);
            print_row(step, exponent, time_ms, m2424::compare(reference, decoded, 1e-5), adapter.info(current),
                      adapter.serializedSize(current), "ok");
        } catch (const std::exception& error) {
            std::printf("%zu,%zu,0.000000,0.000000e+00,0.000000e+00,0.000000e+00,0,0,0,failed:%s\n",
                        step, exponent * 2, error.what());
            return 0;
        }
    }

    return 0;
}
