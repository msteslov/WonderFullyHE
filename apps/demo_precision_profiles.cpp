#include "m2424/accuracy.hpp"
#include "m2424/profiles.hpp"
#include "m2424/seal_adapter.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <functional>
#include <numeric>
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

std::vector<double> head(std::vector<double> values, std::size_t n) {
    values.resize(std::min(values.size(), n));
    return values;
}

void print_row(const std::string& profile_name, const char* operation,
               double runtime_ms, const m2424::AccuracyReport& accuracy,
               const m2424::CipherInfo& info, std::size_t serialized_bytes) {
    std::printf("%s,%s,%.6f,%.6e,%.6e,%.6e,%zu,%zu,%zu,%zu,%s\n",
                profile_name.c_str(),
                operation,
                runtime_ms,
                accuracy.max_abs_error,
                accuracy.mean_abs_error,
                info.scale,
                info.chain_index,
                info.coeff_modulus_size,
                info.ciphertext_size,
                serialized_bytes,
                accuracy.ok ? "PASS" : "FAIL");
}

} // namespace

int main() {
    constexpr std::size_t payload_size = 64;
    constexpr double tolerance = 1e-6;

    std::vector<double> input;
    input.reserve(payload_size);
    for (std::size_t i = 0; i < payload_size; ++i) {
        input.push_back(0.2 + 0.1 * std::sin(static_cast<double>(i) / 11.0));
    }

    std::vector<double> square_ref;
    square_ref.reserve(payload_size);
    for (double value : input) {
        square_ref.push_back(value * value);
    }

    std::vector<double> add_ref;
    add_ref.reserve(payload_size);
    for (double value : input) {
        add_ref.push_back(value + value);
    }

    std::printf("profile,operation,runtime_ms,max_abs_error,mean_abs_error,scale,chain_index,coeff_modulus_size,ciphertext_size,serialized_bytes,status\n");

    for (const auto& [profile_name, profile] : m2424::profiles::all()) {
        if (profile_name == "fast_demo_ckks") {
            continue;
        }

        auto adapter = m2424::SealAdapter::create(profile);
        adapter.keygen(true, true);

        auto encrypted = adapter.encrypt(adapter.encode(input));

        m2424::Cipher added;
        const double add_ms = elapsed_ms([&] {
            added = adapter.add(encrypted, encrypted);
        });
        const auto added_out = head(adapter.decode(adapter.decrypt(added)), payload_size);
        print_row(profile_name, "add", add_ms, m2424::compare(add_ref, added_out, tolerance),
                  adapter.info(added), adapter.serialized_size(added));

        m2424::Cipher squared;
        const double square_ms = elapsed_ms([&] {
            squared = adapter.mul_relin_rescale(encrypted, encrypted);
        });
        const auto squared_out = head(adapter.decode(adapter.decrypt(squared)), payload_size);
        print_row(profile_name, "mul_square", square_ms, m2424::compare(square_ref, squared_out, tolerance),
                  adapter.info(squared), adapter.serialized_size(squared));
    }

    return 0;
}
