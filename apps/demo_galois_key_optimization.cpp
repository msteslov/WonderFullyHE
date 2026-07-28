#include "m2424/accuracy.hpp"
#include "m2424/linear_transform.hpp"
#include "m2424/profiles.hpp"
#include "m2424/seal_adapter.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <numeric>
#include <string>
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

std::string join_steps(const std::vector<int>& steps) {
    std::string result;
    for (std::size_t i = 0; i < steps.size(); ++i) {
        if (i != 0) {
            result += " ";
        }
        result += std::to_string(steps[i]);
    }
    return result;
}

std::vector<double> head(const std::vector<double>& values, std::size_t n) {
    return std::vector<double>(values.begin(), values.begin() + std::min(values.size(), n));
}

struct RunResult {
    const char* mode{};
    std::size_t rotation_count{};
    std::size_t galoisKeysSize{};
    double keygen_ms{};
    double linear_ms{};
    double sum_slots_ms{};
    double linear_max_error{};
    double sum_slots_error{};
};

RunResult run_case(const char* mode,
                   const m2424::CkksProfile& profile,
                   const std::vector<int>& rotationSteps,
                   const std::vector<double>& input,
                   std::size_t sum_size) {
    auto adapter = m2424::SealAdapter::create(profile);

    double keygen_ms = 0.0;
    if (rotationSteps.empty()) {
        keygen_ms = elapsed_ms([&] {
            adapter.generateKeys(true, true);
        });
    } else {
        keygen_ms = elapsed_ms([&] {
            adapter.generateKeys(rotationSteps, true);
        });
    }

    auto encrypted = adapter.encrypt(adapter.encode(input));

    m2424::LinearTransform transform({
        {0, {0.5}},
        {1, {0.25}},
        {2, {-0.125}}
    });

    m2424::Cipher linear_result;
    const double linear_ms = elapsed_ms([&] {
        linear_result = transform.apply(adapter, encrypted);
    });

    std::vector<double> linear_ref;
    linear_ref.reserve(input.size());
    for (std::size_t i = 0; i < input.size(); ++i) {
        const double r1 = i + 1 < input.size() ? input[i + 1] : 0.0;
        const double r2 = i + 2 < input.size() ? input[i + 2] : 0.0;
        linear_ref.push_back(0.5 * input[i] + 0.25 * r1 - 0.125 * r2);
    }
    const auto linear_out = head(adapter.decode(adapter.decrypt(linear_result)), input.size());
    const auto linear_accuracy = m2424::compare(linear_ref, linear_out, 1e-3);

    m2424::Cipher sum_result;
    const double sum_slots_ms = elapsed_ms([&] {
        sum_result = m2424::sum_slots(adapter, encrypted, sum_size);
    });
    const double expected_sum = std::accumulate(input.begin(), input.end(), 0.0);
    const auto sum_out = adapter.decode(adapter.decrypt(sum_result));
    const double sum_error = std::abs(expected_sum - sum_out.front());

    return RunResult{
        mode,
        rotationSteps.empty() ? 0 : rotationSteps.size(),
        adapter.galoisKeysSize(),
        keygen_ms,
        linear_ms,
        sum_slots_ms,
        linear_accuracy.max_abs_error,
        sum_error
    };
}

} // namespace

int main() {
    const std::size_t payload_size = 64;
    const std::size_t sum_size = 64;
    const auto profile = m2424::profiles::depth_ckks();

    std::vector<double> input;
    input.reserve(payload_size);
    for (std::size_t i = 0; i < payload_size; ++i) {
        input.push_back(0.1 + 0.2 * std::sin(static_cast<double>(i) / 9.0));
    }

    auto limited_steps = m2424::sum_slots_rotation_steps(sum_size);
    limited_steps.push_back(1);
    limited_steps.push_back(2);
    std::sort(limited_steps.begin(), limited_steps.end());
    limited_steps.erase(std::unique(limited_steps.begin(), limited_steps.end()), limited_steps.end());

    const auto full = run_case("full", profile, {}, input, sum_size);
    const auto limited = run_case("limited", profile, limited_steps, input, sum_size);

    std::printf("rotation_steps_limited,%s\n", join_steps(limited_steps).c_str());
    std::printf("mode,rotation_count,galoisKeysSize,keygen_ms,linear_transform_ms,sum_slots_ms,linear_max_error,sum_slots_error\n");
    for (const auto& row : {full, limited}) {
        std::printf("%s,%zu,%zu,%.6f,%.6f,%.6f,%.6e,%.6e\n",
                    row.mode,
                    row.rotation_count,
                    row.galoisKeysSize,
                    row.keygen_ms,
                    row.linear_ms,
                    row.sum_slots_ms,
                    row.linear_max_error,
                    row.sum_slots_error);
    }
    std::printf("galois_key_size_ratio,%.6f\n",
                static_cast<double>(full.galoisKeysSize) / static_cast<double>(limited.galoisKeysSize));
    return 0;
}
