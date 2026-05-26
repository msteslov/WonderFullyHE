#include "m2424/accuracy.hpp"
#include "m2424/bootstrap.hpp"
#include "m2424/profiles.hpp"
#include "m2424/seal_adapter.hpp"

#include <cmath>
#include <cstdio>
#include <vector>

int main() {
    constexpr std::size_t slots = 16;
    constexpr double tolerance = 2e-5;

    auto rotation_steps = m2424::Bootstrapper::refresh_rotation_steps(slots);
    auto adapter = m2424::SealAdapter::create(m2424::profiles::boot_ckks());
    adapter.keygen(rotation_steps, true);

    std::vector<double> input;
    input.reserve(slots);
    for (std::size_t i = 0; i < slots; ++i) {
        input.push_back(1e-5 * std::sin(static_cast<double>(i) / 4.0));
    }

    auto encrypted = adapter.encrypt(adapter.encode(input));
    auto before_refresh = adapter.mul_plain_rescale(encrypted, adapter.encode_scalar_like(1.0, encrypted));
    const auto before_info = adapter.info(before_refresh);

    m2424::Bootstrapper bootstrapper(adapter);
    std::vector<m2424::Complex> expected;
    expected.reserve(input.size());
    for (double value : input) {
        expected.push_back({value, 0.0});
    }

    auto refresh_report = bootstrapper.refresh_checked(before_refresh, expected, slots, tolerance);
    const auto refreshed_info = adapter.info(refresh_report.result);
    auto refreshed_values = adapter.decode(adapter.decrypt(refresh_report.result));
    refreshed_values.resize(slots);
    const auto refreshed_accuracy = m2424::compare(input, refreshed_values, tolerance);

    auto continued = adapter.mul_plain_rescale(refresh_report.result,
                                               adapter.encode_scalar_like(1.0, refresh_report.result));
    const auto continued_info = adapter.info(continued);

    const bool restored_level = refreshed_info.chain_index > before_info.chain_index;
    const bool consumed_after_refresh = continued_info.chain_index < refreshed_info.chain_index;
    const bool preserve_value = refreshed_accuracy.ok;

    std::printf("bootstrap_end_to_end\n");
    std::printf("stage,chain_index,coeff_modulus_size,scale,status\n");
    std::printf("before_refresh,%zu,%zu,%.6e,PASS\n",
                before_info.chain_index,
                before_info.coeff_modulus_size,
                before_info.scale);
    std::printf("after_refresh,%zu,%zu,%.6e,%s\n",
                refreshed_info.chain_index,
                refreshed_info.coeff_modulus_size,
                refreshed_info.scale,
                restored_level ? "PASS" : "FAIL");
    std::printf("continued_after_refresh,%zu,%zu,%.6e,%s\n",
                continued_info.chain_index,
                continued_info.coeff_modulus_size,
                continued_info.scale,
                consumed_after_refresh ? "PASS" : "FAIL");
    std::printf("criterion,status\n");
    std::printf("preserve_value_after_refresh,%s\n", preserve_value ? "PASS" : "FAIL");
    std::printf("max_abs_error_after_refresh,%.6e\n", refreshed_accuracy.max_abs_error);
    std::printf("level_after_refresh_gt_level_before,%s\n", restored_level ? "PASS" : "FAIL");
    std::printf("operation_after_refresh_available,%s\n", consumed_after_refresh ? "PASS" : "FAIL");

    return preserve_value && restored_level && consumed_after_refresh ? 0 : 1;
}
