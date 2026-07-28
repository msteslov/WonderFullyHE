#include "m2424/accuracy.hpp"
#include "m2424/profiles.hpp"
#include "m2424/seal_adapter.hpp"

#include <cmath>
#include <cstdio>
#include <algorithm>
#include <vector>

int main() {
    const auto profile = m2424::profiles::boot_ckks();
    auto adapter = m2424::SealAdapter::create(profile);
    adapter.generateKeys(false, false);

    constexpr std::size_t payload_size = 16;
    std::vector<double> input;
    input.reserve(payload_size);
    for (std::size_t i = 0; i < payload_size; ++i) {
        input.push_back(0.01 * std::sin(static_cast<double>(i) / 5.0));
    }

    auto encrypted = adapter.encrypt(adapter.encode(input));
    auto one = adapter.encodeScalarFor(1.0, encrypted);
    auto lowered = multiplyPlainAndRescale(adapter, encrypted, one);
    auto raised = adapter.modRaiseToFirst(lowered);

    auto lowered_out = adapter.decode(adapter.decrypt(lowered));
    auto raised_out = adapter.decode(adapter.decrypt(raised));
    lowered_out.resize(std::min(lowered_out.size(), input.size()));
    raised_out.resize(std::min(raised_out.size(), input.size()));
    const auto lowered_accuracy = m2424::compare(input, lowered_out, 1e-5);
    const auto raised_decode_delta = m2424::compare(lowered_out, raised_out, 1e-5);

    const auto initial_info = adapter.info(encrypted);
    const auto lowered_info = adapter.info(lowered);
    const auto raised_info = adapter.info(raised);

    std::printf("stage,chainIndex,coeffModulusSize,scale,max_abs_error,status\n");
    std::printf("initial,%zu,%zu,%.6e,0.000000e+00,PASS\n",
                initial_info.chainIndex,
                initial_info.coeffModulusSize,
                initial_info.scale);
    std::printf("lowered,%zu,%zu,%.6e,%.6e,%s\n",
                lowered_info.chainIndex,
                lowered_info.coeffModulusSize,
                lowered_info.scale,
                lowered_accuracy.max_abs_error,
                lowered_accuracy.ok ? "PASS" : "FAIL");
    const bool structural_ok = raised_info.chainIndex > lowered_info.chainIndex
        && raised_info.coeffModulusSize > lowered_info.coeffModulusSize;
    const bool raised_ok = structural_ok && raised_decode_delta.ok;

    std::printf("raised_structural,%zu,%zu,%.6e,%.6e,%s\n",
                raised_info.chainIndex,
                raised_info.coeffModulusSize,
                raised_info.scale,
                raised_decode_delta.max_abs_error,
                raised_ok ? "PASS" : "FAIL");

    return raised_ok ? 0 : 1;
}
