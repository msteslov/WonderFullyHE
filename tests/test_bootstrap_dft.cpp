#include "m2424/bootstrap_dft.hpp"
#include "m2424/profiles.hpp"
#include "m2424/seal_adapter.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>
#include <vector>

namespace {

m2424::ComplexVector make_input(std::size_t slots) {
    m2424::ComplexVector input;
    input.reserve(slots);
    for (std::size_t i = 0; i < slots; ++i) {
        input.push_back({1e-5 * std::sin(static_cast<double>(i) / 4.0),
                         1e-5 * std::cos(static_cast<double>(i) / 5.0)});
    }
    return input;
}

double max_error(const m2424::ComplexVector& expected, const m2424::ComplexVector& actual) {
    double result = 0.0;
    const std::size_t n = std::min(expected.size(), actual.size());
    for (std::size_t i = 0; i < n; ++i) {
        result = std::max(result, std::abs(expected[i] - actual[i]));
    }
    return result;
}

m2424::ComplexVector head(m2424::ComplexVector values, std::size_t n) {
    values.resize(std::min(values.size(), n));
    return values;
}

std::vector<int> merged_steps(const m2424::FactorizedLinearTransform& a,
                              const m2424::FactorizedLinearTransform& b) {
    auto steps = a.rotation_steps();
    const auto b_steps = b.rotation_steps();
    steps.insert(steps.end(), b_steps.begin(), b_steps.end());
    std::sort(steps.begin(), steps.end());
    steps.erase(std::unique(steps.begin(), steps.end()), steps.end());
    return steps;
}

bool cpu_roundtrip_ok() {
    for (std::size_t slots : {4U, 8U, 16U}) {
        auto decode = m2424::FactorizedLinearTransform(m2424::make_bootstrap_dft_plan(
            slots, m2424::BootstrapDftType::HomomorphicDecode, 40.0));
        auto encode = m2424::FactorizedLinearTransform(m2424::make_bootstrap_dft_plan(
            slots, m2424::BootstrapDftType::HomomorphicEncode, 40.0));
        const auto input = make_input(slots);
        const auto roundtrip = decode.apply_plain(encode.apply_plain(input));
        if (max_error(input, roundtrip) > 1e-10) {
            return false;
        }
        if (decode.rotation_steps().empty() || encode.rotation_steps().empty()) {
            return false;
        }
    }
    return true;
}

bool invalid_scale_rejected() {
    try {
        (void)m2424::make_bootstrap_dft_plan(4,
                                             m2424::BootstrapDftType::HomomorphicDecode,
                                             0.0);
    } catch (const std::invalid_argument&) {
        return true;
    }
    return false;
}

bool unsupported_slots_rejected() {
    try {
        (void)m2424::make_bootstrap_dft_plan(32,
                                             m2424::BootstrapDftType::HomomorphicDecode,
                                             40.0);
    } catch (const std::invalid_argument&) {
        return true;
    }
    return false;
}

bool ciphertext_roundtrip_ok() {
    constexpr std::size_t slots = 4;
    auto decode = m2424::FactorizedLinearTransform(m2424::make_bootstrap_dft_plan(
        slots, m2424::BootstrapDftType::HomomorphicDecode, 40.0));
    auto encode = m2424::FactorizedLinearTransform(m2424::make_bootstrap_dft_plan(
        slots, m2424::BootstrapDftType::HomomorphicEncode, 40.0));
    auto adapter = m2424::SealAdapter::create(m2424::profiles::boot_ckks());
    adapter.keygen(merged_steps(decode, encode), true);

    const auto input = make_input(slots);
    auto current = adapter.encrypt(adapter.encode_complex(input));
    current = encode.apply(adapter, current);
    current = decode.apply(adapter, current);
    const auto actual = head(adapter.decode_complex(adapter.decrypt(current)), slots);
    const auto info = adapter.info(current);
    return max_error(input, actual) < 1e-3
        && std::isfinite(info.scale)
        && info.chain_index > 0;
}

} // namespace

int main() {
    bool ok = true;
    ok = cpu_roundtrip_ok() && ok;
    ok = invalid_scale_rejected() && ok;
    ok = unsupported_slots_rejected() && ok;
    ok = ciphertext_roundtrip_ok() && ok;
    std::printf("[test_bootstrap_dft] %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
