#include "m2424/bootstrap_dft.hpp"
#include "m2424/profiles.hpp"
#include "m2424/seal_adapter.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

namespace {

constexpr std::size_t kSlots = 4;
constexpr double kAmplitude = 1e-5;
constexpr double kDftTolerance = 1e-9;
constexpr const char* kProfileName = "precision_boot_deep_ckks";

m2424::ComplexVector make_expected() {
    m2424::ComplexVector expected;
    expected.reserve(kSlots);
    for (std::size_t i = 0; i < kSlots; ++i) {
        const double x = static_cast<double>(i + 1);
        const double real_sign = (i % 2 == 0) ? 1.0 : -1.0;
        const double imag_sign = (i % 3 == 0) ? -1.0 : 1.0;
        expected.push_back({
            real_sign * kAmplitude * (0.25 + 0.1 * std::sin(x)),
            imag_sign * kAmplitude * (0.15 + 0.05 * std::cos(0.5 * x))
        });
    }
    return expected;
}

m2424::ComplexVector head(m2424::ComplexVector values) {
    values.resize(std::min(values.size(), kSlots));
    return values;
}

double max_abs_value(const m2424::ComplexVector& values) {
    double result = 0.0;
    for (const auto& value : values) {
        result = std::max(result, std::abs(value));
    }
    return result;
}

double max_error(const m2424::ComplexVector& expected, const m2424::ComplexVector& actual) {
    double result = 0.0;
    for (std::size_t i = 0; i < kSlots; ++i) {
        result = std::max(result, std::abs(actual[i] - expected[i]));
    }
    return result;
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

double decoded_error(m2424::SealAdapter& adapter,
                     const m2424::Cipher& cipher,
                     const m2424::ComplexVector& expected) {
    const auto actual = head(adapter.decode_complex(adapter.decrypt(cipher)));
    return max_error(expected, actual);
}

void print_cipher_stage(m2424::SealAdapter& adapter,
                        const char* profile_name,
                        int total_coeff_modulus_bits,
                        double plain_scale_log2,
                        const char* stage,
                        const m2424::Cipher& before,
                        const m2424::Cipher& after,
                        const m2424::ComplexVector& expected) {
    const auto before_info = adapter.info(before);
    const auto after_info = adapter.info(after);
    const auto actual = head(adapter.decode_complex(adapter.decrypt(after)));
    std::printf("[diagnose_bootstrap_invariants] profile=%s coeff_modulus_total_bits=%d plain_scale_log2=%.0f stage=%s chain_before=%zu chain_after=%zu scale_before_log2=%.6f scale_after_log2=%.6f consumed_levels=%zu max_abs_value=%.12e max_error=%.12e\n",
                profile_name,
                total_coeff_modulus_bits,
                plain_scale_log2,
                stage,
                before_info.chain_index,
                after_info.chain_index,
                std::log2(before_info.scale),
                std::log2(after_info.scale),
                before_info.chain_index >= after_info.chain_index ? before_info.chain_index - after_info.chain_index : 0,
                max_abs_value(actual),
                max_error(expected, actual));
}

} // namespace

int main() {
    try {
        const auto expected = make_expected();
        const auto profile = m2424::profiles::precision_boot_deep_ckks();
        const auto total_coeff_modulus_bits = std::accumulate(
            profile.coeff_modulus_bits.begin(), profile.coeff_modulus_bits.end(), 0);

        bool pass = false;
        double best_roundtrip_error = std::numeric_limits<double>::infinity();
        double best_plain_scale_log2 = 0.0;
        for (double plain_scale_log2 : {50.0, 55.0, 60.0}) {
            auto slot_to_coeff = m2424::FactorizedLinearTransform(m2424::make_bootstrap_dft_plan(
                kSlots, m2424::BootstrapDftType::HomomorphicEncode, plain_scale_log2));
            auto coeff_to_slot = m2424::FactorizedLinearTransform(m2424::make_bootstrap_dft_plan(
                kSlots, m2424::BootstrapDftType::HomomorphicDecode, plain_scale_log2));
            const auto plain_coeff = slot_to_coeff.apply_plain(expected);

            auto adapter = m2424::SealAdapter::create(profile);
            adapter.keygen(merged_steps(slot_to_coeff, coeff_to_slot), true);

            auto encrypted = adapter.encrypt(adapter.encode_complex(expected));
            print_cipher_stage(adapter,
                               kProfileName,
                               total_coeff_modulus_bits,
                               plain_scale_log2,
                               "baseline_encrypt_decrypt",
                               encrypted,
                               encrypted,
                               expected);

            auto encrypted_coeff = slot_to_coeff.apply(adapter, encrypted);
            print_cipher_stage(adapter,
                               kProfileName,
                               total_coeff_modulus_bits,
                               plain_scale_log2,
                               "slot_to_coeff",
                               encrypted,
                               encrypted_coeff,
                               plain_coeff);

            auto encrypted_roundtrip = coeff_to_slot.apply(adapter, encrypted_coeff);
            print_cipher_stage(adapter,
                               kProfileName,
                               total_coeff_modulus_bits,
                               plain_scale_log2,
                               "slot_to_coeff_to_coeff_to_slot",
                               encrypted_coeff,
                               encrypted_roundtrip,
                               expected);

            const double roundtrip_error = decoded_error(adapter, encrypted_roundtrip, expected);
            if (roundtrip_error < best_roundtrip_error) {
                best_roundtrip_error = roundtrip_error;
                best_plain_scale_log2 = plain_scale_log2;
            }
            pass = pass || roundtrip_error <= kDftTolerance;
        }

        if (!pass) {
            std::printf("[diagnose_bootstrap_invariants] threshold_fail stage=best_dft_roundtrip best_plain_scale_log2=%.0f max_error=%.12e tolerance=%.12e\n",
                        best_plain_scale_log2,
                        best_roundtrip_error,
                        kDftTolerance);
        }
        std::printf("[diagnose_bootstrap_invariants] best_plain_scale_log2=%.0f best_dft_roundtrip_error=%.12e\n",
                    best_plain_scale_log2,
                    best_roundtrip_error);
        std::printf("[diagnose_bootstrap_invariants] %s\n", pass ? "PASS" : "BLOCKED");
    } catch (const std::exception& error) {
        std::printf("[diagnose_bootstrap_invariants] FAIL: %s\n", error.what());
        return 1;
    }
    return 0;
}
