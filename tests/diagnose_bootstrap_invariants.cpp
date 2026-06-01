#include "m2424/bootstrap_dft.hpp"
#include "m2424/profiles.hpp"
#include "m2424/seal_adapter.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>
#include <numeric>
#include <string>
#include <vector>

namespace {

constexpr std::size_t kSlots = 4;
constexpr double kAmplitude = 1e-5;
constexpr double kBaselineTolerance = 1e-11;
constexpr double kDftTolerance = 1e-9;
constexpr double kPlainScaleLog2 = 50.0;
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

void print_cipher_stage(m2424::SealAdapter& adapter,
                        const char* profile_name,
                        int total_coeff_modulus_bits,
                        const char* stage,
                        const m2424::Cipher& cipher,
                        const m2424::ComplexVector& expected) {
    const auto info = adapter.info(cipher);
    const auto actual = head(adapter.decode_complex(adapter.decrypt(cipher)));
    std::printf("[diagnose_bootstrap_invariants] profile=%s coeff_modulus_total_bits=%d stage=%s chain_index=%zu scale_log2=%.6f max_abs_value=%.12e max_error=%.12e\n",
                profile_name,
                total_coeff_modulus_bits,
                stage,
                info.chain_index,
                std::log2(info.scale),
                max_abs_value(actual),
                max_error(expected, actual));
}

bool report_threshold(const char* stage, double error, double tolerance) {
    if (error > tolerance) {
        std::printf("[diagnose_bootstrap_invariants] threshold_fail stage=%s max_error=%.12e tolerance=%.12e\n",
                    stage,
                    error,
                    tolerance);
        return false;
    }
    return true;
}

} // namespace

int main() {
    bool ok = true;
    try {
        const auto expected = make_expected();
        auto slot_to_coeff = m2424::FactorizedLinearTransform(m2424::make_bootstrap_dft_plan(
            kSlots, m2424::BootstrapDftType::HomomorphicEncode, kPlainScaleLog2));
        auto coeff_to_slot = m2424::FactorizedLinearTransform(m2424::make_bootstrap_dft_plan(
            kSlots, m2424::BootstrapDftType::HomomorphicDecode, kPlainScaleLog2));
        const auto plain_coeff = slot_to_coeff.apply_plain(expected);

        const auto profile = m2424::profiles::precision_boot_deep_ckks();
        const auto total_coeff_modulus_bits = std::accumulate(
            profile.coeff_modulus_bits.begin(), profile.coeff_modulus_bits.end(), 0);

        auto adapter = m2424::SealAdapter::create(profile);
        adapter.keygen(merged_steps(slot_to_coeff, coeff_to_slot), true);

        auto encrypted = adapter.encrypt(adapter.encode_complex(expected));
        print_cipher_stage(adapter, kProfileName, total_coeff_modulus_bits, "baseline_encrypt_decrypt", encrypted, expected);
        bool pass = report_threshold(
            "baseline_encrypt_decrypt",
            max_error(expected, head(adapter.decode_complex(adapter.decrypt(encrypted)))),
            kBaselineTolerance);

        auto encrypted_coeff = slot_to_coeff.apply(adapter, encrypted);
        print_cipher_stage(adapter, kProfileName, total_coeff_modulus_bits, "slot_to_coeff", encrypted_coeff, plain_coeff);
        pass = report_threshold(
            "slot_to_coeff",
            max_error(plain_coeff, head(adapter.decode_complex(adapter.decrypt(encrypted_coeff)))),
            kDftTolerance) && pass;

        auto encrypted_roundtrip = coeff_to_slot.apply(adapter, encrypted_coeff);
        print_cipher_stage(adapter, kProfileName, total_coeff_modulus_bits, "slot_to_coeff_to_coeff_to_slot", encrypted_roundtrip, expected);
        pass = report_threshold(
            "slot_to_coeff_to_coeff_to_slot",
            max_error(expected, head(adapter.decode_complex(adapter.decrypt(encrypted_roundtrip)))),
            kDftTolerance) && pass;

        std::printf("[diagnose_bootstrap_invariants] %s\n", pass ? "PASS" : "FAIL");
        ok = pass;
    } catch (const std::exception& error) {
        ok = false;
        std::printf("[diagnose_bootstrap_invariants] FAIL: %s\n", error.what());
    }
    return ok ? 0 : 1;
}
