#include "m2424/bootstrap_precision_model.hpp"
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
constexpr double kTargetRoundtripError = 2e-10;
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

double max_error(const m2424::ComplexVector& expected, const m2424::ComplexVector& actual) {
    double result = 0.0;
    for (std::size_t i = 0; i < kSlots; ++i) {
        result = std::max(result, std::abs(actual[i] - expected[i]));
    }
    return result;
}

double decoded_error(m2424::SealAdapter& adapter,
                     const m2424::Cipher& cipher,
                     const m2424::ComplexVector& expected) {
    return max_error(expected, head(adapter.decode_complex(adapter.decrypt(cipher))));
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

void print_steps(const char* label, const std::vector<int>& steps) {
    std::printf("[diagnose_dft_precision_model] %s=", label);
    for (std::size_t i = 0; i < steps.size(); ++i) {
        std::printf("%s%d", i == 0 ? "" : ",", steps[i]);
    }
    std::printf("\n");
}

void print_cost(const char* label, const m2424::BootstrapDftCost& cost) {
    std::printf("[diagnose_dft_precision_model] cost=%s layers=%zu diagonal_terms=%zu rotations=%zu plaintext_muls=%zu rescales=%zu\n",
                label,
                cost.layer_count,
                cost.diagonal_term_count,
                cost.rotation_count,
                cost.plaintext_multiplication_count,
                cost.rescale_count);
}

} // namespace

int main() {
    try {
        const auto expected = make_expected();
        const auto profile = m2424::profiles::precision_boot_deep_ckks();
        const auto total_bits = std::accumulate(profile.coeff_modulus_bits.begin(),
                                                profile.coeff_modulus_bits.end(),
                                                0);
        (void)total_bits;

        const auto run_mode = [&](const char* mode, bool small_slots4) {
        std::vector<m2424::DftPrecisionMeasurement> measurements;
        for (double plain_scale_log2 : {50.0, 55.0, 60.0}) {
            auto slot_to_coeff = m2424::FactorizedLinearTransform(
                small_slots4
                    ? m2424::make_small_slots4_stc_plan(plain_scale_log2)
                    : m2424::make_bootstrap_dft_plan(
                        kSlots, m2424::BootstrapDftType::HomomorphicEncode, plain_scale_log2));
            auto coeff_to_slot = m2424::FactorizedLinearTransform(
                small_slots4
                    ? m2424::make_small_slots4_cts_plan(plain_scale_log2)
                    : m2424::make_bootstrap_dft_plan(
                        kSlots, m2424::BootstrapDftType::HomomorphicDecode, plain_scale_log2));
            const auto plain_coeff = slot_to_coeff.apply_plain(expected);
            const auto rotations = merged_steps(slot_to_coeff, coeff_to_slot);

            if (plain_scale_log2 == 50.0) {
                std::printf("[diagnose_dft_precision_model] mode=%s\n", mode);
                print_steps("slot_to_coeff_rotation_steps", slot_to_coeff.rotation_steps());
                print_steps("coeff_to_slot_rotation_steps", coeff_to_slot.rotation_steps());
                if (!small_slots4) {
                    print_cost("slot_to_coeff",
                               m2424::estimate_bootstrap_dft_cost(kSlots,
                                                                  m2424::BootstrapDftType::HomomorphicEncode,
                                                                  plain_scale_log2));
                    print_cost("coeff_to_slot",
                               m2424::estimate_bootstrap_dft_cost(kSlots,
                                                                  m2424::BootstrapDftType::HomomorphicDecode,
                                                                  plain_scale_log2));
                }
            }

            auto adapter = m2424::SealAdapter::create(profile);
            adapter.keygen(rotations, true);

            const auto encrypted = adapter.encrypt(adapter.encode_complex(expected));
            const auto input_info = adapter.info(encrypted);
            const double baseline_error = decoded_error(adapter, encrypted, expected);

            const auto encrypted_coeff = slot_to_coeff.apply(adapter, encrypted);
            const auto coeff_info = adapter.info(encrypted_coeff);
            const double stc_error = decoded_error(adapter, encrypted_coeff, plain_coeff);

            const auto encrypted_roundtrip = coeff_to_slot.apply(adapter, encrypted_coeff);
            const auto output_info = adapter.info(encrypted_roundtrip);
            const double roundtrip_error = decoded_error(adapter, encrypted_roundtrip, expected);

            m2424::DftPrecisionMeasurement measurement;
            measurement.profile_name = kProfileName;
            measurement.slots = kSlots;
            measurement.ciphertext_scale_log2 = std::log2(input_info.scale);
            measurement.transform_plain_scale_log2 = plain_scale_log2;
            measurement.rotation_count = rotations.size();
            measurement.input_chain_index = input_info.chain_index;
            measurement.output_chain_index = output_info.chain_index;
            measurement.consumed_levels = input_info.chain_index >= output_info.chain_index
                ? input_info.chain_index - output_info.chain_index
                : 0;
            measurement.baseline_error = baseline_error;
            measurement.slot_to_coeff_error = stc_error;
            measurement.roundtrip_error = roundtrip_error;
            measurements.push_back(measurement);

            std::printf("[diagnose_dft_precision_model] mode=%s profile=%s slots=%zu ciphertext_scale_log2=%.6f transform_plain_scale_log2=%.0f rotation_count=%zu baseline_error=%.12e slot_to_coeff_error=%.12e roundtrip_error=%.12e input_chain=%zu output_chain=%zu consumed_levels=%zu status=%s\n",
                        mode,
                        measurement.profile_name.c_str(),
                        measurement.slots,
                        measurement.ciphertext_scale_log2,
                        measurement.transform_plain_scale_log2,
                        measurement.rotation_count,
                        measurement.baseline_error,
                        measurement.slot_to_coeff_error,
                        measurement.roundtrip_error,
                        measurement.input_chain_index,
                        measurement.output_chain_index,
                        measurement.consumed_levels,
                        measurement.roundtrip_error <= kTargetRoundtripError ? "PASS" : "BLOCKED");
            (void)coeff_info;
        }

        const auto fit = m2424::fit_dft_precision_floor(measurements, kTargetRoundtripError);
        std::printf("[diagnose_dft_precision_model] mode=%s fit A=%.12e noise_floor=%.12e predicted_best_error=%.12e best_plain_scale_log2=%.0f target_roundtrip_error=%.12e floor_blocks_target=%s blocker=%s\n",
                    mode,
                    fit.quantization_coefficient,
                    fit.noise_floor,
                    fit.predicted_best_error,
                    fit.best_plain_scale_log2,
                    kTargetRoundtripError,
                    fit.floor_blocks_target ? "true" : "false",
                    fit.blocker.c_str());
        };

        run_mode("standard", false);
        run_mode("small_slots4", true);
        return 0;
    } catch (const std::exception& error) {
        std::printf("[diagnose_dft_precision_model] FAIL: %s\n", error.what());
        return 1;
    }
}
