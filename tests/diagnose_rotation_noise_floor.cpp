#include "m2424/bootstrap_dft.hpp"
#include "m2424/bootstrap_precision_model.hpp"
#include "m2424/profiles.hpp"
#include "m2424/seal_adapter.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>
#include <numeric>
#include <set>
#include <string>
#include <vector>

namespace {

constexpr std::size_t kSlots = 4;
constexpr double kAmplitude = 1e-5;
constexpr double kPlainScaleLog2 = 60.0;
constexpr double kMeasuredScaleLog2 = 50.0;
constexpr double kMeasuredWorstRotationError = 7.1e-9;
constexpr double kTargetRotationError = 1e-10;

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

m2424::ComplexVector rotate_plain_head(const m2424::ComplexVector& values, int step) {
    m2424::ComplexVector result(values.size());
    for (std::size_t i = 0; i < values.size(); ++i) {
        const int source = static_cast<int>(i) + step;
        if (source >= 0 && source < static_cast<int>(values.size())) {
            result[i] = values[static_cast<std::size_t>(source)];
        } else {
            result[i] = m2424::Complex{0.0, 0.0};
        }
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

double decoded_error(m2424::SealAdapter& adapter,
                     const m2424::Cipher& cipher,
                     const m2424::ComplexVector& expected) {
    return max_error(expected, head(adapter.decode_complex(adapter.decrypt(cipher))));
}

std::vector<int> required_steps() {
    auto slot_to_coeff = m2424::FactorizedLinearTransform(m2424::make_bootstrap_dft_plan(
        kSlots, m2424::BootstrapDftType::HomomorphicEncode, kPlainScaleLog2));
    auto coeff_to_slot = m2424::FactorizedLinearTransform(m2424::make_bootstrap_dft_plan(
        kSlots, m2424::BootstrapDftType::HomomorphicDecode, kPlainScaleLog2));
    std::set<int> steps;
    for (int step : slot_to_coeff.rotation_steps()) {
        steps.insert(step);
        steps.insert(-step);
    }
    for (int step : coeff_to_slot.rotation_steps()) {
        steps.insert(step);
        steps.insert(-step);
    }
    steps.erase(0);
    return {steps.begin(), steps.end()};
}

std::size_t total_bits(const m2424::CkksProfile& profile) {
    return static_cast<std::size_t>(std::accumulate(
        profile.coeff_modulus_bits.begin(), profile.coeff_modulus_bits.end(), 0));
}

double predicted_rotation_error(double scale_log2) {
    return kMeasuredWorstRotationError * std::exp2(kMeasuredScaleLog2 - scale_log2);
}

void run_profile(const char* profile_name, const m2424::CkksProfile& profile) {
    const auto expected = make_expected();
    const auto steps = required_steps();
    auto adapter = m2424::SealAdapter::create(profile);
    adapter.keygen(steps, true);

    const auto encrypted = adapter.encrypt(adapter.encode_complex(expected));
    const auto input_info = adapter.info(encrypted);
    const double baseline_error = decoded_error(adapter, encrypted, expected);

    double worst_rotate_error = 0.0;
    double worst_rotate_back_error = 0.0;
    for (int step : steps) {
        const auto rotated = adapter.rotate(encrypted, step);
        const auto expected_rotated = rotate_plain_head(expected, step);
        const double rotate_error = decoded_error(adapter, rotated, expected_rotated);
        const auto roundtrip = adapter.rotate(rotated, -step);
        const double rotate_back_error = decoded_error(adapter, roundtrip, expected);
        worst_rotate_error = std::max(worst_rotate_error, rotate_error);
        worst_rotate_back_error = std::max(worst_rotate_back_error, rotate_back_error);
        const auto rotated_info = adapter.info(rotated);
        const auto roundtrip_info = adapter.info(roundtrip);
        std::printf("[diagnose_rotation_noise_floor] profile=%s step=%d rotate_chain=%zu rotate_scale_log2=%.6f rotate_error=%.12e roundtrip_chain=%zu roundtrip_scale_log2=%.6f rotate_back_error=%.12e\n",
                    profile_name,
                    step,
                    rotated_info.chain_index,
                    std::log2(rotated_info.scale),
                    rotate_error,
                    roundtrip_info.chain_index,
                    std::log2(roundtrip_info.scale),
                    rotate_back_error);
    }

    const double predicted = predicted_rotation_error(std::log2(input_info.scale));
    const double ratio = predicted > 0.0 ? worst_rotate_error / predicted : 0.0;
    const bool pass = worst_rotate_error <= kTargetRotationError;
    std::printf("[diagnose_rotation_noise_floor] profile=%s total_coeff_modulus_bits=%zu ciphertext_scale_log2=%.6f baseline_error=%.12e worst_rotate_error=%.12e worst_rotate_back_error=%.12e predicted_worst_rotate_error=%.12e measured_predicted_ratio=%.6f status=%s%s\n",
                profile_name,
                total_bits(profile),
                std::log2(input_info.scale),
                baseline_error,
                worst_rotate_error,
                worst_rotate_back_error,
                predicted,
                ratio,
                pass ? "PASS" : "BLOCKED",
                pass ? "" : " blocker=rotation/key-switch noise is not decode-scale-limited in current SEAL setup; current target is blocked by key-switch noise.");

    if (!steps.empty()) {
        const int step = steps.front();
        auto current = encrypted;
        auto expected_current = expected;
        for (std::size_t count : {1UL, 2UL, 4UL}) {
            current = encrypted;
            expected_current = expected;
            for (std::size_t i = 0; i < count; ++i) {
                current = adapter.rotate(current, step);
                expected_current = rotate_plain_head(expected_current, step);
            }
            const auto info = adapter.info(current);
            std::printf("[diagnose_rotation_noise_floor] profile=%s repeated_step=%d count=%zu chain_index=%zu scale_log2=%.6f max_error=%.12e\n",
                        profile_name,
                        step,
                        count,
                        info.chain_index,
                        std::log2(info.scale),
                        decoded_error(adapter, current, expected_current));
        }
    }
}

} // namespace

int main() {
    try {
        const m2424::CalibratedRotationNoiseModel model{
            kMeasuredScaleLog2,
            kMeasuredWorstRotationError,
            3e-11,
            1.5
        };
        std::printf("[diagnose_rotation_noise_floor] required_ciphertext_scale_log2=%.6f measured_scale_log2=%.6f measured_rotation_error=%.12e target_rotation_error=%.12e safety_bits=%.6f\n",
                    m2424::required_ciphertext_scale_log2(model),
                    model.measured_scale_log2,
                    model.measured_rotation_error,
                    model.target_rotation_error,
                    model.safety_bits);
        run_profile("precision_boot_deep_ckks", m2424::profiles::precision_boot_deep_ckks());
        run_profile("precision_boot_ultra_ckks_59", m2424::profiles::precision_boot_ultra_ckks_59());
        return 0;
    } catch (const std::exception& error) {
        std::printf("[diagnose_rotation_noise_floor] FAIL: %s\n", error.what());
        return 1;
    }
}
