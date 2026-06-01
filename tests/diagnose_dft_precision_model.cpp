#include "m2424/bootstrap_precision_model.hpp"
#include "m2424/profiles.hpp"
#include "m2424/seal_adapter.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <exception>
#include <numeric>
#include <limits>
#include <string>
#include <vector>

namespace {

constexpr std::size_t kSlots = 4;
constexpr double kAmplitude = 1e-5;
constexpr double kPlainScaleLog2 = 59.0;
constexpr double kTargetRoundtripError = 2e-10;

struct PlanCost {
    std::size_t rotation_count{};
    std::size_t plaintext_multiplication_count{};
    std::size_t addition_count{};
    std::size_t rescale_count{};
};

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

m2424::BootstrapDftPlan make_dense_plan(m2424::BootstrapDftType type, double plain_scale_log2) {
    m2424::BootstrapDftPlan plan;
    plan.slots = kSlots;
    plan.type = type;
    plan.levels = {1};
    plan.plain_scale_log2 = plain_scale_log2;
    plan.scaling_log2 = 0.0;
    const auto matrix = type == m2424::BootstrapDftType::HomomorphicDecode
        ? m2424::canonical_embedding_matrix(kSlots)
        : m2424::invert_matrix(m2424::canonical_embedding_matrix(kSlots));
    plan.layers.push_back(m2424::BootstrapDftLayer{
        type == m2424::BootstrapDftType::HomomorphicDecode ? "dense_coeff_to_slots" : "dense_slots_to_coeffs",
        m2424::DiagonalLinearTransform::from_matrix(matrix),
        plain_scale_log2,
        0.0,
        1
    });
    return plan;
}

m2424::FactorizedLinearTransform make_stc(const std::string& backend) {
    if (backend == "FftLike") {
        return m2424::FactorizedLinearTransform(m2424::make_bootstrap_dft_plan(
            kSlots, m2424::BootstrapDftType::HomomorphicEncode, kPlainScaleLog2));
    }
    if (backend == "DenseDiagonal") {
        return m2424::FactorizedLinearTransform(make_dense_plan(
            m2424::BootstrapDftType::HomomorphicEncode, kPlainScaleLog2));
    }
    return m2424::FactorizedLinearTransform(m2424::make_small_slots4_butterfly_stc_plan(kPlainScaleLog2));
}

m2424::FactorizedLinearTransform make_cts(const std::string& backend) {
    if (backend == "FftLike") {
        return m2424::FactorizedLinearTransform(m2424::make_bootstrap_dft_plan(
            kSlots, m2424::BootstrapDftType::HomomorphicDecode, kPlainScaleLog2));
    }
    if (backend == "DenseDiagonal") {
        return m2424::FactorizedLinearTransform(make_dense_plan(
            m2424::BootstrapDftType::HomomorphicDecode, kPlainScaleLog2));
    }
    return m2424::FactorizedLinearTransform(m2424::make_small_slots4_butterfly_cts_plan(kPlainScaleLog2));
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

PlanCost cost_of(const m2424::FactorizedLinearTransform& stc,
                 const m2424::FactorizedLinearTransform& cts) {
    PlanCost cost;
    cost.rotation_count = merged_steps(stc, cts).size();
    for (const auto* plan : {&stc.plan(), &cts.plan()}) {
        for (const auto& layer : plan->layers) {
            const auto terms = layer.transform.terms().size();
            cost.plaintext_multiplication_count += terms;
            cost.addition_count += terms == 0 ? 0 : terms - 1;
            cost.rescale_count += 1;
        }
    }
    return cost;
}

void print_steps(const char* backend, const char* label, const std::vector<int>& steps) {
    std::printf("[diagnose_dft_precision_model] backend=%s %s=", backend, label);
    for (std::size_t i = 0; i < steps.size(); ++i) {
        std::printf("%s%d", i == 0 ? "" : ",", steps[i]);
    }
    std::printf("\n");
}

} // namespace

int main() {
    try {
        const auto expected = make_expected();
        const auto profile = m2424::profiles::precision_boot_ultra_ckks_59();
        const auto ref_stc = m2424::FactorizedLinearTransform(m2424::make_bootstrap_dft_plan(
            kSlots, m2424::BootstrapDftType::HomomorphicEncode, kPlainScaleLog2));
        const auto ref_cts = m2424::FactorizedLinearTransform(m2424::make_bootstrap_dft_plan(
            kSlots, m2424::BootstrapDftType::HomomorphicDecode, kPlainScaleLog2));
        const auto ref_stc_expected = ref_stc.apply_plain(expected);
        const auto ref_cts_expected = ref_cts.apply_plain(ref_stc_expected);
        std::string best_backend;
        double best_roundtrip_error = std::numeric_limits<double>::infinity();
        double small_slots4_roundtrip_error = std::numeric_limits<double>::infinity();
        bool small_slots4_semantic_ok = false;

        for (const std::string backend : {"FftLike", "DenseDiagonal", "SmallSlots4Butterfly"}) {
            auto stc = make_stc(backend);
            auto cts = make_cts(backend);
            const auto rotations = merged_steps(stc, cts);
            const auto cost = cost_of(stc, cts);
            print_steps(backend.c_str(), "slot_to_coeff_rotation_steps", stc.rotation_steps());
            print_steps(backend.c_str(), "coeff_to_slot_rotation_steps", cts.rotation_steps());

            const auto plain_stc = stc.apply_plain(expected);
            const auto plain_cts = cts.apply_plain(ref_stc_expected);
            const double plain_stc_reference_error = max_error(ref_stc_expected, plain_stc);
            const double plain_cts_reference_error = max_error(ref_cts_expected, plain_cts);
            const bool semantic_ok = plain_stc_reference_error <= 1e-12 && plain_cts_reference_error <= 1e-12;

            auto adapter = m2424::SealAdapter::create(profile);
            adapter.keygen(rotations, true);
            const auto start = std::chrono::steady_clock::now();
            const auto encrypted = adapter.encrypt(adapter.encode_complex(expected));
            const auto input_info = adapter.info(encrypted);
            const double baseline_error = decoded_error(adapter, encrypted, expected);

            const auto encrypted_stc = stc.apply(adapter, encrypted);
            const auto stc_info = adapter.info(encrypted_stc);
            const double encrypted_stc_error = decoded_error(adapter, encrypted_stc, ref_stc_expected);

            const auto encrypted_cts_reference_input = adapter.encrypt(adapter.encode_complex(ref_stc_expected));
            const auto encrypted_cts_reference = cts.apply(adapter, encrypted_cts_reference_input);
            const double encrypted_cts_error = decoded_error(adapter, encrypted_cts_reference, ref_cts_expected);

            const auto encrypted_roundtrip = cts.apply(adapter, encrypted_stc);
            const auto output_info = adapter.info(encrypted_roundtrip);
            const double roundtrip_error = decoded_error(adapter, encrypted_roundtrip, expected);
            const auto finish = std::chrono::steady_clock::now();
            const double runtime_ms = std::chrono::duration<double, std::milli>(finish - start).count();

            const bool pass = semantic_ok && roundtrip_error <= kTargetRoundtripError;
            if (roundtrip_error < best_roundtrip_error) {
                best_roundtrip_error = roundtrip_error;
                best_backend = backend;
            }
            if (backend == "SmallSlots4Butterfly") {
                small_slots4_roundtrip_error = roundtrip_error;
                small_slots4_semantic_ok = semantic_ok;
            }
            std::printf("[diagnose_dft_precision_model] backend=%s profile=precision_boot_ultra_ckks_59 ciphertext_scale_log2=%.6f transform_plain_scale_log2=%.0f baseline_error=%.12e plain_stc_reference_error=%.12e plain_cts_reference_error=%.12e encrypted_stc_error=%.12e encrypted_cts_error=%.12e roundtrip_error=%.12e rotation_count=%zu plaintext_muls=%zu additions=%zu rescales=%zu input_chain=%zu stc_chain=%zu output_chain=%zu consumed_levels=%zu input_scale_log2=%.6f output_scale_log2=%.6f runtime_ms=%.3f status=%s\n",
                        backend.c_str(),
                        std::log2(input_info.scale),
                        kPlainScaleLog2,
                        baseline_error,
                        plain_stc_reference_error,
                        plain_cts_reference_error,
                        encrypted_stc_error,
                        encrypted_cts_error,
                        roundtrip_error,
                        cost.rotation_count,
                        cost.plaintext_multiplication_count,
                        cost.addition_count,
                        cost.rescale_count,
                        input_info.chain_index,
                        stc_info.chain_index,
                        output_info.chain_index,
                        input_info.chain_index >= output_info.chain_index ? input_info.chain_index - output_info.chain_index : 0,
                        std::log2(input_info.scale),
                        std::log2(output_info.scale),
                        runtime_ms,
                        pass ? "PASS" : (semantic_ok ? "BLOCKED" : "SEMANTIC_FAIL"));
        }
        const bool small_slots4_pass = small_slots4_semantic_ok && small_slots4_roundtrip_error <= kTargetRoundtripError;
        std::printf("[diagnose_dft_precision_model] summary target_roundtrip_error=%.12e best_backend=%s best_roundtrip_error=%.12e small_slots4_roundtrip_error=%.12e small_slots4_semantic_ok=%s phase5_status=%s%s\n",
                    kTargetRoundtripError,
                    best_backend.c_str(),
                    best_roundtrip_error,
                    small_slots4_roundtrip_error,
                    small_slots4_semantic_ok ? "true" : "false",
                    small_slots4_pass && best_backend == "SmallSlots4Butterfly" ? "PASS" : "BLOCKED",
                    small_slots4_pass && best_backend == "SmallSlots4Butterfly"
                        ? ""
                        : " blocker=SmallSlots4Butterfly is semantically equivalent but is not the lowest-error DFT path; stop before bootstrap integration.");
        return 0;
    } catch (const std::exception& error) {
        std::printf("[diagnose_dft_precision_model] FAIL: %s\n", error.what());
        return 1;
    }
}
