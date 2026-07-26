#include "m2424/bootstrap.hpp"
#include "m2424/bootstrap_dft.hpp"
#include "m2424/bootstrap_stc_modup.hpp"
#include "m2424/bootstrap_stc_reference.hpp"
#include "m2424/profiles.hpp"
#include "m2424/seal_adapter.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <exception>
#include <limits>

namespace {

constexpr std::size_t kSlots = 4;
constexpr double kAmplitude = 1e-5;
constexpr double kScaleLog2 = 59.0;
constexpr std::size_t kStcTargetChain = 2;
constexpr double kUsefulGain = 0.1;

const char* to_string(m2424::BootstrapModUpVariant variant) {
    switch (variant) {
    case m2424::BootstrapModUpVariant::CenteredLift:
        return "CenteredLift";
    case m2424::BootstrapModUpVariant::UncenteredLift:
        return "UncenteredLift";
    }
    return "unknown";
}

m2424::ComplexVector make_input(bool complex_mode) {
    m2424::ComplexVector input;
    input.reserve(kSlots);
    for (std::size_t i = 0; i < kSlots; ++i) {
        const double x = static_cast<double>(i + 1);
        const double real_sign = (i % 2 == 0) ? 1.0 : -1.0;
        const double imag_sign = (i % 3 == 0) ? -1.0 : 1.0;
        input.push_back({
            real_sign * kAmplitude * (0.25 + 0.1 * std::sin(x)),
            complex_mode ? imag_sign * kAmplitude * (0.15 + 0.05 * std::cos(0.5 * x)) : 0.0
        });
    }
    return input;
}

m2424::ComplexVector head(m2424::ComplexVector values) {
    values.resize(std::min(values.size(), kSlots));
    return values;
}

double max_error(const m2424::ComplexVector& lhs, const m2424::ComplexVector& rhs) {
    double result = 0.0;
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        result = std::max(result, std::abs(lhs[i] - rhs[i]));
    }
    return result;
}

m2424::ComplexVector divide_by_period_and_reduce(const m2424::ComplexVector& values, double period_log2) {
    m2424::ComplexVector scaled;
    scaled.reserve(values.size());
    const double factor = std::exp2(-period_log2);
    for (const auto& value : values) {
        scaled.push_back(value * factor);
    }
    return m2424::reduce_mod_integer_lattice(scaled);
}

m2424::Cipher drop_level_preserving_scale(m2424::SealAdapter& adapter, const m2424::Cipher& input) {
    const auto info = adapter.info(input);
    const auto bits = adapter.coeff_modulus_bits();
    if (info.chain_index == 0 || info.coeff_modulus_size == 0 || info.coeff_modulus_size > bits.size()) {
        throw std::runtime_error("cannot drop level preserving scale");
    }
    const double plain_scale_log2 = static_cast<double>(bits[info.coeff_modulus_size - 1]);
    return adapter.mul_plain_rescale(
        input,
        adapter.encode_scalar_at_scale_like(1.0, std::exp2(plain_scale_log2), input));
}

void run_case(const char* input_mode, bool complex_mode) {
    const auto input = make_input(complex_mode);
    const auto reference_plan = m2424::make_bootstrap_stc_reference_plan(kSlots, kScaleLog2, 1e-9);
    const auto reference = m2424::evaluate_bootstrap_stc_reference(reference_plan, input);

    auto adapter = m2424::SealAdapter::create(m2424::profiles::precision_boot_ultra_ckks_59());
    adapter.keygen(m2424::Bootstrapper::scalable_refresh_rotation_steps(kSlots), true);
    const auto stc = m2424::FactorizedLinearTransform(m2424::make_bootstrap_dft_plan(
        kSlots, m2424::BootstrapDftType::HomomorphicEncode, kScaleLog2));
    const auto cts = m2424::FactorizedLinearTransform(m2424::make_bootstrap_dft_plan(
        kSlots, m2424::BootstrapDftType::HomomorphicDecode, kScaleLog2));

    auto current = adapter.encrypt(adapter.encode_complex(input));
    while (adapter.info(current).chain_index > kStcTargetChain) {
        current = drop_level_preserving_scale(adapter, current);
    }
    current = stc.apply(adapter, current);
    const auto after_stc = adapter.info(current);
    for (const auto target_q_size : {std::size_t{1}, std::size_t{2}}) {
    for (const auto variant : {m2424::BootstrapModUpVariant::CenteredLift,
                               m2424::BootstrapModUpVariant::UncenteredLift}) {
        m2424::BootstrapScaleDownToQPlan scale_down_plan;
        scale_down_plan.message_scale_log2 = reference_plan.message_scale_log2;
        scale_down_plan.target_scale_log2 = reference_plan.message_scale_log2;
        scale_down_plan.target_coeff_modulus_size = target_q_size;
        scale_down_plan.preserve_scale_on_level_drop = true;
        const auto scaled = m2424::bootstrap_scale_down_to_q(adapter, current, scale_down_plan);
        auto modup = adapter.bootstrap_modup_to_first(scaled.result, variant);
        const auto after_modup = adapter.info(modup);
        auto after_cts_cipher = cts.apply(adapter, modup);
        const auto after_cts = adapter.info(after_cts_cipher);
        const auto encrypted_z = head(adapter.decode_complex(adapter.decrypt(after_cts_cipher)));

        const auto encrypted_scan = m2424::scan_bootstrap_lattice_periods(encrypted_z, input, 80.0, 115.0, kUsefulGain);
        const auto encrypted_fractional = divide_by_period_and_reduce(encrypted_z, reference_plan.modup_period_log2);
        const auto reference_fractional = divide_by_period_and_reduce(reference.after_coeff_to_slots,
                                                                      reference_plan.modup_period_log2);

        const double err_absolute = max_error(encrypted_z, reference.after_coeff_to_slots);
        const double err_fractional = max_error(encrypted_fractional, reference_fractional);

        const char* classification = "missing ScaleDown/ModUp semantics";
        if (!complex_mode && reference.lattice_invariant_passed && err_fractional <= 1e-9) {
            classification = "matches reference";
        } else if (complex_mode && !reference.lattice_invariant_passed) {
            classification = "unsupported complex handling";
        } else if (encrypted_scan.has_useful_gain
                   && encrypted_scan.best_useful_mod_gain.err_mod_gain <= 1e-9
                   && std::abs(encrypted_scan.best_useful_mod_gain.gamma - reference.best_gamma) > 1e-3) {
            classification = "scalar gain error only";
        } else if (encrypted_scan.best_mod_direct.err_mod_direct <= 1e-9) {
            classification = "period error only";
        } else if (err_fractional > 1e-9) {
            classification = "transform semantic mismatch";
        }

        std::printf("[diagnose_bootstrap_encrypted_vs_reference] input_mode=%s target_coeff_modulus_size=%zu modup_variant=%s scale_down_to_q=true reference_pass=%s reference_best_p=%.0f reference_gamma_abs=%.12e reference_gamma_arg=%.12e reference_err_mod_gain=%.12e encrypted_has_useful=%s encrypted_best_p=%.0f encrypted_gamma_abs=%.12e encrypted_gamma_arg=%.12e encrypted_err_mod_direct=%.12e encrypted_err_mod_gain=%.12e encrypted_vs_reference_abs_error=%.12e encrypted_vs_reference_fractional_error=%.12e after_stc_chain=%zu after_stc_scale_log2=%.6f after_scale_down_chain=%zu after_scale_down_coeff_modulus_size=%zu after_scale_down_scale_log2=%.6f after_modup_chain=%zu after_modup_coeff_modulus_size=%zu after_modup_scale_log2=%.6f after_cts_chain=%zu after_cts_scale_log2=%.6f classification=\"%s\"\n",
                    input_mode,
                    target_q_size,
                    to_string(variant),
                    reference.lattice_invariant_passed ? "true" : "false",
                    reference.best_period_log2,
                    std::abs(reference.best_gamma),
                    std::arg(reference.best_gamma),
                    reference.err_mod_gain,
                    encrypted_scan.has_useful_gain ? "true" : "false",
                    encrypted_scan.has_useful_gain ? encrypted_scan.best_useful_mod_gain.period_log2 : std::nan(""),
                    encrypted_scan.has_useful_gain ? std::abs(encrypted_scan.best_useful_mod_gain.gamma) : std::nan(""),
                    encrypted_scan.has_useful_gain ? std::arg(encrypted_scan.best_useful_mod_gain.gamma) : std::nan(""),
                    encrypted_scan.best_mod_direct.err_mod_direct,
                    encrypted_scan.has_useful_gain ? encrypted_scan.best_useful_mod_gain.err_mod_gain : std::numeric_limits<double>::infinity(),
                    err_absolute,
                    err_fractional,
                    after_stc.chain_index,
                    std::log2(after_stc.scale),
                    scaled.chain_after,
                    scaled.coeff_modulus_size_after,
                    scaled.scale_after_log2,
                    after_modup.chain_index,
                    after_modup.coeff_modulus_size,
                    std::log2(after_modup.scale),
                    after_cts.chain_index,
                    std::log2(after_cts.scale),
                    classification);
    }
    }
}

} // namespace

int main() {
    try {
        run_case("real_only", false);
        run_case("complex", true);
        return 0;
    } catch (const std::exception& error) {
        std::printf("[diagnose_bootstrap_encrypted_vs_reference] FAIL: %s\n", error.what());
        return 1;
    }
}
