#include "m2424/bootstrap.hpp"
#include "m2424/bootstrap_dft.hpp"
#include "m2424/bootstrap_stc_modup.hpp"
#include "m2424/profiles.hpp"
#include "m2424/seal_adapter.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <exception>
#include <limits>
#include <string>
#include <vector>

namespace {

constexpr std::size_t kSlots = 4;
constexpr double kAmplitude = 1e-5;
constexpr double kMessageScaleLog2 = 59.0;
constexpr std::size_t kStcTargetChain = 2;
constexpr double kMinimumUsefulGainAbs = 0.1;

const char* to_string(m2424::BootstrapModUpVariant variant) {
    switch (variant) {
    case m2424::BootstrapModUpVariant::CenteredLift:
        return "CenteredLift";
    case m2424::BootstrapModUpVariant::UncenteredLift:
        return "UncenteredLift";
    }
    return "unknown";
}

struct PeriodScanResult {
    double period_log2{};
    std::complex<double> gamma{};
    double err_mod_direct{};
    double err_mod_gain{};
    double max_abs_u{};
    double max_abs_r{};
};

struct BestScan {
    PeriodScanResult mod_direct;
    PeriodScanResult mod_gain;
    PeriodScanResult useful_mod_gain;
};

m2424::ComplexVector make_expected(bool complex_mode) {
    m2424::ComplexVector expected;
    expected.reserve(kSlots);
    for (std::size_t i = 0; i < kSlots; ++i) {
        const double x = static_cast<double>(i + 1);
        const double real_sign = (i % 2 == 0) ? 1.0 : -1.0;
        const double imag_sign = (i % 3 == 0) ? -1.0 : 1.0;
        const double real = real_sign * kAmplitude * (0.25 + 0.1 * std::sin(x));
        const double imag = complex_mode ? imag_sign * kAmplitude * (0.15 + 0.05 * std::cos(0.5 * x)) : 0.0;
        expected.push_back({real, imag});
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

double max_error(const m2424::ComplexVector& lhs, const m2424::ComplexVector& rhs) {
    double result = 0.0;
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        result = std::max(result, std::abs(lhs[i] - rhs[i]));
    }
    return result;
}

std::complex<double> best_gain(const m2424::ComplexVector& actual,
                               const m2424::ComplexVector& expected) {
    std::complex<double> numerator{0.0, 0.0};
    double denominator = 0.0;
    for (std::size_t i = 0; i < actual.size(); ++i) {
        numerator += std::conj(expected[i]) * actual[i];
        denominator += std::norm(expected[i]);
    }
    return denominator == 0.0 ? std::complex<double>{0.0, 0.0} : numerator / denominator;
}

m2424::ComplexVector scale_expected(const m2424::ComplexVector& expected,
                                    std::complex<double> gamma) {
    m2424::ComplexVector result;
    result.reserve(expected.size());
    for (const auto& value : expected) {
        result.push_back(gamma * value);
    }
    return result;
}

m2424::ComplexVector divide_by_period(const m2424::ComplexVector& values, double period_log2) {
    const double factor = std::exp2(-period_log2);
    m2424::ComplexVector result;
    result.reserve(values.size());
    for (const auto& value : values) {
        result.push_back(value * factor);
    }
    return result;
}

m2424::ComplexVector reduce_mod_integer_lattice(const m2424::ComplexVector& values) {
    m2424::ComplexVector result;
    result.reserve(values.size());
    for (const auto& value : values) {
        result.push_back({
            value.real() - std::round(value.real()),
            value.imag() - std::round(value.imag())
        });
    }
    return result;
}

PeriodScanResult scan_one(double period_log2,
                          const m2424::ComplexVector& z,
                          const m2424::ComplexVector& expected) {
    const auto u = divide_by_period(z, period_log2);
    const auto r = reduce_mod_integer_lattice(u);
    const auto gamma_mod = best_gain(r, expected);

    PeriodScanResult result;
    result.period_log2 = period_log2;
    result.gamma = gamma_mod;
    result.err_mod_direct = max_error(r, expected);
    result.err_mod_gain = max_error(r, scale_expected(expected, gamma_mod));
    result.max_abs_u = max_abs_value(u);
    result.max_abs_r = max_abs_value(r);
    return result;
}

BestScan scan_periods(const m2424::ComplexVector& z, const m2424::ComplexVector& expected) {
    BestScan best;
    best.mod_direct.err_mod_direct = std::numeric_limits<double>::infinity();
    best.mod_gain.err_mod_gain = std::numeric_limits<double>::infinity();
    best.useful_mod_gain.err_mod_gain = std::numeric_limits<double>::infinity();
    for (int p = 80; p <= 115; ++p) {
        const auto result = scan_one(static_cast<double>(p), z, expected);
        if (result.err_mod_direct < best.mod_direct.err_mod_direct) {
            best.mod_direct = result;
        }
        if (result.err_mod_gain < best.mod_gain.err_mod_gain) {
            best.mod_gain = result;
        }
        if (std::abs(result.gamma) >= kMinimumUsefulGainAbs
            && result.err_mod_gain < best.useful_mod_gain.err_mod_gain) {
            best.useful_mod_gain = result;
        }
    }
    return best;
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

void run_mode(const char* input_mode, bool complex_mode) {
    const auto expected = make_expected(complex_mode);
    auto adapter = m2424::SealAdapter::create(m2424::profiles::precision_boot_ultra_ckks_59());
    adapter.keygen(m2424::Bootstrapper::scalable_refresh_rotation_steps(kSlots), true);

    const auto stc = m2424::FactorizedLinearTransform(m2424::make_bootstrap_dft_plan(
        kSlots, m2424::BootstrapDftType::HomomorphicEncode, kMessageScaleLog2));
    const auto cts = m2424::FactorizedLinearTransform(m2424::make_bootstrap_dft_plan(
        kSlots, m2424::BootstrapDftType::HomomorphicDecode, kMessageScaleLog2));

    auto current = adapter.encrypt(adapter.encode_complex(expected));
    const auto input_info = adapter.info(current);
    while (adapter.info(current).chain_index > kStcTargetChain) {
        current = drop_level_preserving_scale(adapter, current);
    }
    const auto after_burn = adapter.info(current);
    current = stc.apply(adapter, current);
    const auto after_stc = adapter.info(current);

    std::printf("[diagnose_stc_scaled_modup_lattice] input_mode=%s stage=input chain=%zu scale_log2=%.6f coeff_log2=%.0f\n",
                input_mode, input_info.chain_index, std::log2(input_info.scale), input_info.coeff_modulus_log2);
    std::printf("[diagnose_stc_scaled_modup_lattice] input_mode=%s stage=after_burn chain=%zu scale_log2=%.6f coeff_log2=%.0f\n",
                input_mode, after_burn.chain_index, std::log2(after_burn.scale), after_burn.coeff_modulus_log2);
    std::printf("[diagnose_stc_scaled_modup_lattice] input_mode=%s stage=after_slots_to_coeff chain=%zu scale_log2=%.6f coeff_log2=%.0f\n",
                input_mode, after_stc.chain_index, std::log2(after_stc.scale), after_stc.coeff_modulus_log2);
    for (const auto target_q_size : {std::size_t{1}, std::size_t{2}}) {
    for (const auto variant : {m2424::BootstrapModUpVariant::CenteredLift,
                               m2424::BootstrapModUpVariant::UncenteredLift}) {
        m2424::BootstrapScaleDownToQPlan scale_down_plan;
        scale_down_plan.message_scale_log2 = kMessageScaleLog2;
        scale_down_plan.target_scale_log2 = kMessageScaleLog2;
        scale_down_plan.target_coeff_modulus_size = target_q_size;
        scale_down_plan.preserve_scale_on_level_drop = true;
        auto scaled = m2424::bootstrap_scale_down_to_q(adapter, current, scale_down_plan);
        auto modup = adapter.bootstrap_modup_to_first(scaled.result, variant);
        const auto after_modup = adapter.info(modup);
        auto after_cts_cipher = cts.apply(adapter, modup);
        const auto after_cts = adapter.info(after_cts_cipher);
        const auto z = head(adapter.decode_complex(adapter.decrypt(after_cts_cipher)));
        const auto best = scan_periods(z, expected);

        std::printf("[diagnose_stc_scaled_modup_lattice] input_mode=%s target_coeff_modulus_size=%zu modup_variant=%s scale_down_to_q=true after_scale_down_chain=%zu after_scale_down_coeff_modulus_size=%zu after_scale_down_scale_log2=%.6f after_scale_down_coeff_log2=%.0f after_modup_chain=%zu after_modup_coeff_modulus_size=%zu after_modup_scale_log2=%.6f after_modup_coeff_log2=%.0f after_coeff_to_slot_chain=%zu after_coeff_to_slot_coeff_modulus_size=%zu after_coeff_to_slot_scale_log2=%.6f after_coeff_to_slot_coeff_log2=%.0f max_abs_z=%.12e\n",
                    input_mode,
                    target_q_size,
                    to_string(variant),
                    scaled.chain_after,
                    scaled.coeff_modulus_size_after,
                    scaled.scale_after_log2,
                    scaled.coeff_modulus_log2_after,
                    after_modup.chain_index,
                    after_modup.coeff_modulus_size,
                    std::log2(after_modup.scale),
                    after_modup.coeff_modulus_log2,
                    after_cts.chain_index,
                    after_cts.coeff_modulus_size,
                    std::log2(after_cts.scale),
                    after_cts.coeff_modulus_log2,
                    max_abs_value(z));

        const bool useful = std::isfinite(best.useful_mod_gain.err_mod_gain);
        const bool pass = useful && best.useful_mod_gain.err_mod_gain <= 1e-9;
        std::printf("[diagnose_stc_scaled_modup_lattice] input_mode=%s target_coeff_modulus_size=%zu modup_variant=%s best_period_log2=%.0f best_gamma_abs=%.12e best_gamma_arg=%.12e err_mod_direct=%.12e err_mod_gain=%.12e useful_best_period_log2=%.0f useful_best_gamma_abs=%.12e useful_best_gamma_arg=%.12e useful_err_mod_gain=%.12e useful_status=%s classification=%s conclusion=%s\n",
                    input_mode,
                    target_q_size,
                    to_string(variant),
                    best.mod_gain.period_log2,
                    std::abs(best.mod_gain.gamma),
                    std::arg(best.mod_gain.gamma),
                    best.mod_direct.err_mod_direct,
                    best.mod_gain.err_mod_gain,
                    useful ? best.useful_mod_gain.period_log2 : std::nan(""),
                    useful ? std::abs(best.useful_mod_gain.gamma) : std::nan(""),
                    useful ? std::arg(best.useful_mod_gain.gamma) : std::nan(""),
                    useful ? best.useful_mod_gain.err_mod_gain : std::numeric_limits<double>::infinity(),
                    useful ? "non_degenerate" : "degenerate",
                    pass ? "matches reference lattice" : "transform semantic mismatch",
                    pass ? "lattice_invariant_pass" : "lattice_invariant_blocked");
    }
    }
}

} // namespace

int main() {
    try {
        run_mode("real_only", false);
        run_mode("complex", true);
        return 0;
    } catch (const std::exception& error) {
        std::printf("[diagnose_stc_scaled_modup_lattice] FAIL: %s\n", error.what());
        return 1;
    }
}
