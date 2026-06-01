#include "m2424/bootstrap.hpp"
#include "m2424/bootstrap_dft.hpp"
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
constexpr double kPlainScaleLog2 = 59.0;
constexpr std::size_t kStcTargetChain = 2;
constexpr double kMinimumUsefulGainAbs = 0.1;

struct ScaleDownResult {
    m2424::Cipher result;
    std::size_t chain_before{};
    std::size_t chain_after{};
    double scale_before_log2{};
    double scale_after_log2{};
    std::size_t levels_consumed{};
};

struct PeriodScanResult {
    double period_log2{};
    std::complex<double> gamma{};
    double err_direct{};
    double err_gain{};
    double err_mod_direct{};
    double err_mod_gain{};
    double max_abs_u{};
    double max_abs_r{};
};

struct BestSummary {
    PeriodScanResult direct;
    PeriodScanResult gain;
    PeriodScanResult mod_direct;
    PeriodScanResult mod_gain;
    PeriodScanResult useful_mod_gain;
};

struct PipelineResult {
    m2424::ComplexVector z;
    std::size_t chain_after_scale_down{};
    std::size_t chain_after_cts{};
    double scale_after_scale_down_log2{};
    double scale_after_cts_log2{};
    std::size_t levels_consumed{};
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
    if (denominator == 0.0) {
        return {0.0, 0.0};
    }
    return numerator / denominator;
}

m2424::ComplexVector scale(const m2424::ComplexVector& values, std::complex<double> factor) {
    m2424::ComplexVector result;
    result.reserve(values.size());
    for (const auto& value : values) {
        result.push_back(factor * value);
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

PeriodScanResult scan_period(const char* input_mode,
                             double scale_down_log2,
                             double period_log2,
                             const m2424::ComplexVector& z,
                             const m2424::ComplexVector& expected,
                             bool print_each) {
    const auto u = divide_by_period(z, period_log2);
    const auto r = reduce_mod_integer_lattice(u);
    const auto gamma = best_gain(u, expected);
    const auto gamma_mod = best_gain(r, expected);

    PeriodScanResult result;
    result.period_log2 = period_log2;
    result.gamma = gamma_mod;
    result.err_direct = max_error(u, expected);
    result.err_gain = max_error(u, scale(expected, gamma));
    result.err_mod_direct = max_error(r, expected);
    result.err_mod_gain = max_error(r, scale(expected, gamma_mod));
    result.max_abs_u = max_abs_value(u);
    result.max_abs_r = max_abs_value(r);

    if (print_each) {
        std::printf("[diagnose_bootstrap_lattice_invariant] input_mode=%s scale_down_log2=%.0f p=%.0f abs_gamma=%.12e arg_gamma=%.12e err_direct=%.12e err_gain=%.12e err_mod_direct=%.12e err_mod_gain=%.12e max_abs_u=%.12e max_abs_r=%.12e\n",
                    input_mode,
                    scale_down_log2,
                    period_log2,
                    std::abs(gamma_mod),
                    std::arg(gamma_mod),
                    result.err_direct,
                    result.err_gain,
                    result.err_mod_direct,
                    result.err_mod_gain,
                    result.max_abs_u,
                    result.max_abs_r);
    }
    return result;
}

BestSummary scan_periods(const char* input_mode,
                         double scale_down_log2,
                         const m2424::ComplexVector& z,
                         const m2424::ComplexVector& expected,
                         bool print_each) {
    BestSummary best;
    best.direct.err_direct = std::numeric_limits<double>::infinity();
    best.gain.err_gain = std::numeric_limits<double>::infinity();
    best.mod_direct.err_mod_direct = std::numeric_limits<double>::infinity();
    best.mod_gain.err_mod_gain = std::numeric_limits<double>::infinity();
    best.useful_mod_gain.err_mod_gain = std::numeric_limits<double>::infinity();
    for (int p = 80; p <= 115; ++p) {
        const auto result = scan_period(
            input_mode, scale_down_log2, static_cast<double>(p), z, expected, print_each);
        if (result.err_direct < best.direct.err_direct) {
            best.direct = result;
        }
        if (result.err_gain < best.gain.err_gain) {
            best.gain = result;
        }
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

void print_best(const char* input_mode,
                double scale_down_log2,
                const char* metric,
                const PeriodScanResult& best) {
    std::printf("[diagnose_bootstrap_lattice_invariant] input_mode=%s scale_down_log2=%.0f best_metric=%s best_p=%.0f best_gamma_abs=%.12e best_gamma_arg=%.12e err_direct=%.12e err_gain=%.12e err_mod_direct=%.12e err_mod_gain=%.12e max_abs_u=%.12e max_abs_r=%.12e\n",
                input_mode,
                scale_down_log2,
                metric,
                best.period_log2,
                std::abs(best.gamma),
                std::arg(best.gamma),
                best.err_direct,
                best.err_gain,
                best.err_mod_direct,
                best.err_mod_gain,
                best.max_abs_u,
                best.max_abs_r);
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

ScaleDownResult apply_bootstrap_scale_down_candidate(m2424::SealAdapter& adapter,
                                                     const m2424::Cipher& input,
                                                     double scale_down_log2) {
    const auto before = adapter.info(input);
    ScaleDownResult result;
    result.chain_before = before.chain_index;
    result.scale_before_log2 = std::log2(before.scale);
    if (scale_down_log2 == 0.0) {
        result.result = input;
        result.chain_after = before.chain_index;
        result.scale_after_log2 = result.scale_before_log2;
        result.levels_consumed = 0;
        return result;
    }
    if (!std::isfinite(scale_down_log2) || scale_down_log2 < 0.0) {
        throw std::invalid_argument("scale_down_log2 must be non-negative and finite");
    }
    if (before.chain_index == 0) {
        throw std::runtime_error("cannot apply scale down at chain index 0");
    }
    const auto bits = adapter.coeff_modulus_bits();
    if (before.coeff_modulus_size == 0 || before.coeff_modulus_size > bits.size()) {
        throw std::runtime_error("cannot infer scale down rescale prime");
    }
    const double plain_scale_log2 = static_cast<double>(bits[before.coeff_modulus_size - 1]);
    result.result = adapter.mul_plain_rescale(
        input,
        adapter.encode_scalar_at_scale_like(std::exp2(-scale_down_log2),
                                            std::exp2(plain_scale_log2),
                                            input));
    const auto after = adapter.info(result.result);
    if (!std::isfinite(after.scale) || after.scale <= 0.0) {
        throw std::runtime_error("scale down produced invalid ciphertext scale");
    }
    result.chain_after = after.chain_index;
    result.scale_after_log2 = std::log2(after.scale);
    result.levels_consumed = before.chain_index - after.chain_index;
    return result;
}

PipelineResult run_pre_evalmod(m2424::SealAdapter& adapter,
                               const m2424::ComplexVector& expected,
                               double scale_down_log2) {
    const auto stc = m2424::FactorizedLinearTransform(m2424::make_bootstrap_dft_plan(
        kSlots, m2424::BootstrapDftType::HomomorphicEncode, kPlainScaleLog2));
    const auto cts = m2424::FactorizedLinearTransform(m2424::make_bootstrap_dft_plan(
        kSlots, m2424::BootstrapDftType::HomomorphicDecode, kPlainScaleLog2));

    auto current = adapter.encrypt(adapter.encode_complex(expected));
    const auto start_chain = adapter.info(current).chain_index;
    while (adapter.info(current).chain_index > kStcTargetChain) {
        current = drop_level_preserving_scale(adapter, current);
    }
    current = stc.apply(adapter, current);
    const auto after_stc = adapter.info(current);
    const auto scaled = apply_bootstrap_scale_down_candidate(adapter, current, scale_down_log2);
    current = adapter.mod_raise_to_first(scaled.result);
    current = cts.apply(adapter, current);
    const auto after_cts = adapter.info(current);

    PipelineResult result;
    result.z = head(adapter.decode_complex(adapter.decrypt(current)));
    result.chain_after_scale_down = scaled.chain_after;
    result.chain_after_cts = after_cts.chain_index;
    result.scale_after_scale_down_log2 = scaled.scale_after_log2;
    result.scale_after_cts_log2 = std::log2(after_cts.scale);
    result.levels_consumed = start_chain - after_cts.chain_index;
    std::printf("[diagnose_bootstrap_lattice_invariant] input_mode=pipeline scale_down_candidate scale_down_log2=%.0f stc_chain=%zu stc_scale_log2=%.6f scale_down_chain_before=%zu scale_down_chain_after=%zu scale_down_scale_before_log2=%.6f scale_down_scale_after_log2=%.6f cts_chain=%zu cts_scale_log2=%.6f levels_consumed=%zu\n",
                scale_down_log2,
                after_stc.chain_index,
                std::log2(after_stc.scale),
                scaled.chain_before,
                scaled.chain_after,
                scaled.scale_before_log2,
                scaled.scale_after_log2,
                after_cts.chain_index,
                result.scale_after_cts_log2,
                result.levels_consumed);
    return result;
}

void run_mode(const char* input_mode, bool complex_mode) {
    const auto expected = make_expected(complex_mode);
    auto adapter = m2424::SealAdapter::create(m2424::profiles::precision_boot_ultra_ckks_59());
    adapter.keygen(m2424::Bootstrapper::scalable_refresh_rotation_steps(kSlots), true);

    const std::vector<double> scale_down_candidates{0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 59.0};
    double best_scale_down = 0.0;
    double best_useful_scale_down = 0.0;
    BestSummary best_overall;
    best_overall.mod_gain.err_mod_gain = std::numeric_limits<double>::infinity();
    PeriodScanResult best_useful;
    best_useful.err_mod_gain = std::numeric_limits<double>::infinity();

    for (double scale_down_log2 : scale_down_candidates) {
        const auto pipeline = run_pre_evalmod(adapter, expected, scale_down_log2);
        const auto best = scan_periods(input_mode, scale_down_log2, pipeline.z, expected, scale_down_log2 == 0.0);
        print_best(input_mode, scale_down_log2, "err_direct", best.direct);
        print_best(input_mode, scale_down_log2, "err_gain", best.gain);
        print_best(input_mode, scale_down_log2, "err_mod_direct", best.mod_direct);
        print_best(input_mode, scale_down_log2, "err_mod_gain", best.mod_gain);
        if (std::isfinite(best.useful_mod_gain.err_mod_gain)) {
            print_best(input_mode, scale_down_log2, "useful_err_mod_gain", best.useful_mod_gain);
        }
        std::printf("[diagnose_bootstrap_lattice_invariant] input_mode=%s scale_down_log2=%.0f summary chain_after_scale_down=%zu scale_after_scale_down_log2=%.6f chain_after_cts=%zu scale_after_cts_log2=%.6f levels_consumed=%zu best_p=%.0f best_err_mod_direct=%.12e best_err_mod_gain=%.12e best_gamma_mod_abs=%.12e best_gamma_mod_arg=%.12e",
                    input_mode,
                    scale_down_log2,
                    pipeline.chain_after_scale_down,
                    pipeline.scale_after_scale_down_log2,
                    pipeline.chain_after_cts,
                    pipeline.scale_after_cts_log2,
                    pipeline.levels_consumed,
                    best.mod_gain.period_log2,
                    best.mod_direct.err_mod_direct,
                    best.mod_gain.err_mod_gain,
                    std::abs(best.mod_gain.gamma),
                    std::arg(best.mod_gain.gamma));
        if (std::isfinite(best.useful_mod_gain.err_mod_gain)) {
            std::printf(" useful_best_p=%.0f useful_best_err_mod_gain=%.12e useful_best_gamma_abs=%.12e useful_best_gamma_arg=%.12e",
                        best.useful_mod_gain.period_log2,
                        best.useful_mod_gain.err_mod_gain,
                        std::abs(best.useful_mod_gain.gamma),
                        std::arg(best.useful_mod_gain.gamma));
        } else {
            std::printf(" useful_best_p=nan useful_best_err_mod_gain=inf useful_best_gamma_abs=nan useful_best_gamma_arg=nan");
        }
        std::printf("\n");
        if (best.mod_gain.err_mod_gain < best_overall.mod_gain.err_mod_gain) {
            best_overall = best;
            best_scale_down = scale_down_log2;
        }
        if (best.useful_mod_gain.err_mod_gain < best_useful.err_mod_gain) {
            best_useful = best.useful_mod_gain;
            best_useful_scale_down = scale_down_log2;
        }
    }

    const char* conclusion = "lattice_invariant_fails";
    if (best_overall.mod_direct.err_mod_direct <= 1e-9) {
        conclusion = "period_only";
    } else if (std::isfinite(best_useful.err_mod_gain) && best_useful.err_mod_gain <= 1e-9) {
        conclusion = "period_plus_gain";
    } else if (best_overall.mod_gain.err_mod_gain <= 1e-9) {
        conclusion = "degenerate_zero_gain_fit";
    }
    std::printf("[diagnose_bootstrap_lattice_invariant] input_mode=%s overall_best_scale_down_log2=%.0f overall_best_p=%.0f overall_best_gamma_abs=%.12e overall_best_gamma_arg=%.12e overall_best_err_mod_direct=%.12e overall_best_err_mod_gain=%.12e",
                input_mode,
                best_scale_down,
                best_overall.mod_gain.period_log2,
                std::abs(best_overall.mod_gain.gamma),
                std::arg(best_overall.mod_gain.gamma),
                best_overall.mod_direct.err_mod_direct,
                best_overall.mod_gain.err_mod_gain);
    if (std::isfinite(best_useful.err_mod_gain)) {
        std::printf(" overall_useful_best_scale_down_log2=%.0f overall_useful_best_p=%.0f overall_useful_best_gamma_abs=%.12e overall_useful_best_gamma_arg=%.12e overall_useful_best_err_mod_gain=%.12e",
                    best_useful_scale_down,
                    best_useful.period_log2,
                    std::abs(best_useful.gamma),
                    std::arg(best_useful.gamma),
                    best_useful.err_mod_gain);
    } else {
        std::printf(" overall_useful_best_scale_down_log2=nan overall_useful_best_p=nan overall_useful_best_gamma_abs=nan overall_useful_best_gamma_arg=nan overall_useful_best_err_mod_gain=inf");
    }
    std::printf(" conclusion=%s\n", conclusion);
}

} // namespace

int main() {
    try {
        run_mode("real_only", false);
        run_mode("complex", true);
        return 0;
    } catch (const std::exception& error) {
        std::printf("[diagnose_bootstrap_lattice_invariant] FAIL: %s\n", error.what());
        return 1;
    }
}
