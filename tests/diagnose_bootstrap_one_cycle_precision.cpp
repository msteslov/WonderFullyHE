#include "m2424/bootstrap.hpp"
#include "m2424/profiles.hpp"
#include "m2424/seal_adapter.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>
#include <string>
#include <vector>

namespace {

constexpr std::size_t kSlots = 4;
constexpr double kAmplitude = 1e-5;
constexpr double kTolerance = 1e-9;

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

double mean_error(const m2424::ComplexVector& expected, const m2424::ComplexVector& actual) {
    double result = 0.0;
    for (std::size_t i = 0; i < kSlots; ++i) {
        result += std::abs(actual[i] - expected[i]);
    }
    return result / static_cast<double>(kSlots);
}

double decoded_max_error(m2424::SealAdapter& adapter,
                         const m2424::Cipher& cipher,
                         const m2424::ComplexVector& expected) {
    return max_error(expected, head(adapter.decode_complex(adapter.decrypt(cipher))));
}

const m2424::BootstrapPrototypeStage* find_stage(const m2424::BootstrapPrototypeReport& report,
                                                 const char* name) {
    for (const auto& stage : report.stages) {
        if (stage.name == name) {
            return &stage;
        }
    }
    return nullptr;
}

double stage_error_or_nan(const m2424::BootstrapPrototypeReport& report, const char* name) {
    const auto* stage = find_stage(report, name);
    return stage == nullptr ? std::nan("") : stage->max_abs_error;
}

void print_stage(const m2424::BootstrapPrototypeStage& stage) {
    std::printf("[diagnose_bootstrap_one_cycle_precision] stage=%s status=%s chain_before=%zu chain_after=%zu scale_before_log2=%.6f scale_after_log2=%.6f max_error=%.12e duration_ms=%.3f\n",
                stage.name.c_str(),
                stage.status.c_str(),
                stage.chain_before,
                stage.chain_after,
                std::log2(stage.scale_before),
                std::log2(stage.scale_after),
                stage.max_abs_error,
                stage.duration_ms);
}

bool has_failed_stage(const m2424::BootstrapPrototypeReport& report, std::string& first_failure) {
    for (const auto& stage : report.stages) {
        if (stage.status == "FAIL") {
            first_failure = stage.name;
            return true;
        }
    }
    return false;
}

std::string classify_exception_stage(const std::string& message) {
    const std::vector<std::pair<std::string, std::string>> stages{
        {"stc_first_level_drop", "stc_first_level_drop"},
        {"slot_to_coeff_first", "slot_to_coeff_first"},
        {"mod_raise", "mod_raise"},
        {"coeff_to_slot_after_raise", "coeff_to_slot_after_raise"},
        {"eval_mod_normalization", "eval_mod_normalization"},
        {"eval_mod", "eval_mod"},
        {"output_scale_repair", "output_scale_repair"},
        {"refresh_result", "refresh_result"}
    };
    for (const auto& [needle, stage] : stages) {
        if (message.find(needle) != std::string::npos) {
            return stage;
        }
    }
    return "unknown";
}

bool run_cycle(std::size_t cycle,
               m2424::SealAdapter& adapter,
               const m2424::ComplexVector& expected,
               m2424::Cipher& current) {
    std::printf("[diagnose_bootstrap_one_cycle_precision] cycle=%zu step=baseline\n", cycle);
    std::fflush(stdout);
    const auto before_info = adapter.info(current);
    const double baseline_error = decoded_max_error(adapter, current, expected);

    m2424::Bootstrapper bootstrapper(adapter);
    std::printf("[diagnose_bootstrap_one_cycle_precision] cycle=%zu step=refresh\n", cycle);
    std::fflush(stdout);
    m2424::BootstrapPrototypeReport report;
    try {
        report = bootstrapper.refresh_slots_to_coeffs_first_checked(
            current, expected, kSlots, kTolerance, m2424::EvalModDegree::P3DoubleAngle);
    } catch (const std::exception& error) {
        const std::string reason = error.what();
        std::printf("[diagnose_bootstrap_one_cycle_precision] cycle=%zu status=BLOCKED first_failing_stage=%s reason=%s\n",
                    cycle,
                    classify_exception_stage(reason).c_str(),
                    reason.c_str());
        return false;
    }

    std::printf("[diagnose_bootstrap_one_cycle_precision] cycle=%zu step=final_decrypt\n", cycle);
    std::fflush(stdout);
    const auto actual = head(adapter.decode_complex(adapter.decrypt(report.result)));
    const double final_max_error = max_error(expected, actual);
    const double final_mean_error = mean_error(expected, actual);
    const auto after_info = adapter.info(report.result);

    std::printf("[diagnose_bootstrap_one_cycle_precision] cycle=%zu baseline_error=%.12e slot_to_coeff_first_error=%.12e mod_raise_drift=%.12e coeff_to_slot_after_raise_error=%.12e eval_mod_normalization_error=%.12e eval_mod_error=%.12e output_scale_repair_error=%.12e refresh_result_error=%.12e final_max_error=%.12e final_mean_error=%.12e input_chain=%zu output_chain=%zu input_scale_log2=%.6f output_scale_log2=%.6f continuation_levels=%zu inside_evalmod_interval=%s preserve_value_criterion=%s restore_level_criterion=%s stages=%zu\n",
                cycle,
                baseline_error,
                stage_error_or_nan(report, "slot_to_coeff_first"),
                stage_error_or_nan(report, "mod_raise"),
                stage_error_or_nan(report, "coeff_to_slot_after_raise"),
                stage_error_or_nan(report, "eval_mod_normalization"),
                stage_error_or_nan(report, "eval_mod"),
                stage_error_or_nan(report, "output_scale_repair"),
                stage_error_or_nan(report, "refresh_result"),
                final_max_error,
                final_mean_error,
                before_info.chain_index,
                after_info.chain_index,
                std::log2(before_info.scale),
                std::log2(after_info.scale),
                report.continuation_levels,
                report.inside_evalmod_interval ? "true" : "false",
                report.preserve_value_criterion ? "true" : "false",
                report.restore_level_criterion ? "true" : "false",
                report.stages.size());

    for (const auto& stage : report.stages) {
        print_stage(stage);
    }

    std::string first_failure;
    if (has_failed_stage(report, first_failure)) {
        std::printf("[diagnose_bootstrap_one_cycle_precision] cycle=%zu status=BLOCKED first_failing_stage=%s reason=stage_status_FAIL\n",
                    cycle,
                    first_failure.c_str());
        return false;
    }
    if (!report.inside_evalmod_interval) {
        std::printf("[diagnose_bootstrap_one_cycle_precision] cycle=%zu status=BLOCKED first_failing_stage=eval_mod_interval\n", cycle);
        return false;
    }
    if (!report.preserve_value_criterion) {
        std::printf("[diagnose_bootstrap_one_cycle_precision] cycle=%zu status=BLOCKED first_failing_stage=preserve_value_criterion\n", cycle);
        return false;
    }
    if (!report.restore_level_criterion) {
        std::printf("[diagnose_bootstrap_one_cycle_precision] cycle=%zu status=BLOCKED first_failing_stage=restore_level_criterion\n", cycle);
        return false;
    }
    if (final_max_error > kTolerance) {
        std::printf("[diagnose_bootstrap_one_cycle_precision] cycle=%zu status=BLOCKED first_failing_stage=final_error max_error=%.12e tolerance=%.12e\n",
                    cycle,
                    final_max_error,
                    kTolerance);
        return false;
    }

    current = report.result;
    std::printf("[diagnose_bootstrap_one_cycle_precision] cycle=%zu status=PASS\n", cycle);
    return true;
}

} // namespace

int main() {
    try {
        const auto expected = make_expected();
        const auto rotation_steps = m2424::Bootstrapper::scalable_refresh_rotation_steps(kSlots);
        auto adapter = m2424::SealAdapter::create(m2424::profiles::precision_boot_ultra_ckks_59());
        adapter.keygen(rotation_steps, true);
        auto current = adapter.encrypt(adapter.encode_complex(expected));

        if (!run_cycle(1, adapter, expected, current)) {
            return 0;
        }
        if (!run_cycle(2, adapter, expected, current)) {
            return 0;
        }
        (void)run_cycle(3, adapter, expected, current);
        return 0;
    } catch (const std::exception& error) {
        std::printf("[diagnose_bootstrap_one_cycle_precision] FAIL: %s\n", error.what());
        return 1;
    }
}
