#include "m2424/bootstrap.hpp"
#include "m2424/profiles.hpp"
#include "m2424/seal_adapter.hpp"

#include <cmath>
#include <cstdio>
#include <exception>
#include <stdexcept>
#include <vector>

namespace {

constexpr double kTargetError = 1e-9;

m2424::ComplexVector make_expected(std::size_t slots, double amplitude) {
    m2424::ComplexVector expected;
    expected.reserve(slots);
    for (std::size_t i = 0; i < slots; ++i) {
        const double x = static_cast<double>(i + 1);
        const double real_sign = (i % 2 == 0) ? 1.0 : -1.0;
        const double imag_sign = (i % 3 == 0) ? -1.0 : 1.0;
        expected.push_back({
            real_sign * amplitude * (0.25 + 0.1 * std::sin(x)),
            imag_sign * amplitude * (0.15 + 0.05 * std::cos(0.5 * x))
        });
    }
    return expected;
}

void require(bool condition, const char* message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

void print_stage_report(const m2424::BootstrapPrototypeReport& report) {
    for (const auto& stage : report.stages) {
        const double before_scale_log2 = stage.scale_before > 0.0 ? std::log2(stage.scale_before) : 0.0;
        const double after_scale_log2 = stage.scale_after > 0.0 ? std::log2(stage.scale_after) : 0.0;
        std::printf("[test_bootstrap_multi_cycle_precision] stage name=%s status=%s chain_before=%zu chain_after=%zu scale_before_log2=%.6f scale_after_log2=%.6f max_abs_error=%.6e\n",
                    stage.name.c_str(),
                    stage.status.c_str(),
                    stage.chain_before,
                    stage.chain_after,
                    before_scale_log2,
                    after_scale_log2,
                    stage.max_abs_error);
    }
}

std::pair<double, double> errors(const m2424::ComplexVector& expected,
                                 const m2424::ComplexVector& actual,
                                 std::size_t slots) {
    double max_error = 0.0;
    double sum_error = 0.0;
    for (std::size_t i = 0; i < slots; ++i) {
        const double error = std::abs(actual[i] - expected[i]);
        max_error = std::max(max_error, error);
        sum_error += error;
    }
    return {max_error, sum_error / static_cast<double>(slots)};
}

void run_case(std::size_t slots, std::size_t cycles, double amplitude) {
    const auto rotations = m2424::Bootstrapper::scalable_refresh_rotation_steps(slots);
    require(!rotations.empty(), "bootstrap rotations must not be empty");

    auto adapter = m2424::SealAdapter::create(m2424::profiles::boot_deep_ckks());
    adapter.keygen(rotations, true);

    const auto expected = make_expected(slots, amplitude);
    auto current = adapter.encrypt(adapter.encode_complex(expected));
    m2424::Bootstrapper bootstrapper(adapter);

    for (std::size_t cycle = 1; cycle <= cycles; ++cycle) {
        while (adapter.info(current).chain_index > 2) {
            current = adapter.mul_plain_rescale(current, adapter.encode_scalar_like(1.0, current));
        }

        const auto input_info = adapter.info(current);
        const auto refresh = bootstrapper.refresh_slots_to_coeffs_first_checked(
            current, expected, slots, kTargetError, m2424::EvalModDegree::P3DoubleAngle);
        const auto output_info = adapter.info(refresh.result);
        const auto actual = adapter.decode_complex(adapter.decrypt(refresh.result));
        const auto [max_error, mean_error] = errors(expected, actual, slots);
        for (std::size_t i = 0; i < slots; ++i) {
            std::printf("[test_bootstrap_multi_cycle_precision] value cycle=%zu index=%zu expected=(%.12e,%.12e) actual=(%.12e,%.12e) abs_error=%.12e\n",
                        cycle,
                        i,
                        expected[i].real(),
                        expected[i].imag(),
                        actual[i].real(),
                        actual[i].imag(),
                        std::abs(actual[i] - expected[i]));
        }

        std::printf("[test_bootstrap_multi_cycle_precision] cycle=%zu slots=%zu max_error=%.12e mean_error=%.12e input_chain=%zu output_chain=%zu continuation_levels=%zu scale_log2=%.6f stages=%zu\n",
                    cycle,
                    slots,
                    max_error,
                    mean_error,
                    input_info.chain_index,
                    output_info.chain_index,
                    refresh.continuation_levels,
                    std::log2(output_info.scale),
                    refresh.stages.size());
        std::printf("[test_bootstrap_multi_cycle_precision] cycle=%zu bootstrap_period_log2=%.6f normalization_factor_log2=%.6f max_abs_input=%.12e max_abs_after_coeff_to_slot=%.12e max_abs_after_normalization=%.12e inside_evalmod=%s preserve=%s restore=%s\n",
                    cycle,
                    refresh.bootstrap_period_log2,
                    refresh.normalization_factor_log2,
                    refresh.max_abs_input,
                    refresh.max_abs_after_coeff_to_slot,
                    refresh.max_abs_after_normalization,
                    refresh.inside_evalmod_interval ? "true" : "false",
                    refresh.preserve_value_criterion ? "true" : "false",
                    refresh.restore_level_criterion ? "true" : "false");
        std::printf("[test_bootstrap_multi_cycle_precision] cycle=%zu max_abs_after_mod_raise_decode=%.12e mod_raise_diagnostic_error=%.12e\n",
                    cycle,
                    refresh.max_abs_after_mod_raise_decode,
                    refresh.mod_raise_diagnostic_error);
        print_stage_report(refresh);

        bool no_failed_stage = true;
        for (const auto& stage : refresh.stages) {
            no_failed_stage = no_failed_stage && stage.status != "FAIL";
        }

        require(no_failed_stage, "refresh contains failing stage");
        require(refresh.preserve_value_criterion, "refresh preserve_value_criterion must hold");
        require(refresh.restore_level_criterion, "refresh restore_level_criterion must hold");
        require(refresh.continuation_levels >= 5, "refresh continuation levels must be at least 5");
        require(refresh.inside_evalmod_interval, "refresh values must stay inside EvalMod interval");
        require(max_error <= kTargetError, "multi-cycle max error exceeded target");

        current = refresh.result;
    }
}

} // namespace

int main() {
    bool ok = true;
    try {
        run_case(4, 3, 1e-5);
    } catch (const std::exception& error) {
        ok = false;
        std::printf("[test_bootstrap_multi_cycle_precision] FAIL: %s\n", error.what());
    }
    return ok ? 0 : 1;
}
