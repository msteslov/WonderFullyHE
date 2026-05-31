#include "m2424/bootstrap_dft.hpp"
#include "m2424/bootstrap_plan.hpp"
#include "m2424/eval_mod.hpp"
#include "m2424/security_report.hpp"
#include "m2424/seal_adapter.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>
#include <string>
#include <vector>

namespace {

constexpr std::size_t kSlots = 16;
constexpr double kInputAmplitude = 1e-5;
constexpr double kTransformScaleLog2 = 40.0;
constexpr double kNormalizationPlainScaleLog2 = 160.0;
constexpr double kTargetScaleLog2 = 60.0;
constexpr double kTargetCoeffToSlotMagnitudeLog2 = 80.0;
constexpr std::size_t kLevelDrop = 2;
constexpr std::size_t kEvalModLevels = 3;
constexpr std::size_t kSlotToCoeffLevels = 1;
constexpr std::size_t kResidualLevels = 1;
constexpr double kEvalModMarginLog2 = 2.0;
constexpr double kTolerance = 2e-5;

struct CoeffToSlotMeasurement {
    bool ok{};
    double max_abs_log2{};
    double max_abs{};
    std::size_t chain_index{};
    double scale_log2{};
    double coeff_modulus_log2{};
    std::size_t layers{};
    double transform_scaling_log2{};
    std::string exception;
};

std::vector<double> make_input(std::size_t slots, double amplitude) {
    std::vector<double> input;
    input.reserve(slots);
    for (std::size_t i = 0; i < slots; ++i) {
        input.push_back(amplitude * std::sin(static_cast<double>(i) / 4.0));
    }
    return input;
}

m2424::ComplexVector head(m2424::ComplexVector values, std::size_t n) {
    values.resize(std::min(values.size(), n));
    return values;
}

double max_abs_value(const m2424::ComplexVector& values) {
    double result = 0.0;
    for (const auto& value : values) {
        result = std::max(result, std::abs(value));
    }
    return result;
}

double max_complex_error(const m2424::ComplexVector& expected, const m2424::ComplexVector& actual) {
    const std::size_t n = std::min(expected.size(), actual.size());
    double result = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        result = std::max(result, std::abs(expected[i] - actual[i]));
    }
    return result;
}

m2424::ComplexVector scaled(const m2424::ComplexVector& values, double factor) {
    m2424::ComplexVector result;
    result.reserve(values.size());
    for (const auto& value : values) {
        result.push_back(value * factor);
    }
    return result;
}

void sanitize(std::string& text) {
    std::replace(text.begin(), text.end(), ',', ';');
}

const char* bool_text(bool value) {
    return value ? "true" : "false";
}

m2424::BootstrapLayoutPlanningRequest base_request(double max_abs_log2, double transform_output_scale_log2) {
    m2424::BootstrapLayoutPlanningRequest request;
    request.slots = kSlots;
    request.poly_modulus_degree = 32768;
    request.security_level = m2424::SecurityLevel::TC128;
    request.transform_backend = m2424::BootstrapTransformBackend::FftLike;
    request.max_abs_after_coeff_to_slot_log2 = max_abs_log2;
    request.transform_output_scale_log2 = transform_output_scale_log2;
    request.normalization_plain_scale_log2 = kNormalizationPlainScaleLog2;
    request.target_scale_log2 = kTargetScaleLog2;
    request.evalmod_capacity_margin_log2 = kEvalModMarginLog2;
    request.coeff_to_slot_levels = 0;
    request.evalmod_levels = kEvalModLevels;
    request.slot_to_coeff_levels = kSlotToCoeffLevels;
    request.residual_levels = kResidualLevels;
    return request;
}

std::vector<int> gate_rotation_steps(const m2424::FactorizedLinearTransform& coeff_to_slot) {
    return coeff_to_slot.rotation_steps();
}

m2424::Cipher make_boundary_cipher(m2424::SealAdapter& adapter) {
    auto current = adapter.encrypt(adapter.encode(make_input(kSlots, kInputAmplitude)));
    for (std::size_t i = 0; i < kLevelDrop; ++i) {
        current = adapter.mul_plain_rescale(current, adapter.encode_scalar_like(1.0, current));
    }
    return adapter.mod_raise_to_first(current);
}

CoeffToSlotMeasurement measure_coeff_to_slot(const m2424::BootstrapLayoutPlanningRequest& request,
                                             double transform_scaling_log2) {
    CoeffToSlotMeasurement measurement;
    measurement.transform_scaling_log2 = transform_scaling_log2;
    const auto layout = m2424::plan_bootstrap_layout(request);
    if (layout.status != m2424::BootstrapLayoutPlanningStatus::Ready) {
        measurement.exception = layout.blocker;
        sanitize(measurement.exception);
        return measurement;
    }

    try {
        auto coeff_to_slot = m2424::FactorizedLinearTransform(m2424::make_bootstrap_dft_plan(
            kSlots,
            m2424::BootstrapDftType::HomomorphicDecode,
            kTransformScaleLog2,
            {1},
            transform_scaling_log2));
        auto adapter = m2424::SealAdapter::create(layout.profile);
        adapter.keygen(gate_rotation_steps(coeff_to_slot), true);
        auto current = make_boundary_cipher(adapter);
        std::vector<m2424::BootstrapDftLayerTrace> trace;
        current = coeff_to_slot.apply(adapter, current, &trace, "CoeffToSlot");
        const auto info = adapter.info(current);
        const auto decoded = head(adapter.decode_complex(adapter.decrypt(current)), kSlots);
        measurement.ok = true;
        measurement.max_abs = max_abs_value(decoded);
        measurement.max_abs_log2 = measurement.max_abs > 0.0 ? std::log2(measurement.max_abs) : -INFINITY;
        measurement.chain_index = info.chain_index;
        measurement.scale_log2 = std::log2(info.scale);
        measurement.coeff_modulus_log2 = info.coeff_modulus_log2;
        measurement.layers = trace.size();
    } catch (const std::exception& e) {
        measurement.exception = e.what();
        sanitize(measurement.exception);
    }
    return measurement;
}

bool run_final_gate(const char* label,
                    const m2424::BootstrapLayoutPlanningRequest& request,
                    double transform_scaling_log2) {
    const auto layout = m2424::plan_bootstrap_layout(request);
    bool p3_ready = false;
    bool evalmod_applied = false;
    bool evalmod_error_ok = false;
    std::string exception;
    const char* status = "PASS";

    double actual_max_abs = 0.0;
    double actual_max_abs_log2 = 0.0;
    double normalized_error = 0.0;
    double evalmod_error = 0.0;
    std::size_t cts_chain = 0;
    std::size_t norm_chain = 0;
    std::size_t squash_chain = 0;
    std::size_t evalmod_chain = 0;
    double cts_scale_log2 = 0.0;
    double norm_scale_log2 = 0.0;
    double squash_scale_log2 = 0.0;
    double evalmod_scale_log2 = 0.0;
    std::size_t cts_layers = 0;
    std::size_t normalization_chunks = 0;
    std::size_t normalization_levels = 0;
    std::size_t squash_levels = 0;

    try {
        if (layout.status != m2424::BootstrapLayoutPlanningStatus::Ready) {
            status = "BLOCKED";
            throw std::runtime_error(layout.blocker);
        }

        auto coeff_to_slot = m2424::FactorizedLinearTransform(m2424::make_bootstrap_dft_plan(
            kSlots,
            m2424::BootstrapDftType::HomomorphicDecode,
            kTransformScaleLog2,
            {1},
            transform_scaling_log2));
        auto adapter = m2424::SealAdapter::create(layout.profile);
        adapter.keygen(gate_rotation_steps(coeff_to_slot), true);

        auto current = make_boundary_cipher(adapter);
        std::vector<m2424::BootstrapDftLayerTrace> trace;
        current = coeff_to_slot.apply(adapter, current, &trace, "CoeffToSlot");
        cts_layers = trace.size();
        const auto cts_info = adapter.info(current);
        cts_chain = cts_info.chain_index;
        cts_scale_log2 = std::log2(cts_info.scale);
        const auto before_normalization = head(adapter.decode_complex(adapter.decrypt(current)), kSlots);
        actual_max_abs = max_abs_value(before_normalization);
        actual_max_abs_log2 = actual_max_abs > 0.0 ? std::log2(actual_max_abs) : -INFINITY;

        const auto expected_normalized = scaled(before_normalization, std::exp2(-layout.period_log2));
        auto normalized = m2424::apply_bootstrap_scalar_decomposed(
            adapter, current, -layout.period_log2, request.normalization_plain_scale_log2);
        normalization_chunks = normalized.chunks;
        normalization_levels = normalized.levels_consumed;
        const auto normalized_info = adapter.info(normalized.result);
        norm_chain = normalized_info.chain_index;
        norm_scale_log2 = std::log2(normalized_info.scale);
        const auto actual_normalized = head(adapter.decode_complex(adapter.decrypt(normalized.result)), kSlots);
        normalized_error = max_complex_error(expected_normalized, actual_normalized);

        constexpr std::size_t min_chain_before_evalmod = kEvalModLevels + kSlotToCoeffLevels + kResidualLevels;
        auto squashed = m2424::squash_bootstrap_scale(
            adapter, normalized.result, kTargetScaleLog2, min_chain_before_evalmod);
        squash_levels = squashed.levels_consumed;
        const auto squash_info = adapter.info(squashed.result);
        squash_chain = squash_info.chain_index;
        squash_scale_log2 = std::log2(squash_info.scale);
        const bool scalar_pass = normalized_error <= kTolerance;
        const bool interval_ready = actual_max_abs * std::exp2(-layout.period_log2)
            <= m2424::EvalModPolynomial::approximation_bound;
        const bool squash_ready = squash_scale_log2 <= kTargetScaleLog2
            && squash_chain >= min_chain_before_evalmod;
        p3_ready = scalar_pass && interval_ready && squash_ready;

        if (p3_ready) {
            m2424::EvalModPolynomial eval_mod;
            const auto expected_evalmod = [&] {
                m2424::ComplexVector result;
                result.reserve(actual_normalized.size());
                for (const auto& value : actual_normalized) {
                    result.push_back(eval_mod.evaluate_plain(value, m2424::EvalModDegree::P3));
                }
                return result;
            }();
            auto after_evalmod = eval_mod.evaluate(adapter, squashed.result, m2424::EvalModDegree::P3);
            const auto evalmod_info = adapter.info(after_evalmod);
            evalmod_chain = evalmod_info.chain_index;
            evalmod_scale_log2 = std::log2(evalmod_info.scale);
            const auto actual_evalmod = head(adapter.decode_complex(adapter.decrypt(after_evalmod)), kSlots);
            evalmod_error = max_complex_error(expected_evalmod, actual_evalmod);
            evalmod_applied = true;
            evalmod_error_ok = evalmod_error <= kTolerance;
        }
    } catch (const std::exception& e) {
        exception = e.what();
        sanitize(exception);
        if (status != std::string("BLOCKED")) {
            status = "FAIL";
        }
    }

    std::printf("%s,%s,%zu,%zu,%d,%d,%.6f,%.6f,%.6f,%zu,%zu,%zu,%zu,%zu,%zu,%.6f,%.6f,%.6e,%.6f,%zu,%.6f,%zu,%zu,%zu,%.6f,%zu,%zu,%.6f,%zu,%.6f,%.6e,%s,%s,%s,%.6e,%s,%s,%s\n",
                label,
                m2424::to_string(layout.status),
                layout.slots,
                layout.poly_modulus_degree,
                layout.total_coeff_modulus_bits,
                layout.security_budget_bits,
                request.max_abs_after_coeff_to_slot_log2,
                transform_scaling_log2,
                layout.period_log2,
                layout.coeff_to_slot_levels,
                layout.normalization_levels,
                layout.scale_squash_levels,
                layout.evalmod_levels,
                layout.slot_to_coeff_levels,
                layout.total_levels,
                layout.scale_after_normalization_log2,
                layout.scale_after_squash_log2,
                actual_max_abs,
                actual_max_abs_log2,
                cts_chain,
                cts_scale_log2,
                cts_layers,
                normalization_chunks,
                normalization_levels,
                norm_scale_log2,
                squash_chain,
                squash_levels,
                squash_scale_log2,
                evalmod_chain,
                evalmod_scale_log2,
                normalized_error,
                bool_text(p3_ready),
                bool_text(evalmod_applied),
                bool_text(evalmod_error_ok),
                evalmod_error,
                exception.c_str(),
                status,
                layout.blocker.c_str());
    return p3_ready && evalmod_applied && evalmod_error_ok;
}

} // namespace

int main() {
    std::printf("case,layout_status,slots,poly_modulus_degree,total_coeff_modulus_bits,security_budget_bits,planned_max_abs_log2,transform_scaling_log2,period_log2,coeff_to_slot_levels,normalization_levels,scale_squash_levels,evalmod_levels,slot_to_coeff_levels,total_levels,planned_scale_after_normalization_log2,planned_scale_after_squash_log2,actual_max_abs_after_coeff_to_slot,actual_max_abs_log2,chain_after_coeff_to_slot,scale_after_coeff_to_slot_log2,coeff_to_slot_layers,normalization_chunks,normalization_levels_consumed,scale_after_normalization_log2,chain_after_squash,scale_squash_levels_consumed,scale_after_squash_log2,chain_after_evalmod,scale_after_evalmod_log2,normalization_error,p3_ready,evalmod_applied,evalmod_error_ok,evalmod_error,exception,status,blocker\n");

    const auto initial = base_request(132.0, kTransformScaleLog2);
    const auto measurement = measure_coeff_to_slot(initial, 0.0);
    if (!measurement.ok) {
        std::fprintf(stderr, "CoeffToSlot measurement failed: %s\n", measurement.exception.c_str());
        return 1;
    }

    const double transform_scaling_log2 =
        std::floor(kTargetCoeffToSlotMagnitudeLog2 - measurement.max_abs_log2);
    const auto scaled_measurement = measure_coeff_to_slot(initial, transform_scaling_log2);
    if (!scaled_measurement.ok) {
        std::fprintf(stderr, "Scaled CoeffToSlot measurement failed: %s\n",
                     scaled_measurement.exception.c_str());
        return 1;
    }

    const double planned_max_abs_log2 = std::ceil(scaled_measurement.max_abs_log2);
    const auto final_request = base_request(planned_max_abs_log2, scaled_measurement.scale_log2);
    (void)run_final_gate("measured_fft_gate", final_request, transform_scaling_log2);
    return 0;
}
