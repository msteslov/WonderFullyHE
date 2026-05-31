#include "m2424/bootstrap_prototype.hpp"

#include "bootstrap_prototype_detail.hpp"
#include "m2424/eval_mod.hpp"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace m2424 {
using namespace bootstrap_prototype_detail;

BootstrapPrototypeReport BootstrapPrototype::refresh_cipher_slots_to_coeffs_first_impl(
    const Cipher& input,
    const ComplexVector* expected) const {
    if (expected && expected->size() != slots_) {
        throw std::invalid_argument("expected size must match BootstrapPrototype slots");
    }
    if (evalmod_degree_ != EvalModDegree::P3 && evalmod_degree_ != EvalModDegree::P3DoubleAngle) {
        throw std::invalid_argument("SlotsToCoeffsFirst prototype supports only EvalModDegree::P3/P3DoubleAngle");
    }

    EvalModPolynomial eval_mod;
    BootstrapPrototypeReport report;
    report.slots = slots_;
    report.tolerance = tolerance_;
    report.normalization_factor = 1.0;
    report.normalization_mode = normalization_mode_;
    report.denormalization_position = denormalization_position_;
    report.evalmod_degree = evalmod_degree_;
    report.period_mode = period_mode_;
    report.circuit_order = circuit_order_;
    report.transform_backend = transform_backend_;
    report.stc_first_target_chain_index = stc_first_target_chain_index_;
    report.stc_first_period_offset_log2 = stc_first_period_offset_log2_;
    report.manual_period_log2 = manual_period_log2_;
    report.plain_scale_log2 = plain_scale_log2_;
    report.checked = expected != nullptr;
    if (expected) {
        report.max_abs_input = max_abs_value(*expected);
    }

    auto slot_to_coeff = FactorizedLinearTransform(make_bootstrap_dft_plan(
        slots_, BootstrapDftType::HomomorphicEncode, plain_scale_log2_));
    auto coeff_to_slot = FactorizedLinearTransform(make_bootstrap_dft_plan(
        slots_, BootstrapDftType::HomomorphicDecode, plain_scale_log2_));

    auto current = input;
    const auto input_info = adapter_.info(input);
    auto before = input_info;
    auto after = before;
    double stage_ms = 0.0;
    ComplexVector current_expected = expected ? *expected : ComplexVector{};

    while (adapter_.info(current).chain_index > stc_first_target_chain_index_) {
        current = adapter_.mul_plain_rescale(current, adapter_.encode_scalar_like(1.0, current));
    }
    after = adapter_.info(current);
    report.stages.push_back(BootstrapPrototypeStage{
        "stc_first_level_drop",
        expected ? "DIAG" : "RUN",
        before.chain_index,
        after.chain_index,
        before.coeff_modulus_size,
        after.coeff_modulus_size,
        before.coeff_modulus_log2,
        after.coeff_modulus_log2,
        before.scale,
        after.scale,
        0.0,
        0.0
    });

    ComplexVector coeff_expected;
    if (expected) {
        coeff_expected = slot_to_coeff.apply_plain(current_expected);
    }
    before = adapter_.info(current);
    stage_ms = elapsed_ms([&] {
        current = slot_to_coeff.apply(adapter_, current);
    });
    after = adapter_.info(current);
    double stage_error = 0.0;
    if (expected) {
        const auto actual = head(adapter_.decode_complex(adapter_.decrypt(current)), slots_);
        stage_error = max_complex_error(coeff_expected, actual);
    }
    auto stc_stage = make_stage("slot_to_coeff_first",
                                before,
                                after,
                                stage_error,
                                tolerance_,
                                stage_ms,
                                expected != nullptr);
    if (expected) {
        mark_stage_diagnostic(stc_stage);
    }
    report.stages.push_back(stc_stage);

    before = adapter_.info(current);
    stage_ms = elapsed_ms([&] {
        current = adapter_.mod_raise_to_first(current);
    });
    after = adapter_.info(current);
    report.bootstrap_period_log2 =
        period_mode_ == BootstrapPeriodMode::NoBootstrapPeriod
            ? 0.0
            : before.coeff_modulus_log2 - stc_first_period_offset_log2_;
    report.bootstrap_period = finite_exp2_or_zero(report.bootstrap_period_log2);
    report.bootstrap_scaling_factor = finite_exp2_or_zero(-report.bootstrap_period_log2);
    report.stages.push_back(BootstrapPrototypeStage{
        "mod_raise",
        "STRUCTURAL",
        before.chain_index,
        after.chain_index,
        before.coeff_modulus_size,
        after.coeff_modulus_size,
        before.coeff_modulus_log2,
        after.coeff_modulus_log2,
        before.scale,
        after.scale,
        0.0,
        stage_ms
    });

    ComplexVector roundtrip_expected;
    if (expected) {
        roundtrip_expected = coeff_to_slot.apply_plain(coeff_expected);
    }
    before = adapter_.info(current);
    stage_ms = elapsed_ms([&] {
        current = coeff_to_slot.apply(adapter_, current);
    });
    after = adapter_.info(current);
    if (expected) {
        const auto actual = head(adapter_.decode_complex(adapter_.decrypt(current)), slots_);
        report.max_abs_after_coeff_to_slot = max_abs_value(actual);
        stage_error = max_complex_error(roundtrip_expected, actual);
        if (period_mode_ == BootstrapPeriodMode::NoBootstrapPeriod
            && report.max_abs_after_coeff_to_slot > EvalModPolynomial::approximation_bound) {
            report.bootstrap_period_log2 = std::ceil(
                std::log2(report.max_abs_after_coeff_to_slot / (0.5 * EvalModPolynomial::approximation_bound)));
        }
    }
    report.bootstrap_period = finite_exp2_or_zero(report.bootstrap_period_log2);
    report.bootstrap_scaling_factor = finite_exp2_or_zero(-report.bootstrap_period_log2);
    auto cts_stage = make_stage("coeff_to_slot_after_raise",
                                before,
                                after,
                                stage_error,
                                tolerance_,
                                stage_ms,
                                expected != nullptr);
    if (expected) {
        mark_stage_diagnostic(cts_stage);
    }
    report.stages.push_back(cts_stage);

    ComplexVector normalized_expected;
    if (expected) {
        normalized_expected = head(adapter_.decode_complex(adapter_.decrypt(current)), slots_);
        for (auto& value : normalized_expected) {
            value *= report.bootstrap_scaling_factor;
        }
        report.max_abs_after_normalization = max_abs_value(normalized_expected);
        report.inside_evalmod_interval =
            report.max_abs_after_normalization <= EvalModPolynomial::approximation_bound;
    }

    report.normalization_factor = report.bootstrap_scaling_factor;
    report.normalization_factor_log2 = -report.bootstrap_period_log2;
    report.factor_times_plain_scale_log2 = report.normalization_factor_log2 + plain_scale_log2_;
    report.normalization_scalar_representable = true;
    report.denormalization_scalar_representable = true;

    before = adapter_.info(current);
    stage_ms = elapsed_ms([&] {
        current = apply_normalization(current, report.normalization_factor);
    });
    after = adapter_.info(current);
    stage_error = 0.0;
    if (expected) {
        const auto actual = head(adapter_.decode_complex(adapter_.decrypt(current)), slots_);
        stage_error = max_complex_error(normalized_expected, actual);
    }
    auto norm_stage = make_stage("eval_mod_normalization",
                                 before,
                                 after,
                                 stage_error,
                                 tolerance_,
                                 stage_ms,
                                 expected != nullptr);
    if (expected) {
        mark_stage_diagnostic(norm_stage);
    }
    report.stages.push_back(norm_stage);

    ComplexVector eval_expected;
    if (expected) {
        eval_expected = evaluate_plain(eval_mod, normalized_expected, evalmod_degree_);
    }
    before = adapter_.info(current);
    stage_ms = elapsed_ms([&] {
        current = eval_mod.evaluate(adapter_, current, evalmod_degree_);
    });
    after = adapter_.info(current);
    stage_error = 0.0;
    if (expected) {
        const auto actual = head(adapter_.decode_complex(adapter_.decrypt(current)), slots_);
        stage_error = max_complex_error(eval_expected, actual);
    }
    auto eval_stage = make_stage("eval_mod",
                                 before,
                                 after,
                                 stage_error,
                                 tolerance_,
                                 stage_ms,
                                 expected != nullptr);
    if (expected) {
        mark_stage_diagnostic(eval_stage);
    }
    report.stages.push_back(eval_stage);

    if (output_correction_factor_ != 1.0) {
        if (expected) {
            eval_expected = scaled(eval_expected, output_correction_factor_);
        }
        before = adapter_.info(current);
        stage_ms = elapsed_ms([&] {
            current = apply_normalization(current, output_correction_factor_);
        });
        after = adapter_.info(current);
        stage_error = 0.0;
        if (expected) {
            const auto actual = head(adapter_.decode_complex(adapter_.decrypt(current)), slots_);
            stage_error = max_complex_error(eval_expected, actual);
        }
        auto correction_stage = make_stage("output_correction",
                                           before,
                                           after,
                                           stage_error,
                                           tolerance_,
                                           stage_ms,
                                           expected != nullptr);
        if (expected) {
            mark_stage_diagnostic(correction_stage);
        }
        report.stages.push_back(correction_stage);
    }

    before = adapter_.info(current);
    stage_ms = elapsed_ms([&] {
        current = apply_output_scale_repair(adapter_, current, plain_scale_log2_);
    });
    after = adapter_.info(current);
    stage_error = 0.0;
    if (expected) {
        const auto actual = head(adapter_.decode_complex(adapter_.decrypt(current)), slots_);
        stage_error = max_complex_error(*expected, actual);
    }
    if (after.chain_index != before.chain_index || std::fabs(std::log2(after.scale) - std::log2(before.scale)) > 0.5) {
        auto repair_stage = make_stage("output_scale_repair",
                                       before,
                                       after,
                                       stage_error,
                                       tolerance_,
                                       stage_ms,
                                       expected != nullptr);
        if (expected) {
            mark_stage_diagnostic(repair_stage);
        }
        report.stages.push_back(repair_stage);
    }

    double preserve_error = 0.0;
    if (expected) {
        const auto actual = head(adapter_.decode_complex(adapter_.decrypt(current)), slots_);
        preserve_error = max_complex_error(*expected, actual);
    }
    const auto final_info = adapter_.info(current);
    report.stages.push_back(BootstrapPrototypeStage{
        "refresh_result",
        expected ? (preserve_error <= tolerance_ ? "PASS" : "FAIL") : "RUN",
        final_info.chain_index,
        final_info.chain_index,
        final_info.coeff_modulus_size,
        final_info.coeff_modulus_size,
        final_info.coeff_modulus_log2,
        final_info.coeff_modulus_log2,
        final_info.scale,
        final_info.scale,
        preserve_error,
        0.0
    });

    report.preserve_value_criterion = expected && preserve_error <= tolerance_;
    report.restore_level_criterion = final_info.chain_index >= input_info.chain_index;
    report.continuation_levels = final_info.chain_index;
    report.result = std::move(current);
    if (post_refresh_mod_raise_enabled_) {
        before = adapter_.info(report.result);
        Cipher post_refresh_result;
        stage_ms = elapsed_ms([&] {
            post_refresh_result = adapter_.mod_raise_to_first(report.result);
        });
        after = adapter_.info(post_refresh_result);
        stage_error = 0.0;
        if (expected) {
            const auto actual = head(adapter_.decode_complex(adapter_.decrypt(post_refresh_result)), slots_);
            stage_error = max_complex_error(*expected, actual);
            preserve_error = stage_error;
        }
        BootstrapPrototypeStage post_stage{
            "post_refresh_mod_raise",
            expected ? (stage_error <= tolerance_ ? "PASS" : "FAIL") : "STRUCTURAL",
            before.chain_index,
            after.chain_index,
            before.coeff_modulus_size,
            after.coeff_modulus_size,
            before.coeff_modulus_log2,
            after.coeff_modulus_log2,
            before.scale,
            after.scale,
            stage_error,
            stage_ms
        };
        report.stages.push_back(post_stage);
        report.result = std::move(post_refresh_result);
        report.preserve_value_criterion = !expected || preserve_error <= tolerance_;
        report.restore_level_criterion = after.chain_index >= input_info.chain_index;
        report.continuation_levels = after.chain_index;
    }
    return report;
}

} // namespace m2424
