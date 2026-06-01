#include "m2424/bootstrap_prototype.hpp"

#include "bootstrap_prototype_detail.hpp"
#include "m2424/eval_mod.hpp"
#include "m2424/mod1_circuit.hpp"

#include <cmath>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace m2424 {
using namespace bootstrap_prototype_detail;

namespace {

BootstrapMod1Model mod1_model_for(EvalModDegree degree, double plain_scale_log2) {
    switch (degree) {
    case EvalModDegree::P3:
        return {BootstrapMod1Type::LegacySineP3, 3, 0, 8, plain_scale_log2};
    case EvalModDegree::P3DoubleAngle:
        return {BootstrapMod1Type::LegacySineP3, 3, 1, 8, plain_scale_log2};
    case EvalModDegree::P5:
    case EvalModDegree::P7:
        throw std::invalid_argument("SlotsToCoeffsFirst prototype supports only P3/P3DoubleAngle until CosDiscrete encrypted Mod1 is implemented");
    }
    throw std::invalid_argument("unknown EvalMod degree");
}

std::string stage_context(const char* stage, const CipherInfo& before) {
    std::ostringstream out;
    out << "STC-first " << stage
        << " failed"
        << "; before_chain=" << before.chain_index
        << "; before_scale_log2=" << std::log2(before.scale)
        << "; before_coeff_modulus_log2=" << before.coeff_modulus_log2;
    return out.str();
}

Cipher drop_level_preserving_scale(SealAdapter& adapter, const Cipher& input) {
    const auto info = adapter.info(input);
    if (info.chain_index == 0) {
        throw std::runtime_error("cannot drop level at the end of the modulus chain");
    }
    const auto bits = adapter.coeff_modulus_bits();
    if (info.coeff_modulus_size == 0 || info.coeff_modulus_size > bits.size()) {
        throw std::runtime_error("cannot infer next rescale modulus size");
    }
    const double current_scale_log2 = std::log2(info.scale);
    // Encode 1.0 near the next dropped prime so rescale preserves the ciphertext scale.
    const double plain_scale_log2 = static_cast<double>(bits[info.coeff_modulus_size - 1]);
    if (!std::isfinite(current_scale_log2)
        || !std::isfinite(plain_scale_log2)
        || plain_scale_log2 <= 0.0) {
        throw std::runtime_error("cannot compute value-preserving level drop scale");
    }
    return adapter.mul_plain_rescale(
        input,
        adapter.encode_scalar_at_scale_like(1.0, std::exp2(plain_scale_log2), input));
}

} // namespace

BootstrapPrototypeReport BootstrapPrototype::refresh_cipher_slots_to_coeffs_first_impl(
    const Cipher& input,
    const ComplexVector* expected) const {
    if (expected && expected->size() != slots_) {
        throw std::invalid_argument("expected size must match BootstrapPrototype slots");
    }
    if (evalmod_degree_ != EvalModDegree::P3 && evalmod_degree_ != EvalModDegree::P3DoubleAngle) {
        throw std::invalid_argument("SlotsToCoeffsFirst prototype supports only EvalModDegree::P3/P3DoubleAngle");
    }

    Mod1Circuit mod1(mod1_model_for(evalmod_degree_, plain_scale_log2_));
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
    auto decode_head_for_stage = [&](const Cipher& cipher, const char* stage) {
        try {
            return head(adapter_.decode_complex(adapter_.decrypt(cipher)), slots_);
        } catch (const std::exception& e) {
            throw std::runtime_error(stage_context(stage, adapter_.info(cipher)) + "; diagnostic_decrypt=" + e.what());
        }
    };

    try {
        while (adapter_.info(current).chain_index > stc_first_target_chain_index_) {
            current = drop_level_preserving_scale(adapter_, current);
        }
    } catch (const std::exception& e) {
        throw std::runtime_error(stage_context("stc_first_level_drop", adapter_.info(current)) + "; reason=" + e.what());
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
    try {
        stage_ms = elapsed_ms([&] {
            current = slot_to_coeff.apply(adapter_, current);
        });
    } catch (const std::exception& e) {
        throw std::runtime_error(stage_context("slot_to_coeff_first", before) + "; reason=" + e.what());
    }
    after = adapter_.info(current);
    double stage_error = 0.0;
    if (expected) {
        const auto actual = decode_head_for_stage(current, "slot_to_coeff_first");
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
    ComplexVector coeff_before_mod_raise;
    if (expected) {
        coeff_before_mod_raise = decode_head_for_stage(current, "mod_raise_before");
    }
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
    stage_error = 0.0;
    if (expected) {
        const auto coeff_after_mod_raise = decode_head_for_stage(current, "mod_raise");
        stage_error = max_complex_error(coeff_before_mod_raise, coeff_after_mod_raise);
    }
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
        stage_error,
        stage_ms
    });

    ComplexVector roundtrip_expected;
    if (expected) {
        roundtrip_expected = coeff_to_slot.apply_plain(coeff_expected);
    }
    before = adapter_.info(current);
    try {
        stage_ms = elapsed_ms([&] {
            current = coeff_to_slot.apply(adapter_, current);
        });
    } catch (const std::exception& e) {
        throw std::runtime_error(stage_context("coeff_to_slot_after_raise", before) + "; reason=" + e.what());
    }
    after = adapter_.info(current);
    if (expected) {
        const auto actual = decode_head_for_stage(current, "coeff_to_slot_after_raise");
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
        normalized_expected = decode_head_for_stage(current, "eval_mod_normalization_reference");
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
    try {
        stage_ms = elapsed_ms([&] {
            current = apply_normalization(current, report.normalization_factor);
        });
    } catch (const std::exception& e) {
        throw std::runtime_error(stage_context("eval_mod_normalization", before) + "; reason=" + e.what());
    }
    after = adapter_.info(current);
    stage_error = 0.0;
    if (expected) {
        const auto actual = decode_head_for_stage(current, "eval_mod_normalization");
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
        eval_expected = mod1.evaluate_plain(normalized_expected);
    }
    before = adapter_.info(current);
    try {
        stage_ms = elapsed_ms([&] {
            current = mod1.evaluate(adapter_, current);
        });
    } catch (const std::exception& e) {
        throw std::runtime_error(stage_context("eval_mod", before) + "; reason=" + e.what());
    }
    after = adapter_.info(current);
    stage_error = 0.0;
    if (expected) {
        const auto actual = decode_head_for_stage(current, "eval_mod");
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
        try {
            stage_ms = elapsed_ms([&] {
                current = apply_normalization(current, output_correction_factor_);
            });
        } catch (const std::exception& e) {
            throw std::runtime_error(stage_context("output_correction", before) + "; reason=" + e.what());
        }
        after = adapter_.info(current);
        stage_error = 0.0;
        if (expected) {
            const auto actual = decode_head_for_stage(current, "output_correction");
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
    try {
        stage_ms = elapsed_ms([&] {
            current = apply_output_scale_repair(adapter_, current, plain_scale_log2_);
        });
    } catch (const std::exception& e) {
        throw std::runtime_error(stage_context("output_scale_repair", before) + "; reason=" + e.what());
    }
    after = adapter_.info(current);
    stage_error = 0.0;
    if (expected) {
        const auto actual = decode_head_for_stage(current, "output_scale_repair");
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
        const auto actual = decode_head_for_stage(current, "refresh_result");
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
            const auto actual = decode_head_for_stage(post_refresh_result, "post_refresh_mod_raise");
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
