#include "m2424/bootstrap_prototype.hpp"

#include "m2424/eval_mod.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace m2424 {
namespace {

double max_complex_error(const ComplexVector& expected, const ComplexVector& actual) {
    if (expected.size() != actual.size()) {
        throw std::invalid_argument("vectors must have equal size");
    }
    const std::size_t n = expected.size();
    double result = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        result = std::max(result, std::abs(expected[i] - actual[i]));
    }
    return result;
}

ComplexVector head(ComplexVector values, std::size_t n) {
    values.resize(std::min(values.size(), n));
    return values;
}

ComplexVector scaled(const ComplexVector& values, double factor) {
    ComplexVector result;
    result.reserve(values.size());
    for (const auto& value : values) {
        result.push_back(value * factor);
    }
    return result;
}

ComplexVector evaluate_plain(const EvalModPolynomial& eval_mod,
                             const ComplexVector& input,
                             EvalModDegree degree) {
    ComplexVector result;
    result.reserve(input.size());
    for (const auto& value : input) {
        result.push_back(eval_mod.evaluate_plain(value, degree));
    }
    return result;
}

double max_abs_value(const ComplexVector& values) {
    double result = 0.0;
    for (const auto& value : values) {
        result = std::max(result, std::abs(value));
    }
    return result;
}

double normalization_factor_for(const ComplexVector& values) {
    const double max_abs = max_abs_value(values);
    if (max_abs == 0.0) {
        return 1.0;
    }
    constexpr double target = EvalModPolynomial::approximation_bound * 0.5;
    return max_abs > target ? target / max_abs : 1.0;
}

BootstrapPrototypeStage make_stage(const std::string& name,
                                   const CipherInfo& before,
                                   const CipherInfo& after,
                                   double max_error,
                                   double tolerance,
                                   double duration_ms,
                                   bool checked = true) {
    return BootstrapPrototypeStage{
        name,
        checked ? (max_error <= tolerance ? "PASS" : "FAIL") : "RUN",
        before.chain_index,
        after.chain_index,
        before.coeff_modulus_size,
        after.coeff_modulus_size,
        before.scale,
        after.scale,
        max_error,
        duration_ms
    };
}

BootstrapPrototypeStage make_harness_stage(const CipherInfo& after) {
    return BootstrapPrototypeStage{
        "mod_raise_harness",
        "PASS",
        after.chain_index,
        after.chain_index,
        after.coeff_modulus_size,
        after.coeff_modulus_size,
        after.scale,
        after.scale,
        0.0,
        0.0
    };
}

template <typename Fn>
double elapsed_ms(Fn&& fn) {
    const auto started = std::chrono::steady_clock::now();
    fn();
    const auto finished = std::chrono::steady_clock::now();
    return std::chrono::duration<double, std::milli>(finished - started).count();
}

} // namespace

const char* to_string(BootstrapNormalizationMode mode) noexcept {
    switch (mode) {
    case BootstrapNormalizationMode::ScaleReinterpretation:
        return "ScaleReinterpretation";
    case BootstrapNormalizationMode::PlainMultiplyRescale:
        return "PlainMultiplyRescale";
    }
    return "unknown";
}

const char* to_string(BootstrapDenormalizationPosition position) noexcept {
    switch (position) {
    case BootstrapDenormalizationPosition::BeforeSlotToCoeff:
        return "BeforeSlotToCoeff";
    case BootstrapDenormalizationPosition::AfterSlotToCoeff:
        return "AfterSlotToCoeff";
    }
    return "unknown";
}

BootstrapPrototype::BootstrapPrototype(SealAdapter& adapter, std::size_t slots, double tolerance)
    : BootstrapPrototype(adapter, slots, tolerance, 1.0) {}

BootstrapPrototype::BootstrapPrototype(SealAdapter& adapter,
                                       std::size_t slots,
                                       double tolerance,
                                       double normalization_factor)
    : adapter_(adapter),
      slots_(slots),
      tolerance_(tolerance),
      normalization_factor_(normalization_factor),
      coeff_to_slot_(DiagonalLinearTransform::from_matrix(canonical_embedding_matrix(slots))),
      slot_to_coeff_(DiagonalLinearTransform::from_matrix(invert_matrix(canonical_embedding_matrix(slots)))) {
    if (slots_ == 0) {
        throw std::invalid_argument("slots must be positive");
    }
    if (!std::isfinite(tolerance_) || tolerance_ < 0.0) {
        throw std::invalid_argument("tolerance must be a non-negative finite value");
    }
    if (!std::isfinite(normalization_factor_) || normalization_factor_ <= 0.0) {
        throw std::invalid_argument("normalization factor must be a positive finite value");
    }
}

std::vector<int> BootstrapPrototype::required_rotation_steps(std::size_t slots) {
    if (slots == 0) {
        throw std::invalid_argument("slots must be positive");
    }
    auto coeff_to_slot = DiagonalLinearTransform::from_matrix(canonical_embedding_matrix(slots));
    auto slot_to_coeff = DiagonalLinearTransform::from_matrix(invert_matrix(canonical_embedding_matrix(slots)));

    std::vector<int> steps = coeff_to_slot.rotation_steps();
    auto inverse_steps = slot_to_coeff.rotation_steps();
    steps.insert(steps.end(), inverse_steps.begin(), inverse_steps.end());
    std::sort(steps.begin(), steps.end());
    steps.erase(std::unique(steps.begin(), steps.end()), steps.end());
    return steps;
}

std::vector<int> BootstrapPrototype::rotation_steps() const {
    return required_rotation_steps(slots_);
}

BootstrapPrototypeReport BootstrapPrototype::refresh_harness(const ComplexVector& input) const {
    return refresh_impl(input, true);
}

BootstrapPrototypeReport BootstrapPrototype::refresh_fast(const ComplexVector& input) const {
    return refresh_impl(input, false);
}

BootstrapPrototypeReport BootstrapPrototype::refresh_cipher_fast(const Cipher& input) const {
    return refresh_cipher_impl(input, nullptr);
}

BootstrapPrototypeReport BootstrapPrototype::refresh_cipher_checked(const Cipher& input,
                                                                    const ComplexVector& expected) const {
    return refresh_cipher_impl(input, &expected);
}

void BootstrapPrototype::set_normalization_mode(BootstrapNormalizationMode mode) noexcept {
    normalization_mode_ = mode;
}

void BootstrapPrototype::set_denormalization_position(BootstrapDenormalizationPosition position) noexcept {
    denormalization_position_ = position;
}

void BootstrapPrototype::set_evalmod_degree(EvalModDegree degree) noexcept {
    evalmod_degree_ = degree;
}

void BootstrapPrototype::set_period_mode(BootstrapPeriodMode mode) noexcept {
    period_mode_ = mode;
}

void BootstrapPrototype::set_manual_period_log2(double value) {
    if (!std::isfinite(value) || value < 0.0) {
        throw std::invalid_argument("manual period log2 must be a non-negative finite value");
    }
    manual_period_log2_ = value;
}

void BootstrapPrototype::set_plain_scale_log2(double value) {
    if (!std::isfinite(value) || value <= 0.0) {
        throw std::invalid_argument("plain scale log2 must be a positive finite value");
    }
    plain_scale_log2_ = value;
}

void BootstrapPrototype::set_post_refresh_mod_raise_enabled(bool enabled) noexcept {
    post_refresh_mod_raise_enabled_ = enabled;
}

Cipher BootstrapPrototype::apply_normalization(const Cipher& input, double factor) const {
    if (normalization_mode_ == BootstrapNormalizationMode::ScaleReinterpretation) {
        return adapter_.multiply_decoded_value(input, factor);
    }
    if (!std::isfinite(factor) || factor == 0.0) {
        throw std::invalid_argument("bootstrap normalization factor must be non-zero and finite");
    }
    const double factor_log2 = std::log2(std::abs(factor));
    const double factor_times_plain_scale_log2 = factor_log2 + plain_scale_log2_;
    if (factor_times_plain_scale_log2 < 0.0) {
        throw std::runtime_error("bootstrap normalization scalar is not representable");
    }
    const auto plain = adapter_.encode_scalar_at_scale_like(
        factor, std::exp2(plain_scale_log2_), input);
    return adapter_.mul_plain_rescale(input, plain);
}

void mark_stage_structural(BootstrapPrototypeStage& stage) {
    stage.status = "STRUCTURAL";
    stage.max_abs_error = 0.0;
}

void mark_stage_diagnostic(BootstrapPrototypeStage& stage) {
    stage.status = "DIAG";
}

double finite_exp2_or_zero(double exponent) {
    if (exponent < -1074.0) {
        return 0.0;
    }
    if (exponent > 1023.0) {
        return std::numeric_limits<double>::infinity();
    }
    return std::exp2(exponent);
}

BootstrapPrototypeReport BootstrapPrototype::refresh_impl(const ComplexVector& input, bool checked) const {
    if (input.size() != slots_) {
        throw std::invalid_argument("input size must match BootstrapPrototype slots");
    }

    EvalModPolynomial eval_mod;
    BootstrapPrototypeReport report;
    report.slots = slots_;
    report.tolerance = tolerance_;
    report.checked = checked;
    report.normalization_mode = normalization_mode_;
    report.denormalization_position = denormalization_position_;
    report.evalmod_degree = evalmod_degree_;
    report.period_mode = period_mode_;
    report.manual_period_log2 = manual_period_log2_;
    report.plain_scale_log2 = plain_scale_log2_;
    report.max_abs_input = max_abs_value(input);

    ComplexVector packed(adapter_.slot_count(), Complex{0.0, 0.0});
    for (std::size_t i = 0; i < slots_; ++i) {
        packed[i] = input[i];
    }

    auto current = adapter_.encrypt(adapter_.encode_complex(packed));
    const auto initial_info = adapter_.info(current);
    report.stages.push_back(make_harness_stage(initial_info));

    const auto coeff_expected_raw = coeff_to_slot_.apply_plain(input);
    const double normalization_factor = normalization_factor_for(coeff_expected_raw);
    report.normalization_factor = normalization_factor;
    report.normalization_factor_log2 = std::log2(normalization_factor);
    report.factor_times_plain_scale_log2 = report.normalization_factor_log2 + plain_scale_log2_;
    report.normalization_scalar_representable = report.factor_times_plain_scale_log2 >= 0.0;
    report.denormalization_scalar_representable = true;
    report.max_abs_after_coeff_to_slot = max_abs_value(coeff_expected_raw);
    const auto coeff_expected = scaled(coeff_expected_raw, normalization_factor);
    report.max_abs_after_normalization = max_abs_value(coeff_expected);
    report.inside_evalmod_interval =
        report.max_abs_after_normalization <= EvalModPolynomial::approximation_bound;

    auto before = adapter_.info(current);
    double stage_ms = elapsed_ms([&] {
        current = coeff_to_slot_.apply(adapter_, current);
    });
    auto after = adapter_.info(current);
    ComplexVector coeff_actual;
    double coeff_error = 0.0;
    if (checked) {
        coeff_actual = head(adapter_.decode_complex(adapter_.decrypt(current)), slots_);
        coeff_error = max_complex_error(coeff_expected_raw, coeff_actual);
    }
    report.stages.push_back(make_stage("coeff_to_slot", before, after,
                                       coeff_error, tolerance_, stage_ms, checked));

    before = after;
    stage_ms = elapsed_ms([&] {
        current = apply_normalization(current, normalization_factor);
    });
    after = adapter_.info(current);
    report.stages.push_back(BootstrapPrototypeStage{
        "eval_mod_normalization",
        checked ? "PASS" : "RUN",
        before.chain_index,
        after.chain_index,
        before.coeff_modulus_size,
        after.coeff_modulus_size,
        before.scale,
        after.scale,
        0.0,
        stage_ms
    });

    ComplexVector eval_expected;
    if (checked) {
        eval_expected = evaluate_plain(eval_mod, coeff_expected, evalmod_degree_);
    }
    before = adapter_.info(current);
    stage_ms = elapsed_ms([&] {
        current = eval_mod.evaluate(adapter_, current, evalmod_degree_);
    });
    after = adapter_.info(current);
    ComplexVector eval_actual;
    double eval_error = 0.0;
    if (checked) {
        eval_actual = head(adapter_.decode_complex(adapter_.decrypt(current)), slots_);
        eval_error = max_complex_error(eval_expected, eval_actual);
    }
    report.stages.push_back(make_stage("eval_mod", before, after,
                                       eval_error, tolerance_, stage_ms, checked));

    double result_error = 0.0;
    double preserve_error = 0.0;

    if (denormalization_position_ == BootstrapDenormalizationPosition::BeforeSlotToCoeff &&
        normalization_factor != 1.0) {
        report.denormalization_scalar_representable =
            std::log2(1.0 / normalization_factor) + plain_scale_log2_ >= 0.0;
        if (checked) {
            eval_expected = scaled(eval_expected, 1.0 / normalization_factor);
        }
        before = adapter_.info(current);
        stage_ms = elapsed_ms([&] {
            current = apply_normalization(current, 1.0 / normalization_factor);
        });
        after = adapter_.info(current);
        result_error = 0.0;
        if (checked) {
            const auto actual = head(adapter_.decode_complex(adapter_.decrypt(current)), slots_);
            result_error = max_complex_error(eval_expected, actual);
        }
        report.stages.push_back(make_stage("refresh_denormalization",
                                           before,
                                           after,
                                           result_error,
                                           tolerance_,
                                           stage_ms,
                                           checked));
    }

    ComplexVector result_expected;
    if (checked) {
        result_expected = slot_to_coeff_.apply_plain(eval_expected);
    }
    before = adapter_.info(current);
    stage_ms = elapsed_ms([&] {
        current = slot_to_coeff_.apply(adapter_, current);
    });
    after = adapter_.info(current);
    ComplexVector result_actual;
    if (checked) {
        result_actual = head(adapter_.decode_complex(adapter_.decrypt(current)), slots_);
        result_error = max_complex_error(result_expected, result_actual);
        preserve_error = max_complex_error(input, result_actual);
    }
    report.stages.push_back(make_stage("slot_to_coeff", before, after, result_error, tolerance_, stage_ms, checked));

    if (denormalization_position_ == BootstrapDenormalizationPosition::AfterSlotToCoeff &&
        normalization_factor != 1.0) {
        report.denormalization_scalar_representable =
            std::log2(1.0 / normalization_factor) + plain_scale_log2_ >= 0.0;
        if (checked) {
            result_expected = scaled(result_expected, 1.0 / normalization_factor);
        }
        before = adapter_.info(current);
        stage_ms = elapsed_ms([&] {
            current = apply_normalization(current, 1.0 / normalization_factor);
        });
        after = adapter_.info(current);
        result_error = 0.0;
        if (checked) {
            result_actual = head(adapter_.decode_complex(adapter_.decrypt(current)), slots_);
            result_error = max_complex_error(result_expected, result_actual);
            preserve_error = max_complex_error(input, result_actual);
        }
        report.stages.push_back(make_stage("refresh_denormalization",
                                           before,
                                           after,
                                           result_error,
                                           tolerance_,
                                           stage_ms,
                                           checked));
    }

    const auto final_info = adapter_.info(current);
    report.stages.push_back(BootstrapPrototypeStage{
        "refresh_result",
        checked ? (preserve_error <= tolerance_ ? "PASS" : "FAIL") : "RUN",
        final_info.chain_index,
        final_info.chain_index,
        final_info.coeff_modulus_size,
        final_info.coeff_modulus_size,
        final_info.scale,
        final_info.scale,
        preserve_error,
        0.0
    });

    report.preserve_value_criterion = !checked || preserve_error <= tolerance_;
    report.restore_level_criterion = final_info.chain_index >= initial_info.chain_index;
    report.result = std::move(current);
    return report;
}

BootstrapPrototypeReport BootstrapPrototype::refresh_cipher_impl(const Cipher& input, const ComplexVector* expected) const {
    if (expected && expected->size() != slots_) {
        throw std::invalid_argument("expected size must match BootstrapPrototype slots");
    }

    EvalModPolynomial eval_mod;
    BootstrapPrototypeReport report;
    report.slots = slots_;
    report.tolerance = tolerance_;
    report.normalization_factor = normalization_factor_;
    report.normalization_mode = normalization_mode_;
    report.denormalization_position = denormalization_position_;
    report.evalmod_degree = evalmod_degree_;
    report.period_mode = period_mode_;
    report.manual_period_log2 = manual_period_log2_;
    report.plain_scale_log2 = plain_scale_log2_;
    report.checked = expected != nullptr;
    if (expected) {
        report.max_abs_input = max_abs_value(*expected);
    }

    const auto input_info = adapter_.info(input);
    report.bootstrap_period_log2 = 0.0;
    report.bootstrap_period = finite_exp2_or_zero(report.bootstrap_period_log2);
    report.bootstrap_scaling_factor = finite_exp2_or_zero(-report.bootstrap_period_log2);
    if (!expected) {
        report.normalization_factor *= report.bootstrap_scaling_factor;
    }
    auto before = input_info;
    Cipher current;
    double stage_ms = elapsed_ms([&] {
        current = adapter_.mod_raise_to_first(input);
    });
    auto after = adapter_.info(current);
    report.bootstrap_period_log2 = bootstrap_period_log2(
        period_mode_, manual_period_log2_, adapter_.coeff_modulus_bits(), input_info, after);
    report.bootstrap_period = finite_exp2_or_zero(report.bootstrap_period_log2);
    report.bootstrap_scaling_factor = finite_exp2_or_zero(-report.bootstrap_period_log2);
    ComplexVector current_expected;
    double stage_error = 0.0;
    if (expected) {
        current_expected = *expected;
        const auto actual = head(adapter_.decode_complex(adapter_.decrypt(current)), slots_);
        report.max_abs_after_mod_raise_decode = max_abs_value(actual);
        stage_error = max_complex_error(current_expected, actual);
        report.mod_raise_diagnostic_error = stage_error;
    }
    BootstrapPrototypeStage mod_raise_stage{
        "mod_raise",
        "STRUCTURAL",
        before.chain_index,
        after.chain_index,
        before.coeff_modulus_size,
        after.coeff_modulus_size,
        before.scale,
        after.scale,
        0.0,
        stage_ms
    };
    report.stages.push_back(mod_raise_stage);

    ComplexVector coeff_expected;
    if (expected) {
        coeff_expected = coeff_to_slot_.apply_plain(current_expected);
        report.max_abs_after_coeff_to_slot = max_abs_value(coeff_expected);
    }
    before = adapter_.info(current);
    stage_ms = elapsed_ms([&] {
        current = coeff_to_slot_.apply(adapter_, current);
    });
    after = adapter_.info(current);
    stage_error = 0.0;
    if (expected) {
        const auto actual = head(adapter_.decode_complex(adapter_.decrypt(current)), slots_);
        stage_error = max_complex_error(coeff_expected, actual);
    }
    auto coeff_to_slot_stage = make_stage("coeff_to_slot",
                                          before,
                                          after,
                                          stage_error,
                                          tolerance_,
                                          stage_ms,
                                          expected != nullptr);
    if (expected) {
        mark_stage_diagnostic(coeff_to_slot_stage);
    }
    report.stages.push_back(coeff_to_slot_stage);

    const double amplitude_normalization_factor = expected
        ? normalization_factor_for(coeff_expected)
        : report.normalization_factor;
    const auto scaling_factors = make_bootstrap_scaling_factors(
        amplitude_normalization_factor, report.bootstrap_period_log2, plain_scale_log2_);
    report.normalization_factor = scaling_factors.factor;
    report.normalization_factor_log2 = scaling_factors.normalization_factor_log2;
    report.plain_scale_log2 = scaling_factors.plain_scale_log2;
    report.factor_times_plain_scale_log2 = scaling_factors.factor_times_plain_scale_log2;
    report.normalization_scalar_representable = scaling_factors.representable;
    report.denormalization_scalar_representable = true;
    const double effective_normalization_factor = report.normalization_factor;
    if (expected) {
        coeff_expected = scaled(coeff_expected, effective_normalization_factor);
        report.max_abs_after_normalization = max_abs_value(coeff_expected);
        report.inside_evalmod_interval =
            report.max_abs_after_normalization <= EvalModPolynomial::approximation_bound;
    }
    before = after;
    stage_ms = elapsed_ms([&] {
        current = apply_normalization(current, effective_normalization_factor);
    });
    after = adapter_.info(current);
    stage_error = 0.0;
    if (expected) {
        const auto actual = head(adapter_.decode_complex(adapter_.decrypt(current)), slots_);
        stage_error = max_complex_error(coeff_expected, actual);
    }
    BootstrapPrototypeStage normalization_stage{
        "eval_mod_normalization",
        expected ? "DIAG" : "RUN",
        before.chain_index,
        after.chain_index,
        before.coeff_modulus_size,
        after.coeff_modulus_size,
        before.scale,
        after.scale,
        stage_error,
        stage_ms
    };
    report.stages.push_back(normalization_stage);

    ComplexVector eval_expected;
    if (expected) {
        eval_expected = evaluate_plain(eval_mod, coeff_expected, evalmod_degree_);
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
    auto eval_mod_stage = make_stage("eval_mod",
                                     before,
                                     after,
                                     stage_error,
                                     tolerance_,
                                     stage_ms,
                                     expected != nullptr);
    if (expected) {
        mark_stage_diagnostic(eval_mod_stage);
    }
    report.stages.push_back(eval_mod_stage);

    double preserve_error = 0.0;
    stage_error = 0.0;

    if (denormalization_position_ == BootstrapDenormalizationPosition::BeforeSlotToCoeff &&
        effective_normalization_factor != 1.0) {
        report.denormalization_scalar_representable =
            std::log2(1.0 / effective_normalization_factor) + plain_scale_log2_ >= 0.0;
        if (expected) {
            eval_expected = scaled(eval_expected, 1.0 / effective_normalization_factor);
        }
        before = adapter_.info(current);
        stage_ms = elapsed_ms([&] {
            current = apply_normalization(current, 1.0 / effective_normalization_factor);
        });
        after = adapter_.info(current);
        stage_error = 0.0;
        if (expected) {
            const auto actual = head(adapter_.decode_complex(adapter_.decrypt(current)), slots_);
            stage_error = max_complex_error(eval_expected, actual);
        }
        auto denorm_stage = make_stage("refresh_denormalization",
                                       before,
                                       after,
                                       stage_error,
                                       tolerance_,
                                       stage_ms,
                                       expected != nullptr);
        if (expected) {
            mark_stage_diagnostic(denorm_stage);
        }
        report.stages.push_back(denorm_stage);
    }

    ComplexVector result_expected;
    if (expected) {
        result_expected = slot_to_coeff_.apply_plain(eval_expected);
    }
    before = adapter_.info(current);
    stage_ms = elapsed_ms([&] {
        current = slot_to_coeff_.apply(adapter_, current);
    });
    after = adapter_.info(current);
    stage_error = 0.0;
    if (expected) {
        auto actual = head(adapter_.decode_complex(adapter_.decrypt(current)), slots_);
        stage_error = max_complex_error(result_expected, actual);
        preserve_error = max_complex_error(*expected, actual);
    }
    auto slot_to_coeff_stage = make_stage("slot_to_coeff",
                                          before,
                                          after,
                                          stage_error,
                                          tolerance_,
                                          stage_ms,
                                          expected != nullptr);
    if (expected) {
        mark_stage_diagnostic(slot_to_coeff_stage);
    }
    report.stages.push_back(slot_to_coeff_stage);

    if (denormalization_position_ == BootstrapDenormalizationPosition::AfterSlotToCoeff &&
        effective_normalization_factor != 1.0) {
        report.denormalization_scalar_representable =
            std::log2(1.0 / effective_normalization_factor) + plain_scale_log2_ >= 0.0;
        if (expected) {
            result_expected = scaled(result_expected, 1.0 / effective_normalization_factor);
        }
        before = adapter_.info(current);
        stage_ms = elapsed_ms([&] {
            current = apply_normalization(current, 1.0 / effective_normalization_factor);
        });
        after = adapter_.info(current);
        stage_error = 0.0;
        if (expected) {
            auto actual = head(adapter_.decode_complex(adapter_.decrypt(current)), slots_);
            stage_error = max_complex_error(result_expected, actual);
            preserve_error = max_complex_error(*expected, actual);
        }
        auto denorm_stage = make_stage("refresh_denormalization",
                                       before,
                                       after,
                                       stage_error,
                                       tolerance_,
                                       stage_ms,
                                       expected != nullptr);
        if (expected) {
            mark_stage_diagnostic(denorm_stage);
        }
        report.stages.push_back(denorm_stage);
    }

    const auto final_info = adapter_.info(current);
    report.stages.push_back(BootstrapPrototypeStage{
        "refresh_result",
        expected ? (preserve_error <= tolerance_ ? "PASS" : "FAIL") : "RUN",
        final_info.chain_index,
        final_info.chain_index,
        final_info.coeff_modulus_size,
        final_info.coeff_modulus_size,
        final_info.scale,
        final_info.scale,
        preserve_error,
        0.0
    });

    report.preserve_value_criterion = expected && preserve_error <= tolerance_;
    report.restore_level_criterion = final_info.chain_index > input_info.chain_index;
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
        }
        BootstrapPrototypeStage post_stage{
            "post_refresh_mod_raise",
            expected ? "DIAG" : "STRUCTURAL",
            before.chain_index,
            after.chain_index,
            before.coeff_modulus_size,
            after.coeff_modulus_size,
            before.scale,
            after.scale,
            stage_error,
            stage_ms
        };
        if (!expected) {
            mark_stage_structural(post_stage);
        }
        report.stages.push_back(post_stage);
    }
    return report;
}

} // namespace m2424
