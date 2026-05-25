#include "m2424/bootstrap_prototype.hpp"

#include "m2424/eval_mod.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
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

BootstrapPrototype::BootstrapPrototype(SealAdapter& adapter, std::size_t slots, double tolerance)
    : adapter_(adapter),
      slots_(slots),
      tolerance_(tolerance),
      coeff_to_slot_(DiagonalLinearTransform::from_matrix(canonical_embedding_matrix(slots))),
      slot_to_coeff_(DiagonalLinearTransform::from_matrix(invert_matrix(canonical_embedding_matrix(slots)))) {
    if (slots_ == 0) {
        throw std::invalid_argument("slots must be positive");
    }
    if (!std::isfinite(tolerance_) || tolerance_ < 0.0) {
        throw std::invalid_argument("tolerance must be a non-negative finite value");
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

BootstrapPrototypeReport BootstrapPrototype::refresh_impl(const ComplexVector& input, bool checked) const {
    if (input.size() != slots_) {
        throw std::invalid_argument("input size must match BootstrapPrototype slots");
    }

    EvalModPolynomial eval_mod;
    BootstrapPrototypeReport report;
    report.slots = slots_;
    report.tolerance = tolerance_;
    report.normalization_factor = 1.0;
    report.checked = checked;

    ComplexVector packed(adapter_.slot_count(), Complex{0.0, 0.0});
    for (std::size_t i = 0; i < packed.size(); ++i) {
        packed[i] = input[i % slots_];
    }

    auto current = adapter_.encrypt(adapter_.encode_complex(packed));
    report.stages.push_back(make_harness_stage(adapter_.info(current)));

    const auto coeff_expected = coeff_to_slot_.apply_plain(input);
    double max_coeff_abs = 0.0;
    for (const auto& value : coeff_expected) {
        max_coeff_abs = std::max(max_coeff_abs, std::abs(value));
    }
    if (max_coeff_abs > EvalModPolynomial::approximation_bound) {
        throw std::invalid_argument("CoeffToSlot output exceeds EvalMod approximation interval");
    }

    auto before = adapter_.info(current);
    double stage_ms = elapsed_ms([&] {
        current = coeff_to_slot_.apply(adapter_, current);
    });
    auto after = adapter_.info(current);
    ComplexVector coeff_actual;
    double coeff_error = 0.0;
    if (checked) {
        coeff_actual = head(adapter_.decode_complex(adapter_.decrypt(current)), slots_);
        coeff_error = max_complex_error(coeff_expected, coeff_actual);
    }
    report.stages.push_back(make_stage("coeff_to_slot", before, after,
                                       coeff_error, tolerance_, stage_ms, checked));

    report.stages.push_back(BootstrapPrototypeStage{
        "eval_mod_normalization",
        "PASS",
        after.chain_index,
        after.chain_index,
        after.scale,
        after.scale,
        0.0,
        0.0
    });

    ComplexVector eval_expected;
    if (checked) {
        eval_expected = eval_mod.evaluate_plain(coeff_expected);
    }
    before = adapter_.info(current);
    stage_ms = elapsed_ms([&] {
        current = eval_mod.evaluate(adapter_, current);
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
    double result_error = 0.0;
    double preserve_error = 0.0;
    if (checked) {
        result_actual = head(adapter_.decode_complex(adapter_.decrypt(current)), slots_);
        result_error = max_complex_error(result_expected, result_actual);
        preserve_error = max_complex_error(input, result_actual);
    }
    report.stages.push_back(make_stage("slot_to_coeff", before, after, result_error, tolerance_, stage_ms, checked));

    report.stages.push_back(BootstrapPrototypeStage{
        "refresh_result",
        checked ? (preserve_error <= tolerance_ ? "PASS" : "FAIL") : "RUN",
        after.chain_index,
        after.chain_index,
        after.scale,
        after.scale,
        preserve_error,
        0.0
    });

    report.preserve_value_criterion = !checked || preserve_error <= tolerance_;
    report.restore_level_criterion = after.chain_index >= 0;
    report.result = std::move(current);
    return report;
}

} // namespace m2424
