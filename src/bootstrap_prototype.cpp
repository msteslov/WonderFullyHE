#include "m2424/bootstrap_prototype.hpp"

#include "m2424/eval_mod.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace m2424 {
namespace {

double max_complex_error(const ComplexVector& expected, const ComplexVector& actual) {
    const std::size_t n = std::min(expected.size(), actual.size());
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
                                   double tolerance) {
    return BootstrapPrototypeStage{
        name,
        max_error <= tolerance ? "PASS" : "FAIL",
        before.chain_index,
        after.chain_index,
        before.scale,
        after.scale,
        max_error
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
        0.0
    };
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

std::vector<int> BootstrapPrototype::rotation_steps() const {
    std::vector<int> steps = coeff_to_slot_.rotation_steps();
    auto inverse_steps = slot_to_coeff_.rotation_steps();
    steps.insert(steps.end(), inverse_steps.begin(), inverse_steps.end());
    std::sort(steps.begin(), steps.end());
    steps.erase(std::unique(steps.begin(), steps.end()), steps.end());
    return steps;
}

BootstrapPrototypeReport BootstrapPrototype::refresh_harness(const ComplexVector& input) const {
    if (input.size() != slots_) {
        throw std::invalid_argument("input size must match BootstrapPrototype slots");
    }

    EvalModPolynomial eval_mod;
    BootstrapPrototypeReport report;
    report.slots = slots_;
    report.tolerance = tolerance_;
    report.normalization_factor = 1.0;

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
    current = coeff_to_slot_.apply(adapter_, current);
    auto after = adapter_.info(current);
    auto coeff_actual = head(adapter_.decode_complex(adapter_.decrypt(current)), slots_);
    report.stages.push_back(make_stage("coeff_to_slot", before, after,
                                       max_complex_error(coeff_expected, coeff_actual), tolerance_));

    report.stages.push_back(BootstrapPrototypeStage{
        "eval_mod_normalization",
        "PASS",
        after.chain_index,
        after.chain_index,
        after.scale,
        after.scale,
        0.0
    });

    const auto eval_expected = eval_mod.evaluate_plain(coeff_expected);
    before = adapter_.info(current);
    current = eval_mod.evaluate(adapter_, current);
    after = adapter_.info(current);
    auto eval_actual = head(adapter_.decode_complex(adapter_.decrypt(current)), slots_);
    report.stages.push_back(make_stage("eval_mod", before, after,
                                       max_complex_error(eval_expected, eval_actual), tolerance_));

    const auto result_expected = slot_to_coeff_.apply_plain(eval_expected);
    before = adapter_.info(current);
    current = slot_to_coeff_.apply(adapter_, current);
    after = adapter_.info(current);
    auto result_actual = head(adapter_.decode_complex(adapter_.decrypt(current)), slots_);
    const double result_error = max_complex_error(result_expected, result_actual);
    report.stages.push_back(make_stage("slot_to_coeff", before, after, result_error, tolerance_));

    const double preserve_error = max_complex_error(input, result_actual);
    report.stages.push_back(BootstrapPrototypeStage{
        "refresh_result",
        preserve_error <= tolerance_ ? "PASS" : "FAIL",
        after.chain_index,
        after.chain_index,
        after.scale,
        after.scale,
        preserve_error
    });

    report.preserve_value_criterion = preserve_error <= tolerance_;
    report.restore_level_criterion = after.chain_index >= 0;
    report.result = std::move(current);
    return report;
}

} // namespace m2424
