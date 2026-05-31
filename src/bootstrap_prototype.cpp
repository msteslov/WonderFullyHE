#include "m2424/bootstrap_prototype.hpp"

#include "bootstrap_prototype_detail.hpp"
#include "m2424/eval_mod.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace m2424 {
using namespace bootstrap_prototype_detail;

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
    if (transform_backend_ == BootstrapTransformBackend::FftLike) {
        auto coeff_to_slot = make_bootstrap_dft_plan(slots_, BootstrapDftType::HomomorphicDecode, 40.0);
        auto slot_to_coeff = make_bootstrap_dft_plan(slots_, BootstrapDftType::HomomorphicEncode, 40.0);
        auto steps = coeff_to_slot.rotation_steps();
        auto inverse_steps = slot_to_coeff.rotation_steps();
        steps.insert(steps.end(), inverse_steps.begin(), inverse_steps.end());
        std::sort(steps.begin(), steps.end());
        steps.erase(std::unique(steps.begin(), steps.end()), steps.end());
        return steps;
    }
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

void BootstrapPrototype::set_circuit_order(BootstrapCircuitOrder order) noexcept {
    circuit_order_ = order;
}

void BootstrapPrototype::set_transform_backend(BootstrapTransformBackend backend) noexcept {
    transform_backend_ = backend;
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

void BootstrapPrototype::set_output_correction_factor(double value) {
    if (!std::isfinite(value) || value <= 0.0) {
        throw std::invalid_argument("bootstrap output correction factor must be positive and finite");
    }
    output_correction_factor_ = value;
}

void BootstrapPrototype::set_post_refresh_mod_raise_enabled(bool enabled) noexcept {
    post_refresh_mod_raise_enabled_ = enabled;
}

void BootstrapPrototype::set_stc_first_target_chain_index(std::size_t value) noexcept {
    stc_first_target_chain_index_ = value;
}

void BootstrapPrototype::set_stc_first_period_offset_log2(double value) {
    if (!std::isfinite(value)) {
        throw std::invalid_argument("STC-first period offset log2 must be finite");
    }
    stc_first_period_offset_log2_ = value;
}

Cipher BootstrapPrototype::apply_normalization(const Cipher& input, double factor) const {
    if (normalization_mode_ == BootstrapNormalizationMode::ScaleReinterpretation) {
        return adapter_.unsafe_reinterpret_scale_for_diagnostics(input, factor);
    }
    if (!std::isfinite(factor) || factor == 0.0) {
        throw std::invalid_argument("bootstrap normalization factor must be non-zero and finite");
    }
    return apply_bootstrap_scalar_decomposed(
        adapter_, input, std::log2(std::abs(factor)), plain_scale_log2_).result;
}



} // namespace m2424
