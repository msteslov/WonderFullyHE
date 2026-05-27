#pragma once

#include "m2424/diagonal_transform.hpp"
#include "m2424/eval_mod.hpp"
#include "m2424/bootstrap_scaling.hpp"
#include "m2424/seal_adapter.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace m2424 {

enum class BootstrapNormalizationMode {
    ScaleReinterpretation,
    PlainMultiplyRescale
};

enum class BootstrapDenormalizationPosition {
    BeforeSlotToCoeff,
    AfterSlotToCoeff
};

struct BootstrapPrototypeStage {
    std::string name;
    std::string status;
    std::size_t chain_before{};
    std::size_t chain_after{};
    std::size_t coeff_modulus_size_before{};
    std::size_t coeff_modulus_size_after{};
    double scale_before{};
    double scale_after{};
    double max_abs_error{};
    double duration_ms{};
};

struct BootstrapPrototypeReport {
    std::size_t slots{};
    double tolerance{};
    double normalization_factor{};
    double bootstrap_period{};
    double bootstrap_period_log2{};
    double bootstrap_scaling_factor{};
    BootstrapPeriodMode period_mode{BootstrapPeriodMode::TotalCoeffModulus};
    double manual_period_log2{};
    double normalization_factor_log2{};
    double plain_scale_log2{40.0};
    double factor_times_plain_scale_log2{};
    bool normalization_scalar_representable{};
    bool denormalization_scalar_representable{};
    double max_abs_after_mod_raise_decode{};
    double mod_raise_diagnostic_error{};
    BootstrapNormalizationMode normalization_mode{BootstrapNormalizationMode::PlainMultiplyRescale};
    BootstrapDenormalizationPosition denormalization_position{BootstrapDenormalizationPosition::AfterSlotToCoeff};
    EvalModDegree evalmod_degree{EvalModDegree::P7};
    double max_abs_input{};
    double max_abs_after_coeff_to_slot{};
    double max_abs_after_normalization{};
    bool inside_evalmod_interval{};
    bool checked{};
    bool preserve_value_criterion{};
    bool restore_level_criterion{};
    std::vector<BootstrapPrototypeStage> stages;
    Cipher result;
};

class BootstrapPrototype {
public:
    BootstrapPrototype(SealAdapter& adapter, std::size_t slots, double tolerance);
    BootstrapPrototype(SealAdapter& adapter, std::size_t slots, double tolerance, double normalization_factor);

    static std::vector<int> required_rotation_steps(std::size_t slots);

    std::vector<int> rotation_steps() const;
    BootstrapPrototypeReport refresh_harness(const ComplexVector& input) const;
    BootstrapPrototypeReport refresh_fast(const ComplexVector& input) const;
    BootstrapPrototypeReport refresh_cipher_fast(const Cipher& input) const;
    BootstrapPrototypeReport refresh_cipher_checked(const Cipher& input, const ComplexVector& expected) const;
    void set_normalization_mode(BootstrapNormalizationMode mode) noexcept;
    void set_denormalization_position(BootstrapDenormalizationPosition position) noexcept;
    void set_evalmod_degree(EvalModDegree degree) noexcept;
    void set_period_mode(BootstrapPeriodMode mode) noexcept;
    void set_manual_period_log2(double value);
    void set_plain_scale_log2(double value);
    void set_post_refresh_mod_raise_enabled(bool enabled) noexcept;

private:
    BootstrapPrototypeReport refresh_impl(const ComplexVector& input, bool checked) const;
    BootstrapPrototypeReport refresh_cipher_impl(const Cipher& input, const ComplexVector* expected) const;
    Cipher apply_normalization(const Cipher& input, double factor) const;

    SealAdapter& adapter_;
    std::size_t slots_{};
    double tolerance_{};
    double normalization_factor_{};
    BootstrapNormalizationMode normalization_mode_{BootstrapNormalizationMode::PlainMultiplyRescale};
    BootstrapDenormalizationPosition denormalization_position_{BootstrapDenormalizationPosition::AfterSlotToCoeff};
    EvalModDegree evalmod_degree_{EvalModDegree::P7};
    BootstrapPeriodMode period_mode_{BootstrapPeriodMode::TotalCoeffModulus};
    double manual_period_log2_{0.0};
    double plain_scale_log2_{40.0};
    bool post_refresh_mod_raise_enabled_{false};
    DiagonalLinearTransform coeff_to_slot_;
    DiagonalLinearTransform slot_to_coeff_;
};

const char* to_string(BootstrapNormalizationMode mode) noexcept;
const char* to_string(BootstrapDenormalizationPosition position) noexcept;

} // namespace m2424
