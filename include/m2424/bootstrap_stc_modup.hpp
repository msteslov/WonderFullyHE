#pragma once

#include "m2424/seal_adapter.hpp"

#include <cstddef>
#include <cstdint>
#include <string>

namespace m2424 {

enum class BootstrapStcDomain {
    Slots,
    Coefficients,
    ScaledCoefficients,
    ModRaisedCoefficients,
    PreEvalModSlots
};

struct BootstrapStcScaleDownPlan {
    double input_scale_log2{};
    double target_scale_log2{};
    double scale_down_log2{};
    double expected_period_log2{};
    double expected_message_gain_log2{};
    std::size_t consumed_levels{};
    std::string note;
};

struct BootstrapStcModUpPlan {
    std::size_t slots{};
    double message_scale_log2{};
    double bootstrap_scale_log2{};
    double coeff_to_slot_log_scale{};
    double slots_to_coeff_log_scale{};
    BootstrapStcScaleDownPlan scale_down;
};

struct BootstrapStcScaleDownResult {
    Cipher result;
    std::size_t chain_before{};
    std::size_t chain_after{};
    double scale_before_log2{};
    double scale_after_log2{};
    double coeff_modulus_log2_before{};
    double coeff_modulus_log2_after{};
    std::size_t levels_consumed{};
};

struct BootstrapScaleDownToQPlan {
    double message_scale_log2{};
    double target_scale_log2{};
    double message_ratio_log2{8.0};
    std::size_t target_coeff_modulus_size{1};
    bool preserve_scale_on_level_drop{true};
};

struct BootstrapScaleDownToQResult {
    Cipher result;
    std::size_t chain_before{};
    std::size_t chain_after{};
    std::size_t coeff_modulus_size_before{};
    std::size_t coeff_modulus_size_after{};
    double scale_before_log2{};
    double scale_after_log2{};
    double coeff_modulus_log2_before{};
    double coeff_modulus_log2_after{};
    double current_message_ratio_log2{};
    double target_message_ratio_log2{};
    double scale_up_log2{};
    std::uint64_t scale_up_integer{1};
    double err_scale_log2{};
    std::size_t levels_consumed{};
    std::string note;
};

const char* to_string(BootstrapStcDomain domain) noexcept;

BootstrapStcModUpPlan plan_stc_first_modup(const CipherInfo& after_slots_to_coeff,
                                           double message_scale_log2,
                                           double target_pre_evalmod_scale_log2);

BootstrapStcScaleDownResult apply_stc_scale_down(SealAdapter& adapter,
                                                 const Cipher& coeff_cipher,
                                                 const BootstrapStcScaleDownPlan& plan);

BootstrapScaleDownToQResult bootstrap_scale_down_to_q(SealAdapter& adapter,
                                                      const Cipher& input,
                                                      const BootstrapScaleDownToQPlan& plan);

} // namespace m2424
