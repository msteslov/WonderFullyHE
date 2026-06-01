#include "m2424/bootstrap_stc_modup.hpp"

#include <cmath>
#include <sstream>
#include <stdexcept>

namespace m2424 {
namespace {

void validate_log2(double value, const char* name) {
    if (!std::isfinite(value)) {
        throw std::invalid_argument(std::string(name) + " must be finite");
    }
}

} // namespace

const char* to_string(BootstrapStcDomain domain) noexcept {
    switch (domain) {
    case BootstrapStcDomain::Slots:
        return "Slots";
    case BootstrapStcDomain::Coefficients:
        return "Coefficients";
    case BootstrapStcDomain::ScaledCoefficients:
        return "ScaledCoefficients";
    case BootstrapStcDomain::ModRaisedCoefficients:
        return "ModRaisedCoefficients";
    case BootstrapStcDomain::PreEvalModSlots:
        return "PreEvalModSlots";
    }
    return "Unknown";
}

BootstrapStcModUpPlan plan_stc_first_modup(const CipherInfo& after_slots_to_coeff,
                                           double message_scale_log2,
                                           double target_pre_evalmod_scale_log2) {
    validate_log2(after_slots_to_coeff.scale > 0.0 ? std::log2(after_slots_to_coeff.scale) : std::nan(""),
                  "after_slots_to_coeff scale log2");
    validate_log2(after_slots_to_coeff.coeff_modulus_log2, "after_slots_to_coeff coeff modulus log2");
    validate_log2(message_scale_log2, "message_scale_log2");
    validate_log2(target_pre_evalmod_scale_log2, "target_pre_evalmod_scale_log2");
    if (after_slots_to_coeff.chain_index == 0) {
        throw std::invalid_argument("ScaleDown before ModUp requires at least one remaining level after SlotsToCoeff");
    }

    BootstrapStcModUpPlan plan;
    plan.message_scale_log2 = message_scale_log2;
    plan.bootstrap_scale_log2 = target_pre_evalmod_scale_log2;
    plan.coeff_to_slot_log_scale = target_pre_evalmod_scale_log2;
    plan.slots_to_coeff_log_scale = std::log2(after_slots_to_coeff.scale);

    plan.scale_down.input_scale_log2 = std::log2(after_slots_to_coeff.scale);
    plan.scale_down.target_scale_log2 = target_pre_evalmod_scale_log2;
    plan.scale_down.scale_down_log2 =
        std::max(0.0, plan.scale_down.input_scale_log2 - target_pre_evalmod_scale_log2);
    plan.scale_down.consumed_levels = plan.scale_down.scale_down_log2 == 0.0 ? 0 : 1;
    plan.scale_down.expected_period_log2 =
        after_slots_to_coeff.coeff_modulus_log2 - plan.scale_down.scale_down_log2;
    plan.scale_down.expected_message_gain_log2 =
        target_pre_evalmod_scale_log2 - message_scale_log2;

    std::ostringstream note;
    note << "diagnostic STC-first ScaleDown plan; preserves message gain when target_pre_evalmod_scale_log2 "
         << "matches message_scale_log2; expected_period_log2 is a model value and must be verified by "
         << "lattice scan before EvalMod";
    plan.scale_down.note = note.str();
    return plan;
}

BootstrapStcScaleDownResult apply_stc_scale_down(SealAdapter& adapter,
                                                 const Cipher& coeff_cipher,
                                                 const BootstrapStcScaleDownPlan& plan) {
    validate_log2(plan.input_scale_log2, "plan.input_scale_log2");
    validate_log2(plan.target_scale_log2, "plan.target_scale_log2");
    validate_log2(plan.scale_down_log2, "plan.scale_down_log2");
    if (plan.scale_down_log2 < 0.0) {
        throw std::invalid_argument("plan.scale_down_log2 must be non-negative");
    }

    const auto before = adapter.info(coeff_cipher);
    BootstrapStcScaleDownResult result;
    result.chain_before = before.chain_index;
    result.scale_before_log2 = std::log2(before.scale);
    result.coeff_modulus_log2_before = before.coeff_modulus_log2;

    if (!std::isfinite(result.scale_before_log2) || before.scale <= 0.0) {
        throw std::runtime_error("ScaleDown input scale is invalid");
    }
    if (std::fabs(result.scale_before_log2 - plan.input_scale_log2) > 2.0) {
        throw std::runtime_error("ScaleDown plan input scale does not match ciphertext scale");
    }
    if (plan.scale_down_log2 == 0.0) {
        result.result = coeff_cipher;
    } else {
        if (before.chain_index == 0) {
            throw std::runtime_error("cannot apply STC ScaleDown at chain index 0");
        }
        const auto bits = adapter.coeff_modulus_bits();
        if (before.coeff_modulus_size == 0 || before.coeff_modulus_size > bits.size()) {
            throw std::runtime_error("cannot infer STC ScaleDown rescale prime");
        }
        const double plain_scale_log2 = static_cast<double>(bits[before.coeff_modulus_size - 1]);
        const double scalar = std::exp2(-plan.scale_down_log2);
        if (!std::isfinite(scalar) || scalar <= 0.0) {
            throw std::runtime_error("STC ScaleDown scalar is not representable");
        }
        result.result = adapter.mul_plain_rescale(
            coeff_cipher,
            adapter.encode_scalar_at_scale_like(scalar, std::exp2(plain_scale_log2), coeff_cipher));
    }

    const auto after = adapter.info(result.result);
    if (!std::isfinite(after.scale) || after.scale <= 0.0) {
        throw std::runtime_error("STC ScaleDown produced invalid scale");
    }
    result.chain_after = after.chain_index;
    result.scale_after_log2 = std::log2(after.scale);
    result.coeff_modulus_log2_after = after.coeff_modulus_log2;
    result.levels_consumed = before.chain_index >= after.chain_index
        ? before.chain_index - after.chain_index
        : 0;
    return result;
}

} // namespace m2424
