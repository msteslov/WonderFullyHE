#include "m2424/bootstrap_scaling.hpp"
#include "m2424/diagonal_transform.hpp"
#include "m2424/eval_mod.hpp"
#include "m2424/profiles.hpp"
#include "m2424/seal_adapter.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>
#include <string>
#include <vector>

namespace {

struct ReferenceCase {
    const char* label{};
    const char* profile_name{};
    m2424::CkksProfile profile{};
    std::size_t slots{};
    m2424::BootstrapPeriodMode period_mode{};
    double manual_period_log2{};
    double plain_scale_log2{};
};

std::vector<double> make_input(std::size_t slots, double amplitude) {
    std::vector<double> input;
    input.reserve(slots);
    for (std::size_t i = 0; i < slots; ++i) {
        input.push_back(amplitude * std::sin(static_cast<double>(i) / 4.0));
    }
    return input;
}

m2424::ComplexVector to_complex(const std::vector<double>& values) {
    m2424::ComplexVector result;
    result.reserve(values.size());
    for (double value : values) {
        result.push_back({value, 0.0});
    }
    return result;
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

m2424::ComplexVector evalmod_reference(const m2424::EvalModPolynomial& eval_mod,
                                       const m2424::ComplexVector& values,
                                       m2424::EvalModDegree degree) {
    m2424::ComplexVector result;
    result.reserve(values.size());
    for (const auto& value : values) {
        result.push_back(eval_mod.evaluate_plain(value, degree));
    }
    return result;
}

double amplitude_normalization_factor(const m2424::ComplexVector& values) {
    constexpr double target = m2424::EvalModPolynomial::approximation_bound * 0.5;
    const double max_abs = max_abs_value(values);
    return max_abs > target ? target / max_abs : 1.0;
}

void sanitize(std::string& text) {
    std::replace(text.begin(), text.end(), ',', ';');
}

const char* bool_text(bool value) {
    return value ? "true" : "false";
}

void print_stage(const ReferenceCase& config,
                 const char* stage,
                 const char* domain,
                 const m2424::CipherInfo& info,
                 double actual_max_abs,
                 double expected_max_abs,
                 double max_error,
                 bool inside_evalmod_interval,
                 std::size_t normalization_chunks,
                 std::size_t normalization_levels_consumed,
                 std::size_t scale_squash_levels_consumed,
                 std::size_t denormalization_chunks,
                 std::size_t denormalization_levels_consumed,
                 std::size_t chain_remaining_before_evalmod,
                 std::size_t chain_remaining_after_evalmod,
                 const std::string& blocker,
                 const std::string& exception,
                 const char* status) {
    std::printf("%s,%s,%zu,%s,%s,%.0f,%.0f,%s,%zu,%.6f,%.6f,%.6e,%.6e,%.6e,%s,%zu,%zu,%zu,%zu,%zu,%zu,%zu,%s,%s,%s\n",
                config.label,
                config.profile_name,
                config.slots,
                stage,
                m2424::to_string(config.period_mode),
                config.period_mode == m2424::BootstrapPeriodMode::ManualPowerOfTwo
                    ? config.manual_period_log2
                    : 0.0,
                config.plain_scale_log2,
                domain,
                info.chain_index,
                std::log2(info.scale),
                info.coeff_modulus_log2,
                actual_max_abs,
                expected_max_abs,
                max_error,
                bool_text(inside_evalmod_interval),
                normalization_chunks,
                normalization_levels_consumed,
                scale_squash_levels_consumed,
                denormalization_chunks,
                denormalization_levels_consumed,
                chain_remaining_before_evalmod,
                chain_remaining_after_evalmod,
                blocker.c_str(),
                exception.c_str(),
                status);
}

void print_cipher_stage(const ReferenceCase& config,
                        const char* stage,
                        const char* domain,
                        m2424::SealAdapter& adapter,
                        const m2424::Cipher& cipher,
                        const m2424::ComplexVector& expected,
                        bool inside_evalmod_interval,
                        std::size_t normalization_chunks,
                        std::size_t normalization_levels_consumed,
                        std::size_t scale_squash_levels_consumed,
                        std::size_t denormalization_chunks,
                        std::size_t denormalization_levels_consumed,
                        std::size_t chain_remaining_before_evalmod,
                        std::size_t chain_remaining_after_evalmod,
                        const char* blocker,
                        const char* status) {
    std::string exception;
    double actual_max_abs = 0.0;
    double max_error = 0.0;
    try {
        const auto actual = head(adapter.decode_complex(adapter.decrypt(cipher)), config.slots);
        actual_max_abs = max_abs_value(actual);
        max_error = max_complex_error(expected, actual);
    } catch (const std::exception& e) {
        exception = e.what();
        sanitize(exception);
        status = "FAIL";
    }

    print_stage(config,
                stage,
                domain,
                adapter.info(cipher),
                actual_max_abs,
                max_abs_value(expected),
                max_error,
                inside_evalmod_interval,
                normalization_chunks,
                normalization_levels_consumed,
                scale_squash_levels_consumed,
                denormalization_chunks,
                denormalization_levels_consumed,
                chain_remaining_before_evalmod,
                chain_remaining_after_evalmod,
                blocker,
                exception,
                status);
}

std::vector<int> rotation_steps_for(const m2424::DiagonalLinearTransform& coeff_to_slot,
                                    const m2424::DiagonalLinearTransform& slot_to_coeff) {
    auto steps = coeff_to_slot.rotation_steps();
    const auto inverse_steps = slot_to_coeff.rotation_steps();
    steps.insert(steps.end(), inverse_steps.begin(), inverse_steps.end());
    std::sort(steps.begin(), steps.end());
    steps.erase(std::unique(steps.begin(), steps.end()), steps.end());
    return steps;
}

void print_blocked_stage(const ReferenceCase& config,
                         const char* stage,
                         const char* domain,
                         const m2424::CipherInfo& info,
                         const m2424::ComplexVector& expected,
                         bool inside_evalmod_interval,
                         std::size_t normalization_chunks,
                         std::size_t normalization_levels_consumed,
                         std::size_t scale_squash_levels_consumed,
                         std::size_t denormalization_chunks,
                         std::size_t denormalization_levels_consumed,
                         std::size_t chain_remaining_before_evalmod,
                         std::size_t chain_remaining_after_evalmod,
                         const char* blocker,
                         const std::exception& e) {
    auto exception = std::string(e.what());
    sanitize(exception);
    print_stage(config,
                stage,
                domain,
                info,
                0.0,
                max_abs_value(expected),
                0.0,
                inside_evalmod_interval,
                normalization_chunks,
                normalization_levels_consumed,
                scale_squash_levels_consumed,
                denormalization_chunks,
                denormalization_levels_consumed,
                chain_remaining_before_evalmod,
                chain_remaining_after_evalmod,
                blocker,
                exception,
                "BLOCKED");
}

void run_case(const ReferenceCase& config) {
    constexpr double amplitude = 1e-5;
    constexpr std::size_t level_drop = 2;
    constexpr double target_evalmod_scale_log2 = 60.0;
    constexpr std::size_t min_chain_before_evalmod = 3;
    constexpr m2424::EvalModDegree degree = m2424::EvalModDegree::P3;

    auto matrix = m2424::canonical_embedding_matrix(config.slots);
    auto coeff_to_slot = m2424::DiagonalLinearTransform::from_matrix(matrix);
    auto slot_to_coeff = m2424::DiagonalLinearTransform::from_matrix(m2424::invert_matrix(matrix));
    auto adapter = m2424::SealAdapter::create(config.profile);
    adapter.keygen(rotation_steps_for(coeff_to_slot, slot_to_coeff), true);

    const auto input = make_input(config.slots, amplitude);
    const auto input_expected = to_complex(input);
    auto current = adapter.encrypt(adapter.encode(input));

    print_cipher_stage(config, "encrypt", "coefficient", adapter, current, input_expected, false,
                       0, 0, 0, 0, 0, 0, 0, "none", "DIAG");

    for (std::size_t i = 0; i < level_drop; ++i) {
        current = adapter.mul_plain_rescale(current, adapter.encode_scalar_like(1.0, current));
    }
    print_cipher_stage(config, "level_drop", "coefficient", adapter, current, input_expected, false,
                       0, 0, 0, 0, 0, 0, 0, "none", "DIAG");

    const auto before_mod_raise = adapter.info(current);
    current = adapter.mod_raise_to_first(current);
    const auto after_mod_raise = adapter.info(current);
    print_cipher_stage(config, "ModRaise", "coefficient", adapter, current, input_expected, false,
                       0, 0, 0, 0, 0, 0, 0, "structural_only", "STRUCTURAL");

    const auto coeff_to_slot_reference = coeff_to_slot.apply_plain(input_expected);
    current = coeff_to_slot.apply(adapter, current);
    print_cipher_stage(config, "CoeffToSlot", "slot", adapter, current, coeff_to_slot_reference, false,
                       0, 0, 0, 0, 0, 0, 0, "diagnostic_after_structural_modraise", "DIAG");

    const auto coeff_to_slot_actual_base = head(adapter.decode_complex(adapter.decrypt(current)), config.slots);
    const double period_log2 = m2424::bootstrap_period_log2(config.period_mode,
                                                            config.manual_period_log2,
                                                            config.profile.coeff_modulus_bits,
                                                            before_mod_raise,
                                                            after_mod_raise);
    const auto scaling = m2424::make_bootstrap_scaling_factors(
        amplitude_normalization_factor(coeff_to_slot_actual_base),
        period_log2,
        config.plain_scale_log2);

    std::size_t normalization_chunks = 0;
    std::size_t normalization_levels_consumed = 0;
    std::size_t scale_squash_levels_consumed = 0;
    std::size_t denormalization_chunks = 0;
    std::size_t denormalization_levels_consumed = 0;
    std::size_t chain_remaining_before_evalmod = 0;
    std::size_t chain_remaining_after_evalmod = 0;

    auto normalized_expected = scaled(coeff_to_slot_actual_base, scaling.factor);
    const bool inside_evalmod_interval =
        max_abs_value(normalized_expected) <= m2424::EvalModPolynomial::approximation_bound;

    try {
        auto normalized = m2424::apply_bootstrap_scalar_decomposed(
            adapter, current, scaling.normalization_factor_log2, config.plain_scale_log2);
        current = normalized.result;
        normalization_chunks = normalized.chunks;
        normalization_levels_consumed = normalized.levels_consumed;
        chain_remaining_before_evalmod = adapter.info(current).chain_index;
        print_cipher_stage(config, "Normalization", "slot", adapter, current, normalized_expected,
                           inside_evalmod_interval, normalization_chunks, normalization_levels_consumed,
                           scale_squash_levels_consumed, 0, 0, chain_remaining_before_evalmod, 0,
                           inside_evalmod_interval ? "none" : "period_model_blocked", "DIAG");
    } catch (const std::exception& e) {
        print_blocked_stage(config, "Normalization", "slot", adapter.info(current), normalized_expected,
                            inside_evalmod_interval, 0, 0, 0, 0, 0, 0, 0,
                            "scale_strategy_blocked", e);
        return;
    }

    try {
        auto squashed = m2424::squash_bootstrap_scale(
            adapter, current, target_evalmod_scale_log2, min_chain_before_evalmod);
        current = squashed.result;
        scale_squash_levels_consumed = squashed.levels_consumed;
        chain_remaining_before_evalmod = adapter.info(current).chain_index;
        print_cipher_stage(config, "ScaleSquash", "slot", adapter, current, normalized_expected,
                           inside_evalmod_interval, normalization_chunks, normalization_levels_consumed,
                           scale_squash_levels_consumed, 0, 0, chain_remaining_before_evalmod, 0,
                           inside_evalmod_interval ? "none" : "period_model_blocked", "DIAG");
    } catch (const std::exception& e) {
        print_blocked_stage(config, "ScaleSquash", "slot", adapter.info(current), normalized_expected,
                            inside_evalmod_interval, normalization_chunks, normalization_levels_consumed,
                            scale_squash_levels_consumed, 0, 0, chain_remaining_before_evalmod, 0,
                            "scale_strategy_blocked", e);
        return;
    }

    if (!inside_evalmod_interval) {
        print_stage(config, "EvalMod", "slot", adapter.info(current), 0.0, max_abs_value(normalized_expected),
                    0.0, false, normalization_chunks, normalization_levels_consumed,
                    scale_squash_levels_consumed, 0, 0, chain_remaining_before_evalmod, 0,
                    "period_model_blocked", "input outside EvalMod interval", "BLOCKED");
        return;
    }

    m2424::EvalModPolynomial eval_mod;
    const auto evalmod_expected = evalmod_reference(eval_mod, normalized_expected, degree);
    try {
        current = eval_mod.evaluate(adapter, current, degree);
        chain_remaining_after_evalmod = adapter.info(current).chain_index;
        print_cipher_stage(config, "EvalMod", "slot", adapter, current, evalmod_expected, true,
                           normalization_chunks, normalization_levels_consumed,
                           scale_squash_levels_consumed, 0, 0,
                           chain_remaining_before_evalmod, chain_remaining_after_evalmod,
                           "none", "DIAG");
    } catch (const std::exception& e) {
        print_blocked_stage(config, "EvalMod", "slot", adapter.info(current), evalmod_expected, true,
                            normalization_chunks, normalization_levels_consumed,
                            scale_squash_levels_consumed, 0, 0,
                            chain_remaining_before_evalmod, 0,
                            "evalmod_capacity_blocked", e);
        return;
    }

    const auto denormalized_expected = scaled(evalmod_expected, 1.0 / scaling.factor);
    try {
        auto denormalized = m2424::apply_bootstrap_scalar_decomposed(
            adapter, current, -scaling.normalization_factor_log2, config.plain_scale_log2);
        current = denormalized.result;
        denormalization_chunks = denormalized.chunks;
        denormalization_levels_consumed = denormalized.levels_consumed;
        print_cipher_stage(config, "Denormalization", "slot", adapter, current, denormalized_expected, true,
                           normalization_chunks, normalization_levels_consumed,
                           scale_squash_levels_consumed,
                           denormalization_chunks, denormalization_levels_consumed,
                           chain_remaining_before_evalmod, chain_remaining_after_evalmod,
                           "none", "DIAG");
    } catch (const std::exception& e) {
        print_blocked_stage(config, "Denormalization", "slot", adapter.info(current), denormalized_expected,
                            true, normalization_chunks, normalization_levels_consumed,
                            scale_squash_levels_consumed, 0, 0,
                            chain_remaining_before_evalmod, chain_remaining_after_evalmod,
                            "denormalization_blocked", e);
        return;
    }

    const auto slot_to_coeff_expected = slot_to_coeff.apply_plain(denormalized_expected);
    try {
        current = slot_to_coeff.apply(adapter, current);
        print_cipher_stage(config, "SlotToCoeff", "coefficient", adapter, current, slot_to_coeff_expected,
                           true, normalization_chunks, normalization_levels_consumed,
                           scale_squash_levels_consumed,
                           denormalization_chunks, denormalization_levels_consumed,
                           chain_remaining_before_evalmod, chain_remaining_after_evalmod,
                           "ready_for_evalmod_p3_path", "DIAG");
    } catch (const std::exception& e) {
        print_blocked_stage(config, "SlotToCoeff", "coefficient", adapter.info(current),
                            slot_to_coeff_expected, true, normalization_chunks,
                            normalization_levels_consumed, scale_squash_levels_consumed,
                            denormalization_chunks, denormalization_levels_consumed,
                            chain_remaining_before_evalmod, chain_remaining_after_evalmod,
                            "slot_to_coeff_blocked", e);
    }
}

} // namespace

int main() {
    std::printf("case,profile,slots,stage,period_mode,manual_period_log2,plain_scale_log2,domain,chain_index,scale_log2,coeff_modulus_log2,actual_max_abs,expected_max_abs,max_error,inside_evalmod_interval,normalization_chunks,normalization_levels_consumed,scale_squash_levels_consumed,denormalization_chunks,denormalization_levels_consumed,chain_remaining_before_evalmod,chain_remaining_after_evalmod,blocker,exception,status\n");

    for (std::size_t slots : {std::size_t{4}, std::size_t{8}, std::size_t{16}}) {
        run_case({"diagnostic_no_period", "boot_ckks", m2424::profiles::boot_ckks(), slots,
                  m2424::BootstrapPeriodMode::NoBootstrapPeriod, 0.0, 50.0});
        run_case({"manual_220_scalar_mechanics", "boot_ckks", m2424::profiles::boot_ckks(), slots,
                  m2424::BootstrapPeriodMode::ManualPowerOfTwo, 220.0, 240.0});
        run_case({"manual_256_period_probe", "boot_ckks", m2424::profiles::boot_ckks(), slots,
                  m2424::BootstrapPeriodMode::ManualPowerOfTwo, 256.0, 160.0});
    }

    run_case({"deep_manual_736_p3_probe", "boot_deep_ckks", m2424::profiles::boot_deep_ckks(), 16,
              m2424::BootstrapPeriodMode::ManualPowerOfTwo, 736.0, 50.0});
    return 0;
}
