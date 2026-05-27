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

struct TraceCase {
    const char* label{};
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

m2424::ComplexVector evaluate_plain(const m2424::EvalModPolynomial& eval_mod,
                                    const m2424::ComplexVector& values,
                                    m2424::EvalModDegree degree) {
    m2424::ComplexVector result;
    result.reserve(values.size());
    for (const auto& value : values) {
        result.push_back(eval_mod.evaluate_plain(value, degree));
    }
    return result;
}

double normalization_factor_for(const m2424::ComplexVector& values) {
    constexpr double evalmod_target = m2424::EvalModPolynomial::approximation_bound * 0.5;
    const double max_abs = max_abs_value(values);
    return max_abs > evalmod_target ? evalmod_target / max_abs : 1.0;
}

void sanitize(std::string& text) {
    std::replace(text.begin(), text.end(), ',', ';');
}

void print_stage(const char* label,
                 const char* stage,
                 const TraceCase& config,
                 const m2424::CipherInfo& info,
                 double max_abs,
                 double expected_max_abs,
                 double max_error,
                 bool inside_evalmod_interval,
                 std::size_t normalization_chunks,
                 std::size_t normalization_levels_consumed,
                 std::size_t denormalization_chunks,
                 std::size_t denormalization_levels_consumed,
                 std::size_t chain_remaining_before_evalmod,
                 std::size_t chain_remaining_after_evalmod,
                 const std::string& exception,
                 const char* status) {
    std::printf("%s,%s,%s,%.0f,%.0f,%zu,%.6e,%.6e,%.6e,%.6e,%.6e,%s,%zu,%zu,%zu,%zu,%zu,%zu,%s,%s\n",
                label,
                stage,
                m2424::to_string(config.period_mode),
                config.period_mode == m2424::BootstrapPeriodMode::ManualPowerOfTwo
                    ? config.manual_period_log2
                    : 0.0,
                config.plain_scale_log2,
                info.chain_index,
                std::log2(info.scale),
                info.coeff_modulus_log2,
                max_abs,
                expected_max_abs,
                max_error,
                inside_evalmod_interval ? "true" : "false",
                normalization_chunks,
                normalization_levels_consumed,
                denormalization_chunks,
                denormalization_levels_consumed,
                chain_remaining_before_evalmod,
                chain_remaining_after_evalmod,
                exception.c_str(),
                status);
}

void print_cipher_stage(const char* label,
                        const char* stage,
                        const TraceCase& config,
                        m2424::SealAdapter& adapter,
                        const m2424::Cipher& cipher,
                        const m2424::ComplexVector& expected,
                        std::size_t slots,
                        bool inside_evalmod_interval,
                        std::size_t normalization_chunks,
                        std::size_t normalization_levels_consumed,
                        std::size_t denormalization_chunks,
                        std::size_t denormalization_levels_consumed,
                        std::size_t chain_remaining_before_evalmod,
                        std::size_t chain_remaining_after_evalmod,
                        const char* status) {
    auto exception = std::string{};
    double max_abs = 0.0;
    double max_error = 0.0;
    try {
        const auto actual = head(adapter.decode_complex(adapter.decrypt(cipher)), slots);
        max_abs = max_abs_value(actual);
        max_error = max_complex_error(expected, actual);
    } catch (const std::exception& e) {
        exception = e.what();
        sanitize(exception);
        status = "FAIL";
    }
    print_stage(label,
                stage,
                config,
                adapter.info(cipher),
                max_abs,
                max_abs_value(expected),
                max_error,
                inside_evalmod_interval,
                normalization_chunks,
                normalization_levels_consumed,
                denormalization_chunks,
                denormalization_levels_consumed,
                chain_remaining_before_evalmod,
                chain_remaining_after_evalmod,
                exception,
                status);
}

void run_case(const TraceCase& config) {
    constexpr std::size_t slots = 16;
    constexpr double amplitude = 1e-5;
    constexpr std::size_t level_drop = 2;
    const auto profile = m2424::profiles::boot_ckks();

    auto coeff_to_slot = m2424::DiagonalLinearTransform::from_matrix(
        m2424::canonical_embedding_matrix(slots));
    auto slot_to_coeff = m2424::DiagonalLinearTransform::from_matrix(
        m2424::invert_matrix(m2424::canonical_embedding_matrix(slots)));
    auto rotation_steps = coeff_to_slot.rotation_steps();
    const auto slot_to_coeff_steps = slot_to_coeff.rotation_steps();
    rotation_steps.insert(rotation_steps.end(), slot_to_coeff_steps.begin(), slot_to_coeff_steps.end());
    std::sort(rotation_steps.begin(), rotation_steps.end());
    rotation_steps.erase(std::unique(rotation_steps.begin(), rotation_steps.end()), rotation_steps.end());

    auto adapter = m2424::SealAdapter::create(profile);
    adapter.keygen(rotation_steps, true);
    const auto input = make_input(slots, amplitude);
    const auto input_expected = to_complex(input);

    auto current = adapter.encrypt(adapter.encode(input));
    print_cipher_stage(config.label, "encrypt", config, adapter, current, input_expected, slots,
                       false, 0, 0, 0, 0, 0, 0, "DIAG");

    for (std::size_t i = 0; i < level_drop; ++i) {
        current = adapter.mul_plain_rescale(current, adapter.encode_scalar_like(1.0, current));
    }
    print_cipher_stage(config.label, "level_drop", config, adapter, current, input_expected, slots,
                       false, 0, 0, 0, 0, 0, 0, "DIAG");

    const auto before_mod_raise = adapter.info(current);
    current = adapter.mod_raise_to_first(current);
    const auto after_mod_raise = adapter.info(current);
    print_cipher_stage(config.label, "ModRaise", config, adapter, current, input_expected, slots,
                       false, 0, 0, 0, 0, 0, 0, "STRUCTURAL");

    auto coeff_expected = coeff_to_slot.apply_plain(input_expected);
    current = coeff_to_slot.apply(adapter, current);
    print_cipher_stage(config.label, "CoeffToSlot", config, adapter, current, coeff_expected, slots,
                       false, 0, 0, 0, 0, 0, 0, "DIAG");
    auto coeff_scalar_base = head(adapter.decode_complex(adapter.decrypt(current)), slots);

    const double period_log2 = m2424::bootstrap_period_log2(
        config.period_mode,
        config.manual_period_log2,
        profile.coeff_modulus_bits,
        before_mod_raise,
        after_mod_raise);
    const auto scaling = m2424::make_bootstrap_scaling_factors(
        normalization_factor_for(coeff_scalar_base), period_log2, config.plain_scale_log2);

    std::size_t normalization_chunks = 0;
    std::size_t normalization_levels_consumed = 0;
    std::size_t denormalization_chunks = 0;
    std::size_t denormalization_levels_consumed = 0;
    std::size_t chain_remaining_before_evalmod = 0;
    std::size_t chain_remaining_after_evalmod = 0;

    coeff_expected = scaled(coeff_scalar_base, scaling.factor);
    const bool inside_evalmod_interval =
        max_abs_value(coeff_expected) <= m2424::EvalModPolynomial::approximation_bound;
    try {
        auto normalized = m2424::apply_bootstrap_scalar_decomposed(
            adapter, current, scaling.normalization_factor_log2, config.plain_scale_log2);
        current = normalized.result;
        normalization_chunks = normalized.chunks;
        normalization_levels_consumed = normalized.levels_consumed;
        chain_remaining_before_evalmod = adapter.info(current).chain_index;
        print_cipher_stage(config.label, "Normalization", config, adapter, current, coeff_expected, slots,
                           inside_evalmod_interval, normalization_chunks, normalization_levels_consumed,
                           0, 0, chain_remaining_before_evalmod, 0, "DIAG");
    } catch (const std::exception& e) {
        auto exception = std::string(e.what());
        sanitize(exception);
        print_stage(config.label, "Normalization", config, adapter.info(current), 0.0,
                    max_abs_value(coeff_expected), 0.0, inside_evalmod_interval, 0, 0, 0, 0,
                    0, 0, exception, "FAIL");
        return;
    }

    if (!inside_evalmod_interval) {
        print_stage(config.label, "EvalMod", config, adapter.info(current), 0.0,
                    max_abs_value(coeff_expected), 0.0, false,
                    normalization_chunks, normalization_levels_consumed, 0, 0,
                    chain_remaining_before_evalmod, 0,
                    "input outside EvalMod interval", "BLOCKED");
        return;
    }

    m2424::EvalModPolynomial eval_mod;
    auto eval_expected = evaluate_plain(eval_mod, coeff_expected, m2424::EvalModDegree::P3);
    try {
        current = eval_mod.evaluate(adapter, current, m2424::EvalModDegree::P3);
        chain_remaining_after_evalmod = adapter.info(current).chain_index;
        print_cipher_stage(config.label, "EvalMod", config, adapter, current, eval_expected, slots,
                           true, normalization_chunks, normalization_levels_consumed, 0, 0,
                           chain_remaining_before_evalmod, chain_remaining_after_evalmod, "DIAG");
    } catch (const std::exception& e) {
        auto exception = std::string(e.what());
        sanitize(exception);
        print_stage(config.label, "EvalMod", config, adapter.info(current), 0.0,
                    max_abs_value(eval_expected), 0.0, true,
                    normalization_chunks, normalization_levels_consumed, 0, 0,
                    chain_remaining_before_evalmod, 0, exception, "FAIL");
        return;
    }

    eval_expected = scaled(eval_expected, 1.0 / scaling.factor);
    try {
        auto denormalized = m2424::apply_bootstrap_scalar_decomposed(
            adapter, current, -scaling.normalization_factor_log2, config.plain_scale_log2);
        current = denormalized.result;
        denormalization_chunks = denormalized.chunks;
        denormalization_levels_consumed = denormalized.levels_consumed;
        print_cipher_stage(config.label, "Denormalization", config, adapter, current, eval_expected, slots,
                           true, normalization_chunks, normalization_levels_consumed,
                           denormalization_chunks, denormalization_levels_consumed,
                           chain_remaining_before_evalmod, chain_remaining_after_evalmod, "DIAG");
    } catch (const std::exception& e) {
        auto exception = std::string(e.what());
        sanitize(exception);
        print_stage(config.label, "Denormalization", config, adapter.info(current), 0.0,
                    max_abs_value(eval_expected), 0.0, true,
                    normalization_chunks, normalization_levels_consumed, 0, 0,
                    chain_remaining_before_evalmod, chain_remaining_after_evalmod, exception, "FAIL");
        return;
    }

    auto refreshed_expected = slot_to_coeff.apply_plain(eval_expected);
    try {
        current = slot_to_coeff.apply(adapter, current);
        print_cipher_stage(config.label, "SlotToCoeff", config, adapter, current, refreshed_expected, slots,
                           true, normalization_chunks, normalization_levels_consumed,
                           denormalization_chunks, denormalization_levels_consumed,
                           chain_remaining_before_evalmod, chain_remaining_after_evalmod, "DIAG");
    } catch (const std::exception& e) {
        auto exception = std::string(e.what());
        sanitize(exception);
        print_stage(config.label, "SlotToCoeff", config, adapter.info(current), 0.0,
                    max_abs_value(refreshed_expected), 0.0, true,
                    normalization_chunks, normalization_levels_consumed,
                    denormalization_chunks, denormalization_levels_consumed,
                    chain_remaining_before_evalmod, chain_remaining_after_evalmod, exception, "FAIL");
    }
}

} // namespace

int main() {
    std::printf("case,stage,period_mode,manual_period_log2,plain_scale_log2,chain_index,scale_log2,coeff_modulus_log2,max_abs,expected_max_abs,max_error,inside_evalmod_interval,normalization_chunks,normalization_levels_consumed,denormalization_chunks,denormalization_levels_consumed,chain_remaining_before_evalmod,chain_remaining_after_evalmod,exception,status\n");
    run_case({"diagnostic_no_period", m2424::BootstrapPeriodMode::NoBootstrapPeriod, 0.0, 50.0});
    run_case({"manual_220_scalar_mechanics", m2424::BootstrapPeriodMode::ManualPowerOfTwo, 220.0, 240.0});
    run_case({"manual_256_evalmod_ready", m2424::BootstrapPeriodMode::ManualPowerOfTwo, 256.0, 160.0});
    return 0;
}
