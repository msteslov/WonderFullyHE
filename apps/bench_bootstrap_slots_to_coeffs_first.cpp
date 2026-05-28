#include "m2424/bootstrap_dft.hpp"
#include "m2424/bootstrap_scaling.hpp"
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

m2424::ComplexVector make_input(std::size_t slots, double amplitude) {
    m2424::ComplexVector input;
    input.reserve(slots);
    for (std::size_t i = 0; i < slots; ++i) {
        input.push_back({amplitude * std::sin(static_cast<double>(i) / 4.0), 0.0});
    }
    return input;
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
    double result = 0.0;
    const std::size_t n = std::min(expected.size(), actual.size());
    for (std::size_t i = 0; i < n; ++i) {
        result = std::max(result, std::abs(expected[i] - actual[i]));
    }
    return result;
}

m2424::ComplexVector evalmod_reference(const m2424::ComplexVector& input) {
    m2424::EvalModPolynomial eval_mod;
    m2424::ComplexVector result;
    result.reserve(input.size());
    for (const auto& value : input) {
        result.push_back(eval_mod.evaluate_plain(value, m2424::EvalModDegree::P3));
    }
    return result;
}

std::vector<int> merged_steps(const m2424::FactorizedLinearTransform& a,
                              const m2424::FactorizedLinearTransform& b) {
    auto steps = a.rotation_steps();
    const auto b_steps = b.rotation_steps();
    steps.insert(steps.end(), b_steps.begin(), b_steps.end());
    std::sort(steps.begin(), steps.end());
    steps.erase(std::unique(steps.begin(), steps.end()), steps.end());
    return steps;
}

const char* bool_text(bool value) {
    return value ? "true" : "false";
}

void sanitize(std::string& text) {
    std::replace(text.begin(), text.end(), ',', ';');
}

void print_stage(const char* profile_name,
                 std::size_t slots,
                 const char* stage,
                 const m2424::CipherInfo& info,
                 double actual_max_abs,
                 double expected_max_abs,
                 double max_error,
                 bool inside_evalmod_interval,
                 const char* exception,
                 const char* status,
                 std::size_t normalization_chunks = 0,
                 std::size_t normalization_levels_consumed = 0,
                 std::size_t chain_remaining_before_evalmod = 0) {
    std::printf("%s,%zu,%s,%s,P3,%zu,%.6f,%.6f,%.6e,%.6e,%.6e,%s,%s,%s,%zu,%zu,%zu\n",
                profile_name,
                slots,
                m2424::to_string(m2424::BootstrapCircuitOrder::SlotsToCoeffsFirst),
                stage,
                info.chain_index,
                std::log2(info.scale),
                info.coeff_modulus_log2,
                actual_max_abs,
                expected_max_abs,
                max_error,
                bool_text(inside_evalmod_interval),
                exception,
                status,
                normalization_chunks,
                normalization_levels_consumed,
                chain_remaining_before_evalmod);
}

void print_cipher_stage(const char* profile_name,
                        std::size_t slots,
                        const char* stage,
                        m2424::SealAdapter& adapter,
                        const m2424::Cipher& cipher,
                        const m2424::ComplexVector& expected,
                        bool inside_evalmod_interval,
                        const char* status,
                        std::size_t normalization_chunks = 0,
                        std::size_t normalization_levels_consumed = 0,
                        std::size_t chain_remaining_before_evalmod = 0) {
    std::string exception;
    double actual_max_abs = 0.0;
    double max_error = 0.0;
    const char* final_status = status;
    try {
        const auto actual = head(adapter.decode_complex(adapter.decrypt(cipher)), slots);
        actual_max_abs = max_abs_value(actual);
        max_error = max_complex_error(expected, actual);
    } catch (const std::exception& e) {
        exception = e.what();
        sanitize(exception);
        final_status = "FAIL";
    }
    print_stage(profile_name,
                slots,
                stage,
                adapter.info(cipher),
                actual_max_abs,
                max_abs_value(expected),
                max_error,
                inside_evalmod_interval,
                exception.c_str(),
                final_status,
                normalization_chunks,
                normalization_levels_consumed,
                chain_remaining_before_evalmod);
}

bool run_case(const char* profile_name,
              const m2424::CkksProfile& profile,
              std::size_t slots,
              double normalization_plain_scale_log2) {
    constexpr double amplitude = 1e-5;
    constexpr std::size_t level_drop = 5;
    constexpr std::size_t min_chain_before_evalmod = 3;
    constexpr double tolerance = 1e-3;

    auto slot_to_coeff = m2424::FactorizedLinearTransform(m2424::make_bootstrap_dft_plan(
        slots, m2424::BootstrapDftType::HomomorphicEncode, 40.0));
    auto coeff_to_slot = m2424::FactorizedLinearTransform(m2424::make_bootstrap_dft_plan(
        slots, m2424::BootstrapDftType::HomomorphicDecode, 40.0));

    auto adapter = m2424::SealAdapter::create(profile);
    adapter.keygen(merged_steps(slot_to_coeff, coeff_to_slot), true);

    const auto input = make_input(slots, amplitude);
    auto current = adapter.encrypt(adapter.encode_complex(input));
    print_cipher_stage(profile_name, slots, "encrypt", adapter, current, input, true, "DIAG");

    for (std::size_t i = 0; i < level_drop; ++i) {
        current = adapter.mul_plain_rescale(current, adapter.encode_scalar_like(1.0, current));
    }
    print_cipher_stage(profile_name, slots, "level_drop", adapter, current, input, true, "DIAG");

    const auto coeff_expected = slot_to_coeff.apply_plain(input);
    current = slot_to_coeff.apply(adapter, current);
    print_cipher_stage(profile_name, slots, "SlotToCoeff", adapter, current, coeff_expected, false, "DIAG");

    const auto before_mod_raise = adapter.info(current);
    current = adapter.mod_raise_to_first(current);
    print_cipher_stage(profile_name, slots, "ModRaise", adapter, current, coeff_expected, false, "STRUCTURAL");

    const auto slot_expected = coeff_to_slot.apply_plain(coeff_expected);
    current = coeff_to_slot.apply(adapter, current);
    const auto after_coeff_to_slot_actual = head(adapter.decode_complex(adapter.decrypt(current)), slots);
    constexpr double prototype_period_offset_log2 = 3.0;
    const double period_log2 = before_mod_raise.coeff_modulus_log2 - prototype_period_offset_log2;
    auto normalized_expected = after_coeff_to_slot_actual;
    for (auto& value : normalized_expected) {
        value *= std::exp2(-period_log2);
    }
    bool inside_evalmod_interval =
        max_abs_value(normalized_expected) <= m2424::EvalModPolynomial::approximation_bound;
    print_cipher_stage(profile_name,
                       slots,
                       "CoeffToSlot",
                       adapter,
                       current,
                       slot_expected,
                       inside_evalmod_interval,
                       "DIAG");

    try {
        auto normalized = m2424::apply_bootstrap_scalar_decomposed(
            adapter, current, -period_log2, normalization_plain_scale_log2);
        current = normalized.result;
        const auto chain_remaining_before_evalmod = adapter.info(current).chain_index;
        print_cipher_stage(profile_name,
                           slots,
                           "Normalization",
                           adapter,
                           current,
                           normalized_expected,
                           inside_evalmod_interval,
                           inside_evalmod_interval ? "PASS" : "BLOCKED",
                           normalized.chunks,
                           normalized.levels_consumed,
                           chain_remaining_before_evalmod);
    } catch (const std::exception& e) {
        auto exception = std::string(e.what());
        sanitize(exception);
        print_stage(profile_name,
                    slots,
                    "Normalization",
                    adapter.info(current),
                    0.0,
                    max_abs_value(normalized_expected),
                    0.0,
                    inside_evalmod_interval,
                    exception.c_str(),
                    "FAIL");
        return false;
    }

    if (!inside_evalmod_interval || adapter.info(current).chain_index < min_chain_before_evalmod) {
        return false;
    }

    m2424::EvalModPolynomial eval_mod;
    const auto eval_expected = evalmod_reference(normalized_expected);
    try {
        current = eval_mod.evaluate(adapter, current, m2424::EvalModDegree::P3);
    } catch (const std::exception& e) {
        auto exception = std::string(e.what());
        sanitize(exception);
        print_stage(profile_name,
                    slots,
                    "EvalMod",
                    adapter.info(current),
                    0.0,
                    max_abs_value(eval_expected),
                    0.0,
                    inside_evalmod_interval,
                    exception.c_str(),
                    "FAIL");
        return false;
    }
    print_cipher_stage(profile_name, slots, "EvalMod", adapter, current, eval_expected, true, "PASS");

    const auto actual = head(adapter.decode_complex(adapter.decrypt(current)), slots);
    const double final_error = max_complex_error(input, actual);
    print_stage(profile_name,
                slots,
                "refresh_result",
                adapter.info(current),
                max_abs_value(actual),
                max_abs_value(input),
                final_error,
                true,
                "",
                final_error <= tolerance ? "PASS" : "FAIL");
    return final_error <= tolerance;
}

} // namespace

int main() {
    std::printf("profile,slots,circuit_order,stage,evalmod_degree,chain_index,scale_log2,coeff_modulus_log2,actual_max_abs,expected_max_abs,max_error,inside_evalmod_interval,exception,status,normalization_chunks,normalization_levels_consumed,chain_remaining_before_evalmod\n");
    bool any_pass = false;
    for (std::size_t slots : {4U, 8U, 16U}) {
        try {
            any_pass = run_case("boot_ckks", m2424::profiles::boot_ckks(), slots, 40.0) || any_pass;
        } catch (const std::exception& e) {
            auto exception = std::string(e.what());
            sanitize(exception);
            std::printf("boot_ckks,%zu,%s,case,P3,0,0,0,0,0,0,false,%s,FAIL,0,0,0\n",
                        slots,
                        m2424::to_string(m2424::BootstrapCircuitOrder::SlotsToCoeffsFirst),
                        exception.c_str());
        }
    }
    return any_pass ? 0 : 1;
}
