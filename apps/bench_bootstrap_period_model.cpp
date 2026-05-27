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

struct ProfileCase {
    const char* name{};
    m2424::CkksProfile profile{};
    double min_period_log2{};
    double max_period_log2{};
};

std::vector<double> make_input(std::size_t slots, double amplitude) {
    std::vector<double> input;
    input.reserve(slots);
    for (std::size_t i = 0; i < slots; ++i) {
        input.push_back(amplitude * std::sin(static_cast<double>(i) / 4.0));
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

void sanitize(std::string& text) {
    std::replace(text.begin(), text.end(), ',', ';');
}

std::size_t run_profile(const ProfileCase& profile_case) {
    constexpr std::size_t slots = 16;
    constexpr double tolerance = 2e-5;
    constexpr double amplitude = 1e-5;
    constexpr std::size_t level_drop = 2;
    constexpr std::size_t required_levels_for_p3 = 3;
    const auto& profile = profile_case.profile;

    auto coeff_to_slot = m2424::DiagonalLinearTransform::from_matrix(
        m2424::canonical_embedding_matrix(slots));
    auto adapter = m2424::SealAdapter::create(profile);
    adapter.keygen(coeff_to_slot.rotation_steps(), true);

    auto current = adapter.encrypt(adapter.encode(make_input(slots, amplitude)));
    for (std::size_t i = 0; i < level_drop; ++i) {
        current = adapter.mul_plain_rescale(current, adapter.encode_scalar_like(1.0, current));
    }
    const auto before_mod_raise = adapter.info(current);
    current = adapter.mod_raise_to_first(current);
    const auto after_mod_raise = adapter.info(current);
    current = coeff_to_slot.apply(adapter, current);
    const auto before_normalization = head(adapter.decode_complex(adapter.decrypt(current)), slots);
    const double before_max_abs = max_abs_value(before_normalization);

    std::size_t pass_count = 0;
    double best_period_log2 = 0.0;
    double best_plain_scale_log2 = 0.0;
    std::size_t best_chain_remaining = 0;
    double best_error = 0.0;
    double best_evalmod_ready_period_log2 = 0.0;
    double best_evalmod_ready_plain_scale_log2 = 0.0;
    std::size_t best_evalmod_ready_chain_remaining = 0;
    double best_evalmod_ready_scale_log2 = 0.0;
    double best_evalmod_ready_error = 0.0;
    std::size_t total_cases = 0;
    for (double period_log2 = profile_case.min_period_log2;
         period_log2 <= profile_case.max_period_log2;
         period_log2 += 2.0) {
        for (double plain_scale_log2 : {40.0, 50.0, 60.0, 80.0, 100.0, 120.0, 160.0, 200.0, 220.0, 240.0, 260.0}) {
            ++total_cases;
            const auto info_before = adapter.info(current);
            const double factor_log2 = -period_log2;
            const double factor = std::exp2(factor_log2);
            const auto expected = scaled(before_normalization, factor);
            const double expected_max = max_abs_value(expected);
            double actual_max = 0.0;
            double normalization_error = 0.0;
            std::size_t chunks = 0;
            std::size_t levels_consumed = 0;
            std::size_t chain_remaining = 0;
            double scale_log2_before_evalmod = 0.0;
            double coeff_modulus_log2_before_evalmod = 0.0;
            bool scalar_pass = false;
            bool evalmod_ready = false;
            bool levels_ready = false;
            bool evalmod_scale_ready = false;
            bool p3_ready = false;
            std::string exception;
            const char* status = "FAIL";

            try {
                auto normalized = m2424::apply_bootstrap_scalar_decomposed(
                    adapter, current, factor_log2, plain_scale_log2);
                chunks = normalized.chunks;
                levels_consumed = normalized.levels_consumed;
                const auto normalized_info = adapter.info(normalized.result);
                chain_remaining = normalized_info.chain_index;
                scale_log2_before_evalmod = std::log2(normalized_info.scale);
                coeff_modulus_log2_before_evalmod = normalized_info.coeff_modulus_log2;
                const auto actual = head(adapter.decode_complex(adapter.decrypt(normalized.result)), slots);
                actual_max = max_abs_value(actual);
                normalization_error = max_complex_error(expected, actual);
                scalar_pass = normalization_error <= tolerance;
                evalmod_ready = scalar_pass && expected_max <= m2424::EvalModPolynomial::approximation_bound;
                levels_ready = chain_remaining >= required_levels_for_p3;
                evalmod_scale_ready = scale_log2_before_evalmod <= 60.0;
                p3_ready = evalmod_ready && levels_ready && evalmod_scale_ready;
                if (evalmod_ready && levels_ready) {
                    if (best_evalmod_ready_period_log2 == 0.0 ||
                        scale_log2_before_evalmod < best_evalmod_ready_scale_log2 ||
                        (scale_log2_before_evalmod == best_evalmod_ready_scale_log2 &&
                         chain_remaining > best_evalmod_ready_chain_remaining)) {
                        best_evalmod_ready_period_log2 = period_log2;
                        best_evalmod_ready_plain_scale_log2 = plain_scale_log2;
                        best_evalmod_ready_chain_remaining = chain_remaining;
                        best_evalmod_ready_scale_log2 = scale_log2_before_evalmod;
                        best_evalmod_ready_error = normalization_error;
                    }
                }
                status = p3_ready ? "PASS" : "FAIL";
                if (p3_ready) {
                    ++pass_count;
                    if (best_period_log2 == 0.0 ||
                        chain_remaining > best_chain_remaining ||
                        (chain_remaining == best_chain_remaining && normalization_error < best_error)) {
                        best_period_log2 = period_log2;
                        best_plain_scale_log2 = plain_scale_log2;
                        best_chain_remaining = chain_remaining;
                        best_error = normalization_error;
                    }
                }
            } catch (const std::exception& e) {
                exception = e.what();
                sanitize(exception);
            }

            std::printf("%s,%.0f,%.0f,%zu,%.6e,%.6e,%.6e,%.6e,%.6e,%zu,%zu,%zu,%.6e,%.6e,%s,%s,%s,%s,%s,%s,%s\n",
                        profile_case.name,
                        period_log2,
                        plain_scale_log2,
                        info_before.chain_index,
                        std::log2(info_before.scale),
                        info_before.coeff_modulus_log2,
                        expected_max,
                        actual_max,
                        normalization_error,
                        chunks,
                        levels_consumed,
                        chain_remaining,
                        scale_log2_before_evalmod,
                        coeff_modulus_log2_before_evalmod,
                        scalar_pass ? "true" : "false",
                        evalmod_ready ? "true" : "false",
                        levels_ready ? "true" : "false",
                        evalmod_scale_ready ? "true" : "false",
                        p3_ready ? "true" : "false",
                        exception.c_str(),
                        status);
        }
    }

    const char* blocker = "none";
    if (pass_count == 0) {
        blocker = best_evalmod_ready_period_log2 > 0.0 ? "evalmod_scale" : "period_or_scalar";
    }

    std::printf("summary,%s,%zu,%zu,%.0f,%.0f,%zu,%.6e,%.6e,%.0f,%.0f,%zu,%.6e,%.6e,%s\n",
                profile_case.name,
                total_cases,
                pass_count,
                best_period_log2,
                best_plain_scale_log2,
                best_chain_remaining,
                best_error,
                before_max_abs,
                best_evalmod_ready_period_log2,
                best_evalmod_ready_plain_scale_log2,
                best_evalmod_ready_chain_remaining,
                best_evalmod_ready_scale_log2,
                best_evalmod_ready_error,
                blocker);
    return pass_count;
}

} // namespace

int main() {
    std::printf("profile,period_log2,plain_scale_log2,chain_before,scale_log2_before,coeff_modulus_log2_before,expected_max_abs_after_normalization,actual_max_abs_after_normalization,normalization_error,normalization_chunks,normalization_levels_consumed,chain_remaining_before_evalmod,scale_log2_before_evalmod,coeff_modulus_log2_before_evalmod,scalar_pass,evalmod_ready,levels_ready,evalmod_scale_ready,p3_ready,exception,status\n");
    std::printf("summary_header,profile,total_cases,p3_ready_cases,best_period_log2,best_plain_scale_log2,best_chain_remaining,best_error,before_max_abs,best_evalmod_ready_period_log2,best_evalmod_ready_plain_scale_log2,best_evalmod_ready_chain_remaining,best_evalmod_ready_scale_log2,best_evalmod_ready_error,blocker\n");
    std::size_t total_pass_count = 0;
    total_pass_count += run_profile({"boot_ckks", m2424::profiles::boot_ckks(), 220.0, 280.0});
    total_pass_count += run_profile({"boot_deep_ckks", m2424::profiles::boot_deep_ckks(), 680.0, 800.0});
    std::printf("summary,all_profiles,total_cases,pass_cases\n");
    std::printf("summary,all_profiles,%d,%zu\n", (31 + 61) * 11, total_pass_count);
    return 0;
}
