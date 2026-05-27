#include "m2424/bootstrap_prototype.hpp"
#include "m2424/profiles.hpp"
#include "m2424/seal_adapter.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::vector<double> make_input(const std::string& kind, std::size_t slots, double amplitude) {
    std::vector<double> input;
    input.reserve(slots);
    for (std::size_t i = 0; i < slots; ++i) {
        const double x = static_cast<double>(i);
        if (kind == "sine") {
            input.push_back(amplitude * std::sin(x / 4.0));
        } else if (kind == "alternating") {
            input.push_back((i % 2 == 0 ? amplitude : -amplitude));
        } else if (kind == "ramp") {
            const double centered = (x / static_cast<double>(slots - 1)) - 0.5;
            input.push_back(2.0 * amplitude * centered);
        } else {
            throw std::invalid_argument("unknown input kind");
        }
    }
    return input;
}

std::vector<m2424::Complex> to_complex(const std::vector<double>& input) {
    std::vector<m2424::Complex> result;
    result.reserve(input.size());
    for (double value : input) {
        result.push_back({value, 0.0});
    }
    return result;
}

double max_error(const std::vector<double>& expected, const std::vector<double>& actual) {
    const std::size_t n = std::min(expected.size(), actual.size());
    double result = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        result = std::max(result, std::fabs(expected[i] - actual[i]));
    }
    return result;
}

const m2424::BootstrapPrototypeStage* find_stage(const m2424::BootstrapPrototypeReport& report,
                                                 const char* name) {
    for (const auto& stage : report.stages) {
        if (stage.name == name) {
            return &stage;
        }
    }
    return nullptr;
}

double stage_error(const m2424::BootstrapPrototypeReport& report, const char* name) {
    const auto* stage = find_stage(report, name);
    return stage ? stage->max_abs_error : 0.0;
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

m2424::Cipher apply_scalar_decomposed(m2424::SealAdapter& adapter,
                                      const m2424::Cipher& input,
                                      double factor_log2,
                                      double plain_scale_log2) {
    constexpr double scale_capacity_margin_log2 = 2.0;
    if (factor_log2 + plain_scale_log2 >= 0.0) {
        return adapter.mul_plain_rescale(
            input,
            adapter.encode_scalar_at_scale_like(std::exp2(factor_log2), std::exp2(plain_scale_log2), input));
    }
    if (factor_log2 >= 0.0) {
        throw std::runtime_error("positive bootstrap scalar decomposition is not supported");
    }

    auto current = input;
    double remaining_abs_log2 = -factor_log2;
    while (remaining_abs_log2 > 1e-9) {
        const auto info = adapter.info(current);
        if (info.chain_index == 0) {
            throw std::runtime_error("not enough levels for bootstrap scalar decomposition");
        }
        const double current_scale_log2 = std::log2(info.scale);
        const double capacity_log2 =
            info.coeff_modulus_log2 - current_scale_log2 - scale_capacity_margin_log2;
        const double chunk_plain_scale_log2 = std::min(plain_scale_log2, capacity_log2);
        if (chunk_plain_scale_log2 <= 0.0) {
            throw std::runtime_error("not enough scale capacity for bootstrap scalar decomposition");
        }
        const double chunk_abs_log2 = std::min(remaining_abs_log2, chunk_plain_scale_log2);
        current = adapter.mul_plain_rescale(
            current,
            adapter.encode_scalar_at_scale_like(
                std::exp2(-chunk_abs_log2), std::exp2(chunk_plain_scale_log2), current));
        remaining_abs_log2 -= chunk_abs_log2;
    }
    return current;
}

double normalization_factor_for(const m2424::ComplexVector& values) {
    constexpr double evalmod_target = 0.0009765625 * 0.5;
    const double max_abs = max_abs_value(values);
    if (max_abs == 0.0) {
        return 1.0;
    }
    return max_abs > evalmod_target ? evalmod_target / max_abs : 1.0;
}

template <typename PeriodCases>
bool scaling_gate_has_pass(m2424::SealAdapter& adapter,
                           const m2424::CkksProfile& profile,
                           std::size_t slots,
                           double tolerance,
                           const std::vector<double>& amplitudes,
                           const std::vector<std::size_t>& level_drops,
                           const PeriodCases& period_cases,
                           const std::vector<double>& plain_scale_log2_values,
                           std::string& best_mode,
                           m2424::BootstrapPeriodMode& best_period_mode,
                           double& best_manual_period_log2,
                           double& best_plain_scale_log2) {
    auto coeff_to_slot = m2424::DiagonalLinearTransform::from_matrix(
        m2424::canonical_embedding_matrix(slots));

    for (double amplitude : amplitudes) {
        for (std::size_t level_drop : level_drops) {
            auto current = adapter.encrypt(adapter.encode(make_input("sine", slots, amplitude)));
            for (std::size_t i = 0; i < level_drop; ++i) {
                current = adapter.mul_plain_rescale(current, adapter.encode_scalar_like(1.0, current));
            }
            const auto before_mod_raise = adapter.info(current);
            current = adapter.mod_raise_to_first(current);
            const auto after_mod_raise = adapter.info(current);
            current = coeff_to_slot.apply(adapter, current);

            m2424::ComplexVector before_normalization;
            try {
                before_normalization = head(adapter.decode_complex(adapter.decrypt(current)), slots);
            } catch (...) {
                continue;
            }
            const double amplitude_factor = normalization_factor_for(before_normalization);
            const auto before_normalization_info = adapter.info(current);

            for (const auto& period_case : period_cases) {
                const double period_log2 = m2424::bootstrap_period_log2(
                    period_case.mode,
                    period_case.manual_period_log2,
                    profile.coeff_modulus_bits,
                    before_mod_raise,
                    after_mod_raise);
                for (double plain_scale_log2 : plain_scale_log2_values) {
                    const auto scaling = m2424::make_bootstrap_scaling_factors(
                        amplitude_factor, period_log2, plain_scale_log2);
                    if (!scaling.representable) {
                        continue;
                    }
                    try {
                        const auto expected = scaled(before_normalization, scaling.factor);
                        auto normalized = apply_scalar_decomposed(
                            adapter, current, scaling.normalization_factor_log2, plain_scale_log2);
                        const auto actual = head(adapter.decode_complex(adapter.decrypt(normalized)), slots);
                        const double normalization_error = max_complex_error(expected, actual);
                        if (normalization_error <= tolerance) {
                            best_period_mode = period_case.mode;
                            best_manual_period_log2 = period_case.manual_period_log2;
                            best_plain_scale_log2 = plain_scale_log2;
                            best_mode = std::string(m2424::to_string(period_case.mode)) +
                                "/manual_period_log2=" + std::to_string(period_case.manual_period_log2) +
                                "/plain_scale_log2=" + std::to_string(plain_scale_log2) +
                                "/chain_after=" + std::to_string(before_normalization_info.chain_index);
                            return true;
                        }
                    } catch (...) {
                    }
                }
            }
        }
    }
    return false;
}

} // namespace

int main() {
    constexpr std::size_t slots = 16;
    constexpr double tolerance = 2e-5;
    constexpr std::size_t post_depth = 2;

    const std::vector<double> amplitudes{1e-5, 1e-4, 5e-4, 1e-3, 1e-2, 1e-1};
    const std::vector<std::string> kinds{"sine", "alternating", "ramp"};
    const std::vector<std::size_t> level_drops{1, 2};
    const std::vector<m2424::BootstrapDenormalizationPosition> denorm_positions{
        m2424::BootstrapDenormalizationPosition::BeforeSlotToCoeff,
        m2424::BootstrapDenormalizationPosition::AfterSlotToCoeff
    };
    const std::vector<m2424::EvalModDegree> evalmod_degrees{
        m2424::EvalModDegree::P7,
        m2424::EvalModDegree::P5,
        m2424::EvalModDegree::P3
    };
    struct PeriodCase {
        m2424::BootstrapPeriodMode mode{};
        double manual_period_log2{};
    };
    std::vector<PeriodCase> period_cases{
        {m2424::BootstrapPeriodMode::NoBootstrapPeriod, 0.0},
        {m2424::BootstrapPeriodMode::SourceCoeffModulus, 0.0},
        {m2424::BootstrapPeriodMode::TotalCoeffModulus, 0.0},
        {m2424::BootstrapPeriodMode::LastPrime, 0.0},
        {m2424::BootstrapPeriodMode::DroppedPrimeProduct, 0.0}
    };
    for (double manual : {40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 120.0, 140.0, 220.0, 258.0, 260.0, 300.0}) {
        period_cases.push_back({m2424::BootstrapPeriodMode::ManualPowerOfTwo, manual});
    }
    std::vector<double> plain_scale_log2_values{
        40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 120.0,
        140.0, 160.0, 180.0, 200.0, 240.0, 280.0, 320.0, 400.0, 600.0
    };

    const auto profile = m2424::profiles::boot_ckks();
    auto rotation_steps = m2424::BootstrapPrototype::required_rotation_steps(slots);
    auto adapter = m2424::SealAdapter::create(profile);
    adapter.keygen(rotation_steps, true);

    std::string scaling_best_mode;
    auto scaling_best_period_mode = m2424::BootstrapPeriodMode::NoBootstrapPeriod;
    double scaling_best_manual_period_log2 = 0.0;
    double scaling_best_plain_scale_log2 = profile.scale > 0.0 ? std::log2(profile.scale) : 40.0;
    if (!scaling_gate_has_pass(adapter,
                               profile,
                               slots,
                               tolerance,
                               amplitudes,
                               level_drops,
                               period_cases,
                               plain_scale_log2_values,
                               scaling_best_mode,
                               scaling_best_period_mode,
                               scaling_best_manual_period_log2,
                               scaling_best_plain_scale_log2)) {
        std::printf("summary,total_cases,pass_cases,fail_cases,max_final_error,max_evalmod_error,best_mode\n");
        std::printf("blocked_by_scaling_gate,0,0,0,0.000000e+00,0.000000e+00,none\n");
        return 0;
    }
    period_cases = {{scaling_best_period_mode, scaling_best_manual_period_log2}};
    plain_scale_log2_values = {scaling_best_plain_scale_log2};

    std::printf("case,input_kind,amplitude,level_drop,normalization_mode,denormalization_position,evalmod_degree,period_mode,manual_period_log2,plain_scale_log2,chain_before,chain_after,mod_raise_chain_before,mod_raise_chain_after,mod_raise_coeff_modulus_size_before,mod_raise_coeff_modulus_size_after,mod_raise_scale_before,mod_raise_scale_after,bootstrap_period_log2,bootstrap_period,bootstrap_scaling_factor,normalization_factor_log2,factor_times_plain_scale_log2,normalization_scalar_representable,denormalization_scalar_representable,max_abs_input,mod_raise_decoded_max_abs,mod_raise_diagnostic_error,max_abs_after_coeff_to_slot,normalization_factor,max_abs_after_normalization,inside_evalmod_interval,coeff_to_slot_error,normalization_error,evalmod_error,denormalization_error,slot_to_coeff_error,post_refresh_mod_raise_error,final_error,restore_ok,post_depth,post_ops_ok,post_error,exception,status\n");

    std::size_t case_id = 0;
    std::size_t pass_cases = 0;
    std::size_t fail_cases = 0;
    double max_final_error = 0.0;
    double max_evalmod_error = 0.0;
    double best_final_error = 0.0;
    std::string best_mode = scaling_best_mode;
    for (const auto& kind : kinds) {
        for (double amplitude : amplitudes) {
            for (std::size_t level_drop : level_drops) {
                for (const auto denorm_position : denorm_positions) {
                    for (const auto evalmod_degree : evalmod_degrees) {
                        for (const auto& period_case : period_cases) {
                            for (double plain_scale_log2 : plain_scale_log2_values) {
                                ++case_id;
                                try {
                            const auto input = make_input(kind, slots, amplitude);
                            const auto expected = to_complex(input);

                            auto current = adapter.encrypt(adapter.encode(input));
                            for (std::size_t i = 0; i < level_drop; ++i) {
                                current = adapter.mul_plain_rescale(current, adapter.encode_scalar_like(1.0, current));
                            }
                            const auto before = adapter.info(current);

                            m2424::BootstrapPrototype bootstrapper(adapter, slots, tolerance);
                            bootstrapper.set_normalization_mode(m2424::BootstrapNormalizationMode::PlainMultiplyRescale);
                            bootstrapper.set_denormalization_position(denorm_position);
                            bootstrapper.set_evalmod_degree(evalmod_degree);
                            bootstrapper.set_period_mode(period_case.mode);
                            bootstrapper.set_manual_period_log2(period_case.manual_period_log2);
                            bootstrapper.set_plain_scale_log2(plain_scale_log2);
                            bootstrapper.set_post_refresh_mod_raise_enabled(false);
                            auto checked = bootstrapper.refresh_cipher_checked(current, expected);
                            const auto after = adapter.info(checked.result);

                            bool post_ops_ok = true;
                            auto post = checked.result;
                            for (std::size_t i = 0; i < post_depth; ++i) {
                                try {
                                    post = adapter.mul_plain_rescale(post, adapter.encode_scalar_like(1.0, post));
                                } catch (...) {
                                    post_ops_ok = false;
                                    break;
                                }
                            }

                            double post_error = 0.0;
                            if (post_ops_ok) {
                                try {
                                    auto decoded = adapter.decode(adapter.decrypt(post));
                                    decoded.resize(slots);
                                    post_error = max_error(input, decoded);
                                } catch (...) {
                                    post_ops_ok = false;
                                }
                            }

                            const auto* mod_raise_stage = find_stage(checked, "mod_raise");
                            const double coeff_to_slot_error = stage_error(checked, "coeff_to_slot");
                            const double normalization_error = stage_error(checked, "eval_mod_normalization");
                            const double evalmod_error = stage_error(checked, "eval_mod");
                            const double denormalization_error = stage_error(checked, "refresh_denormalization");
                            const double slot_to_coeff_error = stage_error(checked, "slot_to_coeff");
                            const double post_refresh_mod_raise_error = stage_error(checked, "post_refresh_mod_raise");
                            const double final_error = stage_error(checked, "refresh_result");
                            const bool restore_ok = after.chain_index > before.chain_index;
                            const bool ok =
                                final_error <= tolerance &&
                                post_ops_ok &&
                                post_error <= tolerance;

                            max_final_error = std::max(max_final_error, final_error);
                            max_evalmod_error = std::max(max_evalmod_error, evalmod_error);
                            if (ok) {
                                ++pass_cases;
                                if (best_mode == "none" || final_error < best_final_error) {
                                    best_final_error = final_error;
                                    best_mode = std::string(m2424::to_string(checked.normalization_mode)) + "/" +
                                                m2424::to_string(checked.denormalization_position) + "/" +
                                                m2424::to_string(checked.evalmod_degree) + "/" +
                                                m2424::to_string(checked.period_mode) + "/plain_scale_log2=" +
                                                std::to_string(checked.plain_scale_log2);
                                }
                            } else {
                                ++fail_cases;
                            }

                            std::printf("%zu,%s,%.6e,%zu,%s,%s,%s,%s,%.0f,%.0f,%zu,%zu,%zu,%zu,%zu,%zu,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%s,%s,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%s,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%s,%zu,%s,%.6e,%s,%s\n",
                                        case_id,
                                        kind.c_str(),
                                        amplitude,
                                        level_drop,
                                        m2424::to_string(checked.normalization_mode),
                                        m2424::to_string(checked.denormalization_position),
                                        m2424::to_string(checked.evalmod_degree),
                                        m2424::to_string(checked.period_mode),
                                        checked.manual_period_log2,
                                        checked.plain_scale_log2,
                                        before.chain_index,
                                        after.chain_index,
                                        mod_raise_stage ? mod_raise_stage->chain_before : 0,
                                        mod_raise_stage ? mod_raise_stage->chain_after : 0,
                                        mod_raise_stage ? mod_raise_stage->coeff_modulus_size_before : 0,
                                        mod_raise_stage ? mod_raise_stage->coeff_modulus_size_after : 0,
                                        mod_raise_stage ? mod_raise_stage->scale_before : 0.0,
                                        mod_raise_stage ? mod_raise_stage->scale_after : 0.0,
                                        checked.bootstrap_period_log2,
                                        checked.bootstrap_period,
                                        checked.bootstrap_scaling_factor,
                                        checked.normalization_factor_log2,
                                        checked.factor_times_plain_scale_log2,
                                        checked.normalization_scalar_representable ? "true" : "false",
                                        checked.denormalization_scalar_representable ? "true" : "false",
                                        checked.max_abs_input,
                                        checked.max_abs_after_mod_raise_decode,
                                        checked.mod_raise_diagnostic_error,
                                        checked.max_abs_after_coeff_to_slot,
                                        checked.normalization_factor,
                                        checked.max_abs_after_normalization,
                                        checked.inside_evalmod_interval ? "PASS" : "FAIL",
                                        coeff_to_slot_error,
                                        normalization_error,
                                        evalmod_error,
                                        denormalization_error,
                                        slot_to_coeff_error,
                                        post_refresh_mod_raise_error,
                                        final_error,
                                        restore_ok ? "PASS" : "FAIL",
                                        post_depth,
                                        post_ops_ok ? "PASS" : "FAIL",
                                        post_error,
                                        "",
                                        ok ? "PASS" : "FAIL");
                                } catch (const std::exception& e) {
                            ++fail_cases;
                            std::string exception = e.what();
                            std::replace(exception.begin(), exception.end(), ',', ';');
                            std::printf("%zu,%s,%.6e,%zu,%s,%s,%s,%s,%.0f,%.0f,0,0,0,0,0,0,0,0,0,0,0,0,0,false,false,0,0,0,0,0,0,FAIL,0,0,0,0,0,0,0,FAIL,%zu,FAIL,0,%s,FAIL\n",
                                        case_id,
                                        kind.c_str(),
                                        amplitude,
                                        level_drop,
                                        m2424::to_string(m2424::BootstrapNormalizationMode::PlainMultiplyRescale),
                                        m2424::to_string(denorm_position),
                                        m2424::to_string(evalmod_degree),
                                        m2424::to_string(period_case.mode),
                                        period_case.manual_period_log2,
                                        plain_scale_log2,
                                        post_depth,
                                        exception.c_str());
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    std::printf("summary,total_cases,pass_cases,fail_cases,max_final_error,max_evalmod_error,best_mode\n");
    std::printf("summary,%zu,%zu,%zu,%.6e,%.6e,%s\n",
                case_id,
                pass_cases,
                fail_cases,
                max_final_error,
                max_evalmod_error,
                best_mode.c_str());

    return 0;
}
