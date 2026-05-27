#include "m2424/bootstrap_prototype.hpp"
#include "m2424/profiles.hpp"
#include "m2424/seal_adapter.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>
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

    auto rotation_steps = m2424::BootstrapPrototype::required_rotation_steps(slots);
    auto adapter = m2424::SealAdapter::create(m2424::profiles::boot_ckks());
    adapter.keygen(rotation_steps, true);

    std::printf("case,input_kind,amplitude,level_drop,normalization_mode,denormalization_position,evalmod_degree,chain_before,chain_after,mod_raise_chain_before,mod_raise_chain_after,mod_raise_coeff_modulus_size_before,mod_raise_coeff_modulus_size_after,mod_raise_scale_before,mod_raise_scale_after,bootstrap_period_log2,bootstrap_period,bootstrap_scaling_factor,max_abs_input,mod_raise_decoded_max_abs,mod_raise_diagnostic_error,max_abs_after_coeff_to_slot,normalization_factor,max_abs_after_normalization,inside_evalmod_interval,coeff_to_slot_error,normalization_error,evalmod_error,denormalization_error,slot_to_coeff_error,post_refresh_mod_raise_error,final_error,restore_ok,post_depth,post_ops_ok,post_error,status\n");

    std::size_t case_id = 0;
    std::size_t pass_cases = 0;
    std::size_t fail_cases = 0;
    double max_final_error = 0.0;
    double max_evalmod_error = 0.0;
    double best_final_error = 0.0;
    std::string best_mode = "none";
    for (const auto& kind : kinds) {
        for (double amplitude : amplitudes) {
            for (std::size_t level_drop : level_drops) {
                for (const auto denorm_position : denorm_positions) {
                    for (const auto evalmod_degree : evalmod_degrees) {
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
                                restore_ok &&
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
                                                m2424::to_string(checked.evalmod_degree);
                                }
                            } else {
                                ++fail_cases;
                            }

                            std::printf("%zu,%s,%.6e,%zu,%s,%s,%s,%zu,%zu,%zu,%zu,%zu,%zu,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%s,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%s,%zu,%s,%.6e,%s\n",
                                        case_id,
                                        kind.c_str(),
                                        amplitude,
                                        level_drop,
                                        m2424::to_string(checked.normalization_mode),
                                        m2424::to_string(checked.denormalization_position),
                                        m2424::to_string(checked.evalmod_degree),
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
                                        ok ? "PASS" : "FAIL");
                        } catch (const std::exception&) {
                            ++fail_cases;
                            std::printf("%zu,%s,%.6e,%zu,%s,%s,%s,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,FAIL,0,0,0,0,0,0,0,FAIL,%zu,FAIL,0,FAIL\n",
                                        case_id,
                                        kind.c_str(),
                                        amplitude,
                                        level_drop,
                                        m2424::to_string(m2424::BootstrapNormalizationMode::PlainMultiplyRescale),
                                        m2424::to_string(denorm_position),
                                        m2424::to_string(evalmod_degree),
                                        post_depth);
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
