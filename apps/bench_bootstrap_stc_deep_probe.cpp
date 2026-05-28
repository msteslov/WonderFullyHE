#include "m2424/bootstrap_plan.hpp"
#include "m2424/bootstrap_prototype.hpp"
#include "m2424/profiles.hpp"
#include "m2424/seal_adapter.hpp"
#include "m2424/security_report.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <string>
#include <vector>

namespace {

std::vector<double> make_input(std::size_t slots) {
    std::vector<double> input;
    input.reserve(slots);
    for (std::size_t i = 0; i < slots; ++i) {
        input.push_back(1e-5 * std::sin(static_cast<double>(i) / 4.0));
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

double post_error(m2424::SealAdapter& adapter,
                  const m2424::Cipher& input,
                  const std::vector<double>& expected) {
    auto post = adapter.mul_plain_rescale(input, adapter.encode_scalar_like(1.0, input));
    auto decoded = adapter.decode(adapter.decrypt(post));
    decoded.resize(expected.size());
    double error = 0.0;
    for (std::size_t i = 0; i < expected.size(); ++i) {
        error = std::max(error, std::abs(decoded[i] - expected[i]));
    }
    return error;
}

double final_error_from_report(const m2424::BootstrapPrototypeReport& report) {
    for (const auto& stage : report.stages) {
        if (stage.name == "refresh_result") {
            return stage.max_abs_error;
        }
    }
    return 0.0;
}

void sanitize(std::string& text) {
    std::replace(text.begin(), text.end(), ',', ';');
}

bool is_sweep_mode(int argc, char** argv) {
    return argc > 1 && std::string(argv[1]) == "sweep";
}

std::size_t arg_size_or(int argc, char** argv, int index, std::size_t fallback) {
    if (argc <= index) {
        return fallback;
    }
    const auto value = std::strtoull(argv[index], nullptr, 10);
    return value == 0 ? fallback : static_cast<std::size_t>(value);
}

std::vector<std::size_t> target_chain_values(bool sweep, std::size_t single) {
    return sweep ? std::vector<std::size_t>{10, 11, 12, 13} : std::vector<std::size_t>{single};
}

std::vector<double> period_offset_values(bool sweep, double single) {
    return sweep ? std::vector<double>{-4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0} : std::vector<double>{single};
}

double arg_double_or(int argc, char** argv, int index, double fallback) {
    if (argc <= index) {
        return fallback;
    }
    char* end = nullptr;
    const double value = std::strtod(argv[index], &end);
    return end == argv[index] ? fallback : value;
}

void print_stage(const char* stage,
                 const char* status,
                 const m2424::SecurityReport& security,
                 std::size_t target_chain,
                 double plain_scale_log2,
                 double period_offset_log2,
                 const m2424::CipherInfo& info,
                 double max_abs_error,
                 const std::string& exception = {}) {
    std::printf("boot_deep_ckks,%s,%d,%d,%zu,%.0f,%.0f,%s,%s,%zu,%.6f,%.6f,%.6e,%s\n",
                m2424::to_string(security.effective_level),
                security.total_coeff_modulus_bits,
                security.tc128_limit,
                target_chain,
                plain_scale_log2,
                period_offset_log2,
                stage,
                status,
                info.chain_index,
                std::log2(info.scale),
                info.coeff_modulus_log2,
                max_abs_error,
                exception.c_str());
    std::fflush(stdout);
}

void print_report_stage(const m2424::BootstrapPrototypeStage& stage,
                        const m2424::SecurityReport& security,
                        std::size_t target_chain,
                        double plain_scale_log2,
                        double period_offset_log2) {
    std::printf("boot_deep_ckks,%s,%d,%d,%zu,%.0f,%.0f,%s,%s,%zu,%.6f,%.6f,%.6e,\n",
                m2424::to_string(security.effective_level),
                security.total_coeff_modulus_bits,
                security.tc128_limit,
                target_chain,
                plain_scale_log2,
                period_offset_log2,
                stage.name.c_str(),
                stage.status.c_str(),
                stage.chain_after,
                std::log2(stage.scale_after),
                stage.coeff_modulus_log2_after,
                stage.max_abs_error);
    std::fflush(stdout);
}

} // namespace

int main(int argc, char** argv) {
    constexpr std::size_t slots = 16;
    constexpr double tolerance = 2e-5;
    const bool sweep = is_sweep_mode(argc, argv);
    const std::size_t target_chain = sweep ? 12 : arg_size_or(argc, argv, 1, 12);
    const double plain_scale_log2 = sweep ? 40.0 : arg_double_or(argc, argv, 2, 40.0);
    const double period_offset_log2 = sweep ? 3.0 : arg_double_or(argc, argv, 3, 3.0);
    const auto profile = m2424::profiles::boot_deep_ckks();
    const auto security = m2424::analyze_security("boot_deep_ckks", profile);
    auto rotations = m2424::bootstrap_plan_rotation_steps(m2424::make_scalable_bootstrap_plan(slots));

    std::printf("profile,effective_security,total_coeff_modulus_bits,tc128_limit,target_chain,plain_scale_log2,period_offset_log2,stage,status,chain_index,scale_log2,coeff_modulus_log2,max_abs_error,exception\n");
    std::fflush(stdout);

    auto adapter = m2424::SealAdapter::create(profile);
    std::printf("boot_deep_ckks,%s,%d,%d,%zu,%.0f,%.0f,keygen_start,RUN,0,0,0,0,\n",
                m2424::to_string(security.effective_level),
                security.total_coeff_modulus_bits,
                security.tc128_limit,
                target_chain,
                plain_scale_log2,
                period_offset_log2);
    std::fflush(stdout);
    adapter.keygen(rotations, true);
    std::printf("boot_deep_ckks,%s,%d,%d,%zu,%.0f,%.0f,keygen_done,PASS,0,0,0,0,\n",
                m2424::to_string(security.effective_level),
                security.total_coeff_modulus_bits,
                security.tc128_limit,
                target_chain,
                plain_scale_log2,
                period_offset_log2);
    std::fflush(stdout);

    const auto input = make_input(slots);
    const auto expected = to_complex(input);
    for (const auto target : target_chain_values(sweep, target_chain)) {
        for (const double offset : period_offset_values(sweep, period_offset_log2)) {
            try {
                auto current = adapter.encrypt(adapter.encode(input));
                if (!sweep) {
                    print_stage("encrypt", "PASS", security, target, plain_scale_log2, offset,
                                adapter.info(current), 0.0);
                }

                m2424::BootstrapPrototype bootstrapper(adapter, slots, tolerance);
                bootstrapper.set_circuit_order(m2424::BootstrapCircuitOrder::SlotsToCoeffsFirst);
                bootstrapper.set_transform_backend(m2424::BootstrapTransformBackend::FftLike);
                bootstrapper.set_evalmod_degree(m2424::EvalModDegree::P3);
                bootstrapper.set_period_mode(m2424::BootstrapPeriodMode::SourceCoeffModulus);
                bootstrapper.set_plain_scale_log2(plain_scale_log2);
                bootstrapper.set_stc_first_target_chain_index(target);
                bootstrapper.set_stc_first_period_offset_log2(offset);
                const auto report = bootstrapper.refresh_cipher_checked(current, expected);
                if (!sweep) {
                    for (const auto& stage : report.stages) {
                        print_report_stage(stage, security, target, plain_scale_log2, offset);
                    }
                }

                std::string status = "PASS";
                double error = 0.0;
                try {
                    error = post_error(adapter, report.result, input);
                } catch (const std::exception& e) {
                    status = e.what();
                    sanitize(status);
                }
                const auto result_info = adapter.info(report.result);
                const double final_error = final_error_from_report(report);
                if (sweep) {
                    std::printf("boot_deep_ckks,%s,%d,%d,%zu,%.0f,%.0f,sweep_result,%s,%zu,%.6f,%.6f,%.6e,post_error=%.6e\n",
                                m2424::to_string(security.effective_level),
                                security.total_coeff_modulus_bits,
                                security.tc128_limit,
                                target,
                                plain_scale_log2,
                                offset,
                                final_error <= tolerance && status == "PASS" && error <= tolerance ? "PASS" : "FAIL",
                                result_info.chain_index,
                                std::log2(result_info.scale),
                                result_info.coeff_modulus_log2,
                                final_error,
                                error);
                    std::fflush(stdout);
                } else {
                    print_stage("post_refresh_mul_plain_rescale",
                                status.c_str(),
                                security,
                                target,
                                plain_scale_log2,
                                offset,
                                result_info,
                                error);
                }
            } catch (const std::exception& e) {
                auto exception = std::string(e.what());
                sanitize(exception);
                m2424::CipherInfo empty;
                print_stage("refresh", "FAIL", security, target, plain_scale_log2, offset,
                            empty, 0.0, exception);
            }
        }
    }
    return 0;
}
