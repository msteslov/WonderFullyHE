#include "m2424/bootstrap_plan.hpp"
#include "m2424/bootstrap_prototype.hpp"
#include "m2424/profiles.hpp"
#include "m2424/seal_adapter.hpp"
#include "m2424/security_report.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
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

void sanitize(std::string& text) {
    std::replace(text.begin(), text.end(), ',', ';');
}

} // namespace

int main() {
    constexpr std::size_t slots = 16;
    constexpr double tolerance = 2e-5;
    const auto profile = m2424::profiles::boot_deep_ckks();
    const auto security = m2424::analyze_security("boot_deep_ckks", profile);
    const auto rotations = m2424::bootstrap_plan_rotation_steps(m2424::make_scalable_bootstrap_plan(slots));
    auto adapter = m2424::SealAdapter::create(profile);
    adapter.keygen(rotations, true);

    std::printf("profile,effective_security,total_coeff_modulus_bits,tc128_limit,target_chain,stage,status,chain_before,chain_after,scale_before_log2,scale_after_log2,max_abs_error,exception\n");
    const auto input = make_input(slots);
    const auto expected = to_complex(input);
    for (std::size_t target_chain : {2U, 6U, 10U, 12U}) {
        try {
            auto current = adapter.encrypt(adapter.encode(input));
            m2424::BootstrapPrototype bootstrapper(adapter, slots, tolerance);
            bootstrapper.set_circuit_order(m2424::BootstrapCircuitOrder::SlotsToCoeffsFirst);
            bootstrapper.set_transform_backend(m2424::BootstrapTransformBackend::FftLike);
            bootstrapper.set_evalmod_degree(m2424::EvalModDegree::P3);
            bootstrapper.set_period_mode(m2424::BootstrapPeriodMode::SourceCoeffModulus);
            bootstrapper.set_plain_scale_log2(40.0);
            bootstrapper.set_stc_first_target_chain_index(target_chain);
            const auto report = bootstrapper.refresh_cipher_checked(current, expected);
            for (const auto& stage : report.stages) {
                std::printf("boot_deep_ckks,%s,%d,%d,%zu,%s,%s,%zu,%zu,%.6f,%.6f,%.6e,\n",
                            m2424::to_string(security.effective_level),
                            security.total_coeff_modulus_bits,
                            security.tc128_limit,
                            target_chain,
                            stage.name.c_str(),
                            stage.status.c_str(),
                            stage.chain_before,
                            stage.chain_after,
                            std::log2(stage.scale_before),
                            std::log2(stage.scale_after),
                            stage.max_abs_error);
            }
            std::string status = "PASS";
            double error = 0.0;
            try {
                error = post_error(adapter, report.result, input);
            } catch (const std::exception& e) {
                status = e.what();
                sanitize(status);
            }
            const auto info = adapter.info(report.result);
            std::printf("boot_deep_ckks,%s,%d,%d,%zu,post_refresh_mul_plain_rescale,%s,%zu,%zu,%.6f,%.6f,%.6e,\n",
                        m2424::to_string(security.effective_level),
                        security.total_coeff_modulus_bits,
                        security.tc128_limit,
                        target_chain,
                        status.c_str(),
                        info.chain_index,
                        info.chain_index,
                        std::log2(info.scale),
                        std::log2(info.scale),
                        error);
        } catch (const std::exception& e) {
            auto exception = std::string(e.what());
            sanitize(exception);
            std::printf("boot_deep_ckks,%s,%d,%d,%zu,refresh,FAIL,0,0,0,0,0,%s\n",
                        m2424::to_string(security.effective_level),
                        security.total_coeff_modulus_bits,
                        security.tc128_limit,
                        target_chain,
                        exception.c_str());
        }
    }
    return 0;
}
