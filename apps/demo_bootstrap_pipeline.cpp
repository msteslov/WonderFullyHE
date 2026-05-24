#include "m2424/bootstrap.hpp"
#include "m2424/profiles.hpp"

#include <cstdio>
#include <string>
#include <vector>

static void print_metrics(const char* label, const m2424::BootstrapCipherMetrics& metrics) {
    if (!metrics.available) {
        std::printf("%s,unavailable,0,0,0,0,0\n", label);
        return;
    }
    std::printf("%s,available,%.6e,%zu,%zu,%zu,%zu\n",
                label,
                metrics.scale,
                metrics.chain_index,
                metrics.coeff_modulus_size,
                metrics.ciphertext_size,
                metrics.serialized_bytes);
}

int main() {
    const std::size_t payload_size = 32;
    const auto profile = m2424::profiles::depth_ckks();

    auto adapter = m2424::SealAdapter::create(profile);
    adapter.keygen(true, true);

    std::vector<double> input;
    input.reserve(payload_size);
    for (std::size_t i = 0; i < payload_size; ++i) {
        input.push_back(0.25 + 0.25 * std::sin(static_cast<double>(i) / 7.0));
    }

    m2424::Bootstrapper bootstrapper(adapter);
    auto report = bootstrapper.analyze_depth(input, 8);

    std::printf("bootstrap_pipeline_status\n");
    std::printf("profile,poly_modulus_degree,payload_size,scale,coeff_modulus_bits\n");
    std::printf("depth_ckks,%zu,%zu,%.6e,60-40-40-40-40-60\n",
                profile.poly_modulus_degree, payload_size, profile.scale);
    std::printf("successful_multiplications=%zu\n", report.successful_multiplications);
    std::printf("next_exponent=%zu\n", report.next_exponent * 2);
    std::printf("depth_stop_reason=%s\n", report.stop_reason.c_str());
    std::printf("criterion,status\n");
    std::printf("Dec(c_prime)_approx_Dec(c),%s\n", report.preserve_value_criterion ? "PASS" : "IN_PROGRESS");
    std::printf("level(c_prime)_gt_level(c),%s\n", report.restore_level_criterion ? "PASS" : "IN_PROGRESS");
    std::printf("metrics_label,state,scale,chain_index,coeff_modulus_size,ciphertext_size,serialized_bytes\n");
    print_metrics("input", report.input);
    print_metrics("depth_boundary", report.depth_boundary);
    std::printf("stage,status,before_chain,after_chain,before_scale,after_scale,note\n");
    for (const auto& stage : report.stages) {
        const auto before_chain = stage.before.available ? std::to_string(stage.before.chain_index) : "";
        const auto after_chain = stage.after.available ? std::to_string(stage.after.chain_index) : "";
        const double before_scale = stage.before.available ? stage.before.scale : 0.0;
        const double after_scale = stage.after.available ? stage.after.scale : 0.0;
        std::printf("%s,%s,%s,%s,%.6e,%.6e,%s\n",
                    stage.name.c_str(),
                    m2424::to_string(stage.status),
                    before_chain.c_str(),
                    after_chain.c_str(),
                    before_scale,
                    after_scale,
                    stage.note.c_str());
    }

    return report.successful_multiplications > 0 ? 0 : 1;
}
