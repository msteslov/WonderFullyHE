#include "m2424/bootstrap.hpp"

#include <cmath>
#include <cstdio>
#include <vector>

int main() {
    const std::size_t poly_degree = 16384;
    const std::size_t payload_size = 32;
    m2424::CkksProfile profile{poly_degree, {60, 40, 40, 40, 40, 60}, std::pow(2.0, 40), poly_degree / 2};

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
    std::printf("successful_multiplications=%zu\n", report.successful_multiplications);
    std::printf("next_exponent=%zu\n", report.next_exponent * 2);
    std::printf("depth_stop_reason=%s\n", report.stop_reason.c_str());
    std::printf("stage,status,note\n");
    for (const auto& stage : report.stages) {
        std::printf("%s,%s,%s\n", stage.name.c_str(), m2424::to_string(stage.status), stage.note.c_str());
    }

    return report.successful_multiplications > 0 ? 0 : 1;
}
