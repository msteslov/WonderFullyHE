#include "m2424/experimental.hpp"
#include "m2424/profiles.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

namespace {

m2424::ComplexVector make_input(std::size_t slots) {
    m2424::ComplexVector values;
    values.reserve(slots);
    for (std::size_t index = 0; index < slots; ++index) {
        values.push_back({1e-5 * std::sin(static_cast<double>(index + 1) / 3.0), 0.0});
    }
    return values;
}

} // namespace

int main() {
    constexpr std::size_t slots = 4;
    const auto expected = make_input(slots);

    m2424::BootstrapExperimentConfig fft_p3;
    fft_p3.name = "fft_p3";
    fft_p3.slots = slots;
    fft_p3.transform_backend = m2424::BootstrapTransformBackend::FftLike;
    fft_p3.evalmod_degree = m2424::EvalModDegree::P3;

    auto butterfly_p3 = fft_p3;
    butterfly_p3.name = "butterfly_p3";
    butterfly_p3.transform_backend = m2424::BootstrapTransformBackend::SmallSlots4Butterfly;

    const std::vector<m2424::BootstrapExperimentConfig> configs{fft_p3, butterfly_p3};
    auto adapter = m2424::SealAdapter::create(m2424::profiles::boot_deep_ckks());
    auto rotations = m2424::bootstrap_experiment_rotation_steps(fft_p3);
    const auto butterfly_rotations = m2424::bootstrap_experiment_rotation_steps(butterfly_p3);
    rotations.insert(rotations.end(), butterfly_rotations.begin(), butterfly_rotations.end());
    std::sort(rotations.begin(), rotations.end());
    rotations.erase(std::unique(rotations.begin(), rotations.end()), rotations.end());
    adapter.keygen(rotations, true);

    const auto input = adapter.encrypt(adapter.encode_complex(expected));
    const auto results = m2424::run_bootstrap_experiments(adapter, input, &expected, configs);

    std::printf("name,outcome,blocker,rotations,stages,continuation_levels,preserve_value,restore_levels\n");
    for (const auto& result : results) {
        std::printf("%s,%s,%s,%zu,%zu,%zu,%s,%s\n",
                    result.config.name.c_str(),
                    m2424::to_string(result.outcome),
                    result.blocker.c_str(),
                    result.rotation_steps.size(),
                    result.report.stages.size(),
                    result.report.continuation_levels,
                    result.report.preserve_value_criterion ? "true" : "false",
                    result.report.restore_level_criterion ? "true" : "false");
    }
    return 0;
}
