#include "m2424/experimental.hpp"

#include <cstdio>
#include <exception>
#include <stdexcept>

namespace {

void require(bool condition, const char* message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

} // namespace

int main() {
    try {
        m2424::BootstrapExperimentConfig fft;
        fft.name = "fft";
        fft.slots = 8;
        fft.transform_backend = m2424::BootstrapTransformBackend::FftLike;
        const auto fft_steps = m2424::bootstrap_experiment_rotation_steps(fft);
        require(!fft_steps.empty(), "FFT experiment must require rotations");
        require(fft_steps == m2424::Bootstrapper::scalable_refresh_rotation_steps(fft.slots),
                "FFT experiment rotations must match scalable bootstrap rotations");

        auto dense = fft;
        dense.name = "dense";
        dense.transform_backend = m2424::BootstrapTransformBackend::DenseDiagonal;
        const auto dense_steps = m2424::bootstrap_experiment_rotation_steps(dense);
        require(!dense_steps.empty(), "dense experiment must require rotations");

        auto invalid = fft;
        invalid.name = "invalid";
        invalid.slots = 3;
        bool rejected = false;
        try {
            (void)m2424::bootstrap_experiment_rotation_steps(invalid);
        } catch (const std::invalid_argument&) {
            rejected = true;
        }
        require(rejected, "invalid experiment configuration must be rejected");

        std::printf("[test_bootstrap_experiment] PASS fft_rotations=%zu dense_rotations=%zu\n",
                    fft_steps.size(), dense_steps.size());
    } catch (const std::exception& error) {
        std::printf("[test_bootstrap_experiment] FAIL: %s\n", error.what());
        return 1;
    }
    return 0;
}
