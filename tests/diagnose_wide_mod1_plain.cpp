#include "m2424/mod1_approximation.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>
#include <limits>
#include <vector>

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kTarget = 1e-10;

double reference_mod1(double x) {
    return std::sin(2.0 * kPi * x) / (2.0 * kPi);
}

std::size_t required_double_angle_steps(double input_bound, double compressed_bound) {
    if (input_bound <= compressed_bound) {
        return 0;
    }
    return static_cast<std::size_t>(std::ceil(std::log2(input_bound / compressed_bound)));
}

double run_case(std::size_t degree, double input_bound, double compressed_bound) {
    m2424::BootstrapWideMod1Plan plan;
    plan.input_bound = input_bound;
    plan.compressed_bound = compressed_bound;
    plan.double_angle_steps = required_double_angle_steps(input_bound, compressed_bound);
    plan.polynomial_degree = degree;
    plan.target_error = kTarget;
    plan.evalmod_scale_log2 = 59.0;
    const auto approximation = m2424::make_wide_mod1_approximation(plan);

    std::vector<double> integers{0.0, 1.0, -1.0, 7.0, -11.0};
    if (input_bound >= 1e3) {
        integers.push_back(997.0);
        integers.push_back(-997.0);
    }
    if (input_bound >= 1e6) {
        integers.push_back(999983.0);
        integers.push_back(-999983.0);
    }
    if (input_bound >= 1e9) {
        integers.push_back(999999937.0);
        integers.push_back(-999999937.0);
    }
    const std::vector<double> epsilons{0.0, 1e-6, -1e-6, 3e-6, -3e-6, 1e-5, -1e-5};

    double max_error = 0.0;
    for (double k : integers) {
        for (double epsilon : epsilons) {
            const double x = k + epsilon;
            if (std::abs(x) > input_bound) {
                continue;
            }
            const double actual = m2424::evaluate_wide_mod1_plain(approximation, x);
            const double expected = reference_mod1(x);
            max_error = std::max(max_error, std::abs(actual - expected));
        }
    }

    std::printf("[diagnose_wide_mod1_plain] degree=%zu double_angle_steps=%zu input_bound=%.12e compressed_bound=%.12e max_error=%.12e target=%.12e status=%s note=\"%s\"\n",
                degree,
                plan.double_angle_steps,
                input_bound,
                compressed_bound,
                max_error,
                kTarget,
                max_error <= kTarget ? "PASS" : "FAIL",
                approximation.construction_note.c_str());
    return max_error;
}

} // namespace

int main() {
    try {
        double worst = 0.0;
        for (double input_bound : {16.0, 1024.0}) {
            worst = std::max(worst, run_case(15, input_bound, 0.125));
            worst = std::max(worst, run_case(31, input_bound, 0.125));
        }
        if (worst > kTarget) {
            std::printf("[diagnose_wide_mod1_plain] status=FAIL worst_error=%.12e target=%.12e\n", worst, kTarget);
            return 1;
        }
        std::printf("[diagnose_wide_mod1_plain] status=PASS worst_error=%.12e target=%.12e\n", worst, kTarget);
        return 0;
    } catch (const std::exception& e) {
        std::printf("[diagnose_wide_mod1_plain] FAIL: %s\n", e.what());
        return 1;
    }
}
