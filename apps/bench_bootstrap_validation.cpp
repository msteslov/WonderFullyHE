#include "m2424/bootstrap.hpp"
#include "m2424/profiles.hpp"
#include "m2424/seal_adapter.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
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

} // namespace

int main() {
    constexpr std::size_t slots = 16;
    constexpr double tolerance = 2e-5;
    constexpr std::size_t post_depth = 2;

    const std::vector<double> amplitudes{1e-6, 3e-6, 1e-5};
    const std::vector<std::string> kinds{"sine", "alternating", "ramp"};
    const std::vector<std::size_t> level_drops{1, 2};

    auto rotation_steps = m2424::Bootstrapper::refresh_rotation_steps(slots);
    auto adapter = m2424::SealAdapter::create(m2424::profiles::boot_ckks());
    adapter.keygen(rotation_steps, true);
    m2424::Bootstrapper bootstrapper(adapter);

    std::printf("case,input_kind,amplitude,level_drop,chain_before,chain_after,normalization_factor,preserve_error,restore_ok,post_depth,post_ops_ok,post_error,status\n");

    std::size_t case_id = 0;
    for (const auto& kind : kinds) {
        for (double amplitude : amplitudes) {
            for (std::size_t level_drop : level_drops) {
                ++case_id;
                const auto input = make_input(kind, slots, amplitude);
                const auto expected = to_complex(input);

                auto current = adapter.encrypt(adapter.encode(input));
                for (std::size_t i = 0; i < level_drop; ++i) {
                    current = adapter.mul_plain_rescale(current, adapter.encode_scalar_like(1.0, current));
                }
                const auto before = adapter.info(current);

                auto checked = bootstrapper.refresh_checked(current, expected, slots, tolerance);
                auto refreshed = bootstrapper.refresh(current, slots, tolerance);
                const auto after = adapter.info(refreshed.result);

                bool post_ops_ok = true;
                auto post = refreshed.result;
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
                    auto decoded = adapter.decode(adapter.decrypt(post));
                    decoded.resize(slots);
                    post_error = max_error(input, decoded);
                }

                const auto* final_stage = find_stage(checked, "refresh_result");
                const double preserve_error = final_stage ? final_stage->max_abs_error : 0.0;
                const bool restore_ok = after.chain_index > before.chain_index;
                const bool preserve_ok = preserve_error <= tolerance;
                const bool post_error_ok = post_ops_ok && post_error <= tolerance;
                const bool ok = preserve_ok && restore_ok && post_error_ok;

                std::printf("%zu,%s,%.6e,%zu,%zu,%zu,%.6e,%.6e,%s,%zu,%s,%.6e,%s\n",
                            case_id,
                            kind.c_str(),
                            amplitude,
                            level_drop,
                            before.chain_index,
                            after.chain_index,
                            checked.normalization_factor,
                            preserve_error,
                            restore_ok ? "PASS" : "FAIL",
                            post_depth,
                            post_ops_ok ? "PASS" : "FAIL",
                            post_error,
                            ok ? "PASS" : "FAIL");
            }
        }
    }

    return 0;
}
