#include "m2424/bootstrap_layout_v2.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>

namespace m2424 {
namespace {

std::size_t power_of_two(std::size_t log_value, const char* name) {
    constexpr std::size_t kBits = std::numeric_limits<std::size_t>::digits;
    if (log_value >= kBits) {
        throw std::invalid_argument(std::string(name) + " is too large");
    }
    return std::size_t{1} << log_value;
}

std::size_t sum_bits(const std::vector<int>& bits) {
    return static_cast<std::size_t>(std::accumulate(bits.begin(), bits.end(), 0));
}

std::vector<int> repeated_bits(int bits, std::size_t count) {
    return std::vector<int>(count, bits);
}

BootstrapChainSegment make_segment(std::string name,
                                   std::vector<int> log_q,
                                   double log_scale) {
    BootstrapChainSegment segment;
    segment.name = std::move(name);
    segment.log_q = std::move(log_q);
    segment.log_scale = log_scale;
    segment.levels = segment.log_q.size();
    return segment;
}

void validate_request(std::size_t log_n,
                      std::size_t log_slots,
                      double target_error,
                      std::size_t cycles) {
    if (log_n < 4) {
        throw std::invalid_argument("bootstrap layout log_n must be at least 4");
    }
    if (log_slots >= log_n) {
        throw std::invalid_argument("bootstrap layout log_slots must be smaller than log_n");
    }
    if (!std::isfinite(target_error) || target_error <= 0.0) {
        throw std::invalid_argument("bootstrap layout target_error must be positive and finite");
    }
    if (cycles == 0) {
        throw std::invalid_argument("bootstrap layout cycles must be positive");
    }
}

std::vector<int> concat_log_q(const BootstrapChainSegment& a,
                              const BootstrapChainSegment& b,
                              const BootstrapChainSegment& c,
                              const BootstrapChainSegment& d,
                              const BootstrapChainSegment& e) {
    std::vector<int> result;
    result.reserve(a.log_q.size() + b.log_q.size() + c.log_q.size() + d.log_q.size() + e.log_q.size());
    result.insert(result.end(), a.log_q.begin(), a.log_q.end());
    result.insert(result.end(), b.log_q.begin(), b.log_q.end());
    result.insert(result.end(), c.log_q.begin(), c.log_q.end());
    result.insert(result.end(), d.log_q.begin(), d.log_q.end());
    result.insert(result.end(), e.log_q.begin(), e.log_q.end());
    return result;
}

} // namespace

BootstrapLayoutV2 make_lattigo_like_bootstrap_layout_v2(std::size_t log_n,
                                                        std::size_t log_slots,
                                                        double target_error,
                                                        std::size_t cycles) {
    validate_request(log_n, log_slots, target_error, cycles);

    BootstrapLayoutV2 layout;
    layout.slots_to_coeffs = make_segment("slots_to_coeffs", {60, 39}, 39.0);
    layout.modup = make_segment("modup", {60}, 60.0);
    layout.coeffs_to_slots = make_segment("coeffs_to_slots", {56, 56}, 56.0);

    constexpr std::size_t kCosDiscreteDegree = 30;
    constexpr std::size_t kDoubleAngle = 3;
    const std::size_t evalmod_levels = 2 + kDoubleAngle + static_cast<std::size_t>(std::ceil(std::log2(kCosDiscreteDegree)));
    layout.evalmod = make_segment("evalmod_degree30_double_angle3", repeated_bits(60, evalmod_levels), 60.0);

    const std::size_t residual_levels = std::max<std::size_t>(6, cycles + 3);
    std::vector<int> output_bits;
    output_bits.reserve(residual_levels);
    output_bits.push_back(60);
    for (std::size_t i = 1; i + 1 < residual_levels; ++i) {
        output_bits.push_back(50);
    }
    output_bits.push_back(60);
    layout.output = make_segment("residual_output", std::move(output_bits), 50.0);

    layout.key_switch.log_p = {61, 61, 61};
    layout.key_switch.decomposition_count = layout.key_switch.log_p.size();
    layout.key_switch.uses_special_primes = true;
    layout.key_switch.uses_ephemeral_sparse_secret = true;
    layout.key_switch.ephemeral_secret_weight = 32;

    const auto all_q = concat_log_q(layout.slots_to_coeffs,
                                    layout.modup,
                                    layout.coeffs_to_slots,
                                    layout.evalmod,
                                    layout.output);
    layout.total_log_q = sum_bits(all_q);
    layout.total_log_p = sum_bits(layout.key_switch.log_p);
    layout.minimum_input_level = layout.slots_to_coeffs.levels
        + layout.modup.levels
        + layout.coeffs_to_slots.levels
        + layout.evalmod.levels;
    layout.expected_output_level = layout.output.levels - 1;

    layout.residual_profile.polyModulusDegree = power_of_two(log_n, "log_n");
    layout.residual_profile.coeffModulusBits = layout.output.log_q;
    layout.residual_profile.scale = std::exp2(layout.output.log_scale);
    layout.residual_profile.slots = power_of_two(log_slots, "log_slots");
    return layout;
}

} // namespace m2424
