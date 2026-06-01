#pragma once

#include "m2424/seal_adapter.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace m2424 {

struct BootstrapChainSegment {
    std::string name;
    std::vector<int> log_q;
    double log_scale{};
    std::size_t levels{};
};

struct BootstrapKeySwitchParams {
    std::vector<int> log_p;
    std::size_t decomposition_count{};
    bool uses_special_primes{};
    bool uses_ephemeral_sparse_secret{};
    std::size_t ephemeral_secret_weight{};
};

struct BootstrapLayoutV2 {
    CkksProfile residual_profile;
    BootstrapChainSegment slots_to_coeffs;
    BootstrapChainSegment modup;
    BootstrapChainSegment coeffs_to_slots;
    BootstrapChainSegment evalmod;
    BootstrapChainSegment output;
    BootstrapKeySwitchParams key_switch;
    std::size_t total_log_q{};
    std::size_t total_log_p{};
    std::size_t expected_output_level{};
    std::size_t minimum_input_level{};
};

BootstrapLayoutV2 make_lattigo_like_bootstrap_layout_v2(std::size_t log_n,
                                                        std::size_t log_slots,
                                                        double target_error,
                                                        std::size_t cycles);

} // namespace m2424
