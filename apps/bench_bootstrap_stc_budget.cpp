#include "m2424/bootstrap_plan.hpp"
#include "m2424/eval_mod.hpp"
#include "m2424/profiles.hpp"
#include "m2424/security_report.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace {

int total_bits(const std::vector<int>& bits, std::size_t active_size) {
    int total = 0;
    active_size = std::min(active_size, bits.size());
    for (std::size_t i = 0; i < active_size; ++i) {
        total += bits[i];
    }
    return total;
}

void run_profile(const char* profile_name,
                 const m2424::CkksProfile& profile,
                 std::size_t slots,
                 std::size_t target_chain,
                 double plain_scale_log2,
                 double period_offset_log2) {
    const auto security = m2424::analyze_security(profile_name, profile);
    const std::size_t top_chain = profile.coeff_modulus_bits.size() - 2;
    if (target_chain > top_chain) {
        return;
    }

    constexpr std::size_t slot_to_coeff_levels = 1;
    constexpr std::size_t coeff_to_slot_levels = 1;
    constexpr std::size_t evalmod_p3_levels = 3;
    const std::size_t chain_before_slot_to_coeff = target_chain;
    if (chain_before_slot_to_coeff < slot_to_coeff_levels) {
        return;
    }
    const std::size_t chain_before_mod_raise = chain_before_slot_to_coeff - slot_to_coeff_levels;
    const double period_log2 =
        static_cast<double>(total_bits(profile.coeff_modulus_bits, chain_before_mod_raise + 1)) -
        period_offset_log2;
    const std::size_t normalization_chunks = static_cast<std::size_t>(
        std::ceil(std::max(0.0, period_log2) / plain_scale_log2));
    const std::size_t chain_before_normalization = top_chain - coeff_to_slot_levels;
    const bool normalization_levels_ready = chain_before_normalization >= normalization_chunks;
    const std::size_t chain_before_evalmod = normalization_levels_ready
        ? chain_before_normalization - normalization_chunks
        : 0;
    const bool evalmod_levels_ready = chain_before_evalmod >= evalmod_p3_levels;
    const std::size_t chain_after_evalmod = evalmod_levels_ready
        ? chain_before_evalmod - evalmod_p3_levels
        : 0;
    const bool continuation_ready = chain_after_evalmod >= 1;

    std::printf("%s,%zu,%zu,%zu,%.0f,%.0f,%zu,%zu,%zu,%.6f,%zu,%zu,%zu,%zu,%zu,%s,%s,%s,%s,%s,%d,%d\n",
                profile_name,
                slots,
                top_chain,
                target_chain,
                plain_scale_log2,
                period_offset_log2,
                chain_before_mod_raise,
                top_chain,
                chain_before_normalization,
                period_log2,
                normalization_chunks,
                chain_before_evalmod,
                evalmod_p3_levels,
                chain_after_evalmod,
                chain_after_evalmod,
                normalization_levels_ready ? "true" : "false",
                evalmod_levels_ready ? "true" : "false",
                continuation_ready ? "true" : "false",
                security.passes_tc128 ? "true" : "false",
                m2424::to_string(security.effective_level),
                security.total_coeff_modulus_bits,
                security.tc128_limit);
}

} // namespace

int main() {
    constexpr std::size_t slots = 16;
    constexpr double plain_scale_log2 = 40.0;
    constexpr double period_offset_log2 = 3.0;
    std::printf("profile,slots,top_chain,target_chain,plain_scale_log2,period_offset_log2,chain_before_mod_raise,chain_after_mod_raise,chain_before_normalization,period_log2,normalization_chunks,chain_before_evalmod,evalmod_p3_levels,chain_after_evalmod,chain_remaining_after_refresh,normalization_levels_ready,evalmod_levels_ready,continuation_ready,tc128_ready,effective_security,total_coeff_modulus_bits,tc128_limit\n");
    for (std::size_t target_chain = 2; target_chain <= 12; ++target_chain) {
        run_profile("boot_ckks", m2424::profiles::boot_ckks(), slots, target_chain, plain_scale_log2, period_offset_log2);
    }
    for (std::size_t target_chain = 2; target_chain <= 18; ++target_chain) {
        run_profile("boot_deep_ckks", m2424::profiles::boot_deep_ckks(), slots, target_chain, plain_scale_log2, period_offset_log2);
    }
    return 0;
}
