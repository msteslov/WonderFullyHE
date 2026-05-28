#include <algorithm>
#include <cmath>
#include <cstdio>

namespace {

struct BudgetCase {
    double first_prime_log2{};
    double evalmod_prime_log2{};
    std::size_t evalmod_levels{};
    double active_modulus_log2{};
    double period_offset_log2{};
    double target_scale_log2{};
};

void run_case(const BudgetCase& c) {
    constexpr double security_budget_log2 = 881.0;
    constexpr double pre_evalmod_overhead_log2 = 100.0;
    constexpr double start_scale_log2 = 40.0;
    constexpr double evalmod_capacity_margin_log2 = 2.0;

    const double remaining_modulus_log2 =
        c.first_prime_log2 + static_cast<double>(c.evalmod_levels) * c.evalmod_prime_log2;
    const double period_log2 = c.active_modulus_log2 - c.period_offset_log2;
    const double required_drop_log2 = start_scale_log2 + period_log2 - c.target_scale_log2;
    const double available_drop_log2 = std::max(0.0, c.active_modulus_log2 - remaining_modulus_log2);
    const double remaining_after_strategy_log2 = c.active_modulus_log2 - required_drop_log2;
    const double first_evalmod_product_scale_log2 = 2.0 * c.target_scale_log2;
    const double max_target_scale_for_offset_log2 =
        c.period_offset_log2 - start_scale_log2 - evalmod_capacity_margin_log2;
    const double required_period_offset_log2 =
        start_scale_log2 + c.target_scale_log2 + evalmod_capacity_margin_log2;
    const double period_offset_margin_log2 =
        c.period_offset_log2 - required_period_offset_log2;
    const bool security_ready =
        c.active_modulus_log2 + pre_evalmod_overhead_log2 <= security_budget_log2;
    const bool algebraic_capacity_possible =
        c.target_scale_log2 <= max_target_scale_for_offset_log2;
    const bool scale_strategy_ready =
        required_drop_log2 <= available_drop_log2;
    const bool evalmod_capacity_ready =
        first_evalmod_product_scale_log2 + evalmod_capacity_margin_log2 <= remaining_after_strategy_log2;
    const bool feasible = security_ready && algebraic_capacity_possible &&
        scale_strategy_ready && evalmod_capacity_ready;
    const char* blocker = !security_ready
        ? "security_budget"
        : !algebraic_capacity_possible
        ? "period_offset_too_small_for_evalmod_scale"
        : !scale_strategy_ready
        ? "scale_strategy"
        : !evalmod_capacity_ready
        ? "evalmod_first_multiply"
        : "none";

    std::printf("%.0f,%.0f,%zu,%.0f,%.0f,%.0f,%.0f,%.0f,%.0f,%.0f,%.0f,%.0f,%.0f,%.0f,%.6e,%s,%s,%s,%s,%s,%s\n",
                c.first_prime_log2,
                c.evalmod_prime_log2,
                c.evalmod_levels,
                c.active_modulus_log2,
                c.active_modulus_log2 + pre_evalmod_overhead_log2,
                c.period_offset_log2,
                required_period_offset_log2,
                period_offset_margin_log2,
                max_target_scale_for_offset_log2,
                c.target_scale_log2,
                remaining_modulus_log2,
                required_drop_log2,
                available_drop_log2,
                remaining_after_strategy_log2,
                first_evalmod_product_scale_log2,
                security_ready ? "true" : "false",
                algebraic_capacity_possible ? "true" : "false",
                scale_strategy_ready ? "true" : "false",
                evalmod_capacity_ready ? "true" : "false",
                feasible ? "true" : "false",
                blocker);
}

} // namespace

int main() {
    std::printf("first_prime_log2,evalmod_prime_log2,evalmod_levels,active_modulus_log2,total_modulus_log2,period_offset_log2,required_period_offset_log2,period_offset_margin_log2,max_target_scale_for_offset_log2,target_scale_log2,min_reserved_modulus_log2,required_drop_log2,available_drop_log2,remaining_after_strategy_log2,first_evalmod_product_scale_log2,security_ready,algebraic_capacity_possible,scale_strategy_ready,evalmod_capacity_ready,feasible,blocker\n");
    for (double first_prime_log2 : {20.0, 30.0, 40.0, 50.0, 60.0}) {
        for (double evalmod_prime_log2 : {2.0, 10.0, 20.0, 30.0, 40.0, 50.0}) {
            for (std::size_t evalmod_levels : {1U, 2U, 3U, 4U}) {
                for (double active_modulus_log2 : {300.0, 500.0, 700.0, 780.0}) {
                    for (double period_offset_log2 : {44.0, 60.0, 80.0, 100.0, 102.0, 120.0, 140.0}) {
                        for (double target_scale_log2 : {30.0, 40.0, 50.0, 60.0}) {
                            run_case({first_prime_log2,
                                      evalmod_prime_log2,
                                      evalmod_levels,
                                      active_modulus_log2,
                                      period_offset_log2,
                                      target_scale_log2});
                        }
                    }
                }
            }
        }
    }
    return 0;
}
