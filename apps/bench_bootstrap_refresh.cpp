#include "m2424/bootstrap.hpp"
#include "m2424/diagonal_transform.hpp"
#include "m2424/eval_mod.hpp"
#include "m2424/profiles.hpp"
#include "m2424/seal_adapter.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <exception>
#include <string>
#include <vector>

namespace {

template <typename Fn>
double elapsed_ms(Fn&& fn) {
    const auto started = std::chrono::steady_clock::now();
    fn();
    const auto finished = std::chrono::steady_clock::now();
    return std::chrono::duration<double, std::milli>(finished - started).count();
}

std::vector<double> make_input(std::size_t slots, double amplitude) {
    std::vector<double> input;
    input.reserve(slots);
    for (std::size_t i = 0; i < slots; ++i) {
        const double x = static_cast<double>(i);
        input.push_back(amplitude * std::sin(x / 4.0));
    }
    return input;
}

m2424::ComplexVector head(m2424::ComplexVector values, std::size_t n) {
    values.resize(std::min(values.size(), n));
    return values;
}

double max_abs_value(const m2424::ComplexVector& values) {
    double result = 0.0;
    for (const auto& value : values) {
        result = std::max(result, std::abs(value));
    }
    return result;
}

std::string csv_safe(std::string value) {
    std::replace(value.begin(), value.end(), ',', ';');
    return value;
}

} // namespace

int main() {
    constexpr std::size_t slots = 16;
    constexpr double tolerance = 2e-5;
    constexpr double amplitude = 1e-5;
    const auto profile = m2424::profiles::boot_ckks();

    std::vector<int> rotation_steps;
    const double plan_ms = elapsed_ms([&] {
        rotation_steps = m2424::Bootstrapper::refresh_rotation_steps(slots);
    });

    auto adapter = m2424::SealAdapter::create(profile);
    const double keygen_ms = elapsed_ms([&] {
        adapter.keygen(rotation_steps, true);
    });

    const auto input = make_input(slots, amplitude);
    auto encrypted = adapter.encrypt(adapter.encode(input));
    auto lowered = adapter.mul_plain_rescale(encrypted, adapter.encode_scalar_like(1.0, encrypted));
    const auto before = adapter.info(lowered);

    m2424::Bootstrapper bootstrapper(adapter);

    m2424::CkksOperationBudget continuation_budget;
    continuation_budget.ciphertext_muls = before.chain_index + 1;
    const auto refresh_plan = bootstrapper.plan_refresh_for_budget(lowered,
                                                                   continuation_budget,
                                                                   1e-9,
                                                                   adapter.slot_count());

    std::printf("planner,status,required_levels,available_levels,needs_refresh,planned_work_bits,estimated_abs_error_bound,blocker\n");
    std::printf("refresh_gate,%s,%zu,%zu,%s,%d,%.6e,%s\n",
                m2424::to_string(refresh_plan.status),
                refresh_plan.required_compute_levels,
                refresh_plan.available_compute_levels,
                refresh_plan.needs_refresh ? "true" : "false",
                refresh_plan.parameter_plan.selected_work_bits,
                refresh_plan.parameter_plan.estimated_abs_error_bound,
                refresh_plan.blocker.c_str());

    if (refresh_plan.status != m2424::BootstrapRefreshPlanningStatus::RefreshRequired) {
        std::printf("summary,profile,slots,tolerance,rotation_keys,plan_ms,keygen_ms,refresh_ms,chain_before,chain_after,scale_before,scale_after,normalization_factor,status\n");
        std::printf("refresh,boot_ckks,%zu,%.6e,%zu,%.6f,%.6f,%.6f,%zu,%zu,%.6e,%.6e,%.6e,%s\n",
                    slots,
                    tolerance,
                    rotation_steps.size(),
                    plan_ms,
                    keygen_ms,
                    0.0,
                    before.chain_index,
                    before.chain_index,
                    before.scale,
                    before.scale,
                    1.0,
                    m2424::to_string(refresh_plan.status));
        return 0;
    }

    try {
        auto raised = adapter.mod_raise_to_first(lowered);
        const auto after_mod_raise = adapter.info(raised);
        auto coeff_to_slot = m2424::DiagonalLinearTransform::from_matrix(
            m2424::canonical_embedding_matrix(slots));
        auto slot_domain = coeff_to_slot.apply(adapter, raised);
        const auto slot_domain_info = adapter.info(slot_domain);
        const double max_abs_before_normalization =
            max_abs_value(head(adapter.decode_complex(adapter.decrypt(slot_domain)), slots));
        const auto scale_search = m2424::search_bootstrap_refresh_scale_gate({
            profile,
            before,
            after_mod_raise,
            slot_domain_info,
            max_abs_before_normalization,
            {
                m2424::BootstrapPeriodMode::TotalCoeffModulus,
                m2424::BootstrapPeriodMode::SourceCoeffModulus,
                m2424::BootstrapPeriodMode::DroppedPrimeProduct,
                m2424::BootstrapPeriodMode::LastPrime,
                m2424::BootstrapPeriodMode::ManualPowerOfTwo
            },
            {200.0, 220.0, 240.0, 260.0, 280.0, 300.0, 320.0, 340.0},
            {0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0},
            {40.0, 60.0, 80.0, 120.0, 160.0},
            2.0
        });
        const auto& scale_design = scale_search.best_design;

        std::printf("scale_gate,status,period_log2,coeff_to_slot_prescale_log2,coeff_to_slot_plain_scale_log2,plain_scale_log2,target_scale_log2,max_abs_before_normalization,required_levels,chain_remaining_after_strategy,missing_drop_log2,missing_scalar_levels,missing_total_levels,candidates,ready_candidates,blocker\n");
        std::printf("experimental_refresh,%s,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%zu,%zu,%.6e,%zu,%zu,%zu,%zu,%s\n",
                    m2424::to_string(scale_design.status),
                    scale_design.period_log2,
                    scale_design.coeff_to_slot_prescale_log2,
                    scale_design.coeff_to_slot_plain_scale_log2,
                    scale_design.plain_scale_log2,
                    scale_design.target_scale_log2,
                    max_abs_before_normalization,
                    scale_design.required_levels,
                    scale_design.chain_remaining_after_strategy,
                    scale_design.scale_plan.missing_drop_log2,
                    scale_design.scale_plan.missing_scalar_levels,
                    scale_design.scale_plan.missing_total_levels,
                    scale_search.candidates,
                    scale_search.ready_candidates,
                    csv_safe(scale_design.blocker).c_str());

        if (scale_design.status != m2424::BootstrapScaleDesignStatus::ReadyForEvalModP3
            || scale_design.coeff_to_slot_prescale_log2 != 0.0) {
            std::printf("summary,profile,slots,tolerance,rotation_keys,plan_ms,keygen_ms,refresh_ms,chain_before,chain_after,scale_before,scale_after,normalization_factor,status\n");
            std::printf("refresh,boot_ckks,%zu,%.6e,%zu,%.6f,%.6f,%.6f,%zu,%zu,%.6e,%.6e,%.6e,%s:%s\n",
                        slots,
                        tolerance,
                        rotation_steps.size(),
                        plan_ms,
                        keygen_ms,
                        0.0,
                        before.chain_index,
                        before.chain_index,
                        before.scale,
                        before.scale,
                        1.0,
                        scale_design.coeff_to_slot_prescale_log2 == 0.0
                            ? "blocked_by_scale_gate"
                            : "blocked_by_unapplied_prescale",
                        m2424::to_string(scale_design.status));
            return 0;
        }
    } catch (const std::exception& error) {
        std::printf("scale_gate,status,period_log2,coeff_to_slot_prescale_log2,coeff_to_slot_plain_scale_log2,plain_scale_log2,target_scale_log2,max_abs_before_normalization,required_levels,chain_remaining_after_strategy,missing_drop_log2,missing_scalar_levels,missing_total_levels,candidates,ready_candidates,blocker\n");
        std::printf("experimental_refresh,preflight_exception,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0,0,0.000000e+00,0,0,0,0,%s\n",
                    csv_safe(error.what()).c_str());
        std::printf("summary,profile,slots,tolerance,rotation_keys,plan_ms,keygen_ms,refresh_ms,chain_before,chain_after,scale_before,scale_after,normalization_factor,status\n");
        std::printf("refresh,boot_ckks,%zu,%.6e,%zu,%.6f,%.6f,%.6f,%zu,%zu,%.6e,%.6e,%.6e,blocked_by_scale_gate_exception\n",
                    slots,
                    tolerance,
                    rotation_steps.size(),
                    plan_ms,
                    keygen_ms,
                    0.0,
                    before.chain_index,
                    before.chain_index,
                    before.scale,
                    before.scale,
                    1.0);
        return 0;
    }

    m2424::BootstrapPrototypeReport report;
    double refresh_ms{};
    try {
        refresh_ms = elapsed_ms([&] {
            report = bootstrapper.refresh(lowered, slots, tolerance);
        });
    } catch (const std::exception& error) {
        std::printf("summary,profile,slots,tolerance,rotation_keys,plan_ms,keygen_ms,refresh_ms,chain_before,chain_after,scale_before,scale_after,normalization_factor,status\n");
        std::printf("refresh,boot_ckks,%zu,%.6e,%zu,%.6f,%.6f,%.6f,%zu,%zu,%.6e,%.6e,%.6e,exception:%s\n",
                    slots,
                    tolerance,
                    rotation_steps.size(),
                    plan_ms,
                    keygen_ms,
                    0.0,
                    before.chain_index,
                    before.chain_index,
                    before.scale,
                    before.scale,
                    1.0,
                    error.what());
        return 1;
    }

    const auto after = adapter.info(report.result);

    std::printf("summary,profile,slots,tolerance,rotation_keys,plan_ms,keygen_ms,refresh_ms,chain_before,chain_after,scale_before,scale_after,normalization_factor,status\n");
    std::printf("refresh,boot_ckks,%zu,%.6e,%zu,%.6f,%.6f,%.6f,%zu,%zu,%.6e,%.6e,%.6e,%s\n",
                report.slots,
                report.tolerance,
                rotation_steps.size(),
                plan_ms,
                keygen_ms,
                refresh_ms,
                before.chain_index,
                after.chain_index,
                before.scale,
                after.scale,
                report.normalization_factor,
                report.restore_level_criterion ? "PASS" : "FAIL");

    std::printf("stage,status,chain_before,chain_after,scale_before,scale_after,max_abs_error,duration_ms\n");
    for (const auto& stage : report.stages) {
        std::printf("%s,%s,%zu,%zu,%.6e,%.6e,%.6e,%.6f\n",
                    stage.name.c_str(),
                    stage.status.c_str(),
                    stage.chain_before,
                    stage.chain_after,
                    stage.scale_before,
                    stage.scale_after,
                    stage.max_abs_error,
                    stage.duration_ms);
    }

    return report.restore_level_criterion ? 0 : 1;
}
