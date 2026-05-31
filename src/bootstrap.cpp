#include "m2424/bootstrap.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <utility>
#include <stdexcept>

namespace m2424 {

const char* to_string(BootstrapStageStatus status) noexcept {
    switch (status) {
        case BootstrapStageStatus::Ready:
            return "ready";
        case BootstrapStageStatus::PrimitiveReady:
            return "primitive_ready";
        case BootstrapStageStatus::SpecificationReady:
            return "specification_ready";
        case BootstrapStageStatus::Blocked:
            return "blocked";
    }
    return "unknown";
}

static BootstrapCipherMetrics capture_metrics(const SealAdapter& adapter, const Cipher& cipher) {
    const auto info = adapter.info(cipher);
    return BootstrapCipherMetrics{
        true,
        info.scale,
        info.chain_index,
        info.coeff_modulus_size,
        info.ciphertext_size,
        adapter.serialized_size(cipher)
    };
}

static BootstrapStage make_stage(std::string name,
                                 BootstrapStageStatus status,
                                 const BootstrapCipherMetrics& before,
                                 const BootstrapCipherMetrics& after,
                                 std::string note) {
    return BootstrapStage{
        std::move(name),
        status,
        before,
        after,
        std::move(note)
    };
}

static bool has_failing_stage(const BootstrapPrototypeReport& report) {
    return std::any_of(report.stages.begin(), report.stages.end(), [](const BootstrapPrototypeStage& stage) {
        return stage.status == "FAIL";
    });
}

static std::string guarded_refresh_blocker(const BootstrapPrototypeReport& report, bool require_value_check) {
    if (report.stages.empty()) {
        return "refresh_report_empty";
    }
    if (has_failing_stage(report)) {
        return "refresh_stage_failed";
    }
    if (!report.restore_level_criterion) {
        return "refresh_did_not_restore_levels";
    }
    if (require_value_check && !report.inside_evalmod_interval) {
        return "refresh_input_outside_evalmod_interval";
    }
    if (require_value_check && !report.preserve_value_criterion) {
        return "refresh_did_not_preserve_value";
    }
    return "none";
}

Bootstrapper::Bootstrapper(SealAdapter& adapter) : adapter_(&adapter) {
    stages_ = {
        {"SlotToCoeff", BootstrapStageStatus::PrimitiveReady, {}, {},
         "рабочий scalable refresh сначала переносит slots в coefficient-представление"},
        {"ModRaise", BootstrapStageStatus::PrimitiveReady, {}, {},
         "ciphertext расширяется к первой RNS-базе modulus chain после SlotToCoeff"},
        {"CoeffToSlot", BootstrapStageStatus::PrimitiveReady, {}, {},
         "обратное FFT-like преобразование возвращает значения в slots перед EvalMod"},
        {"eval_mod_normalization", BootstrapStageStatus::Ready, {}, {},
         "амплитуда входа EvalMod приводится к рабочему интервалу полинома"},
        {"EvalMod", BootstrapStageStatus::PrimitiveReady, {}, {},
         "модульная редукция приближается параметризуемым экспериментальным EvalMod"},
        {"post_refresh_mod_raise", BootstrapStageStatus::PrimitiveReady, {}, {},
         "после refresh ciphertext снова поднимается к первой RNS-базе"}
    };
}

std::vector<int> Bootstrapper::refresh_rotation_steps(std::size_t slots) {
    return BootstrapPrototype::required_rotation_steps(slots);
}

std::vector<int> Bootstrapper::scalable_refresh_rotation_steps(std::size_t slots) {
    auto slot_to_coeff = make_bootstrap_dft_plan(slots, BootstrapDftType::HomomorphicEncode, 40.0);
    auto coeff_to_slot = make_bootstrap_dft_plan(slots, BootstrapDftType::HomomorphicDecode, 40.0);
    auto steps = slot_to_coeff.rotation_steps();
    const auto inverse_steps = coeff_to_slot.rotation_steps();
    steps.insert(steps.end(), inverse_steps.begin(), inverse_steps.end());
    std::sort(steps.begin(), steps.end());
    steps.erase(std::unique(steps.begin(), steps.end()), steps.end());
    return steps;
}

BootstrapReport Bootstrapper::analyze_depth(const std::vector<double>& input, std::size_t max_steps) {
    if (!adapter_) {
        throw std::runtime_error("Bootstrapper has no SealAdapter");
    }
    if (input.empty()) {
        throw std::invalid_argument("input must not be empty");
    }

    auto current = adapter_->encrypt(adapter_->encode(input));
    BootstrapReport report;
    report.input = capture_metrics(*adapter_, current);
    report.depth_boundary = report.input;
    report.next_exponent = 1;

    for (std::size_t step = 0; step < max_steps; ++step) {
        try {
            current = adapter_->mul_relin_rescale(current, current);
            ++report.successful_multiplications;
            report.next_exponent *= 2;
            report.depth_boundary = capture_metrics(*adapter_, current);
        } catch (const std::exception& error) {
            report.stop_reason = error.what();
            break;
        }
    }

    if (report.stop_reason.empty()) {
        report.stop_reason = "max_steps reached before depth failure";
    }

    const bool has_depth_boundary = report.depth_boundary.available;
    const bool level_restored = has_depth_boundary && report.depth_boundary.chain_index > report.input.chain_index;
    report.preserve_value_criterion = false;
    report.restore_level_criterion = level_restored;

    report.stages = {
        make_stage("ModRaise",
                   BootstrapStageStatus::PrimitiveReady,
                   report.depth_boundary,
                   report.depth_boundary,
                   "цель этапа: подготовить ciphertext к расширенной цепочке модулей"),
        make_stage("CoeffToSlot",
                   BootstrapStageStatus::PrimitiveReady,
                   report.depth_boundary,
                   report.depth_boundary,
                   "для линейных преобразований доступны ротации CKKS-слотов"),
        make_stage("eval_mod_normalization",
                   BootstrapStageStatus::Ready,
                   report.depth_boundary,
                   report.depth_boundary,
                   "коэффициент нормализации выбирается по амплитуде после CoeffToSlot"),
        make_stage("EvalMod",
                   BootstrapStageStatus::PrimitiveReady,
                   report.depth_boundary,
                   report.depth_boundary,
                   "для приближённого modular reduction доступны multiply/relinearize/rescale"),
        make_stage("SlotToCoeff",
                   BootstrapStageStatus::PrimitiveReady,
                   report.depth_boundary,
                   report.depth_boundary,
                   "этап возвращает данные из slot-представления в coefficient-представление"),
        make_stage("post_refresh_mod_raise",
                   BootstrapStageStatus::PrimitiveReady,
                   report.depth_boundary,
                   report.depth_boundary,
                   "результат refresh поднимается для продолжения вычислений")
    };
    return report;
}

BootstrapPrototypeReport Bootstrapper::refresh(const Cipher& input, std::size_t slots, double tolerance) {
    return refresh(input, slots, tolerance, 1.0);
}

BootstrapPrototypeReport Bootstrapper::refresh(const Cipher& input,
                                               std::size_t slots,
                                               double tolerance,
                                               double normalization_factor) {
    if (!adapter_) {
        throw std::runtime_error("Bootstrapper has no SealAdapter");
    }
    BootstrapPrototype prototype(*adapter_, slots, tolerance, normalization_factor);
    return prototype.refresh_cipher_fast(input);
}

BootstrapPrototypeReport Bootstrapper::refresh_checked(const Cipher& input,
                                                       const ComplexVector& expected,
                                                       std::size_t slots,
                                                       double tolerance) {
    return refresh_checked(input, expected, slots, tolerance, 1.0);
}

BootstrapPrototypeReport Bootstrapper::refresh_checked(const Cipher& input,
                                                       const ComplexVector& expected,
                                                       std::size_t slots,
                                                       double tolerance,
                                                       double normalization_factor) {
    if (!adapter_) {
        throw std::runtime_error("Bootstrapper has no SealAdapter");
    }
    BootstrapPrototype prototype(*adapter_, slots, tolerance, normalization_factor);
    return prototype.refresh_cipher_checked(input, expected);
}

BootstrapPrototypeReport Bootstrapper::refresh_slots_to_coeffs_first(const Cipher& input,
                                                                     std::size_t slots,
                                                                     double tolerance) {
    return refresh_slots_to_coeffs_first(input, slots, tolerance, EvalModDegree::P3);
}

BootstrapPrototypeReport Bootstrapper::refresh_slots_to_coeffs_first(const Cipher& input,
                                                                     std::size_t slots,
                                                                     double tolerance,
                                                                     EvalModDegree evalmod_degree) {
    if (!adapter_) {
        throw std::runtime_error("Bootstrapper has no SealAdapter");
    }
    BootstrapPrototype prototype(*adapter_, slots, tolerance);
    prototype.set_transform_backend(BootstrapTransformBackend::FftLike);
    prototype.set_circuit_order(BootstrapCircuitOrder::SlotsToCoeffsFirst);
    prototype.set_evalmod_degree(evalmod_degree);
    prototype.set_plain_scale_log2(std::log2(adapter_->info(input).scale));
    return prototype.refresh_cipher_fast(input);
}

BootstrapPrototypeReport Bootstrapper::refresh_slots_to_coeffs_first_checked(
    const Cipher& input,
    const ComplexVector& expected,
    std::size_t slots,
    double tolerance) {
    return refresh_slots_to_coeffs_first_checked(input, expected, slots, tolerance, EvalModDegree::P3);
}

BootstrapPrototypeReport Bootstrapper::refresh_slots_to_coeffs_first_checked(
    const Cipher& input,
    const ComplexVector& expected,
    std::size_t slots,
    double tolerance,
    EvalModDegree evalmod_degree) {
    if (!adapter_) {
        throw std::runtime_error("Bootstrapper has no SealAdapter");
    }
    BootstrapPrototype prototype(*adapter_, slots, tolerance);
    prototype.set_transform_backend(BootstrapTransformBackend::FftLike);
    prototype.set_circuit_order(BootstrapCircuitOrder::SlotsToCoeffsFirst);
    prototype.set_evalmod_degree(evalmod_degree);
    prototype.set_plain_scale_log2(std::log2(adapter_->info(input).scale));
    return prototype.refresh_cipher_checked(input, expected);
}

BootstrapGuardedRefreshResult Bootstrapper::refresh_slots_to_coeffs_first_guarded(
    const Cipher& input,
    const CkksOperationBudget& operation_budget,
    double target_error,
    std::size_t slots,
    double tolerance,
    int security_bits,
    ParameterOptimizeFor optimize_for,
    std::size_t min_chain_remaining_after_compute) {
    return refresh_slots_to_coeffs_first_guarded(input,
                                                 operation_budget,
                                                 target_error,
                                                 slots,
                                                 tolerance,
                                                 EvalModDegree::P3,
                                                 security_bits,
                                                 optimize_for,
                                                 min_chain_remaining_after_compute);
}

BootstrapGuardedRefreshResult Bootstrapper::refresh_slots_to_coeffs_first_guarded(
    const Cipher& input,
    const CkksOperationBudget& operation_budget,
    double target_error,
    std::size_t slots,
    double tolerance,
    EvalModDegree evalmod_degree,
    int security_bits,
    ParameterOptimizeFor optimize_for,
    std::size_t min_chain_remaining_after_compute) {
    BootstrapGuardedRefreshResult result;
    result.planning = plan_refresh_for_budget(input,
                                              operation_budget,
                                              target_error,
                                              slots,
                                              security_bits,
                                              optimize_for,
                                              min_chain_remaining_after_compute);
    if (result.planning.status != BootstrapRefreshPlanningStatus::RefreshRequired) {
        result.blocker = result.planning.blocker;
        return result;
    }
    result.refresh = refresh_slots_to_coeffs_first(input, slots, tolerance, evalmod_degree);
    result.refresh_executed = true;
    result.blocker = guarded_refresh_blocker(result.refresh, false);
    return result;
}

BootstrapGuardedRefreshResult Bootstrapper::refresh_slots_to_coeffs_first_checked_guarded(
    const Cipher& input,
    const ComplexVector& expected,
    const CkksOperationBudget& operation_budget,
    double target_error,
    std::size_t slots,
    double tolerance,
    int security_bits,
    ParameterOptimizeFor optimize_for,
    std::size_t min_chain_remaining_after_compute) {
    return refresh_slots_to_coeffs_first_checked_guarded(input,
                                                         expected,
                                                         operation_budget,
                                                         target_error,
                                                         slots,
                                                         tolerance,
                                                         EvalModDegree::P3,
                                                         security_bits,
                                                         optimize_for,
                                                         min_chain_remaining_after_compute);
}

BootstrapGuardedRefreshResult Bootstrapper::refresh_slots_to_coeffs_first_checked_guarded(
    const Cipher& input,
    const ComplexVector& expected,
    const CkksOperationBudget& operation_budget,
    double target_error,
    std::size_t slots,
    double tolerance,
    EvalModDegree evalmod_degree,
    int security_bits,
    ParameterOptimizeFor optimize_for,
    std::size_t min_chain_remaining_after_compute) {
    BootstrapGuardedRefreshResult result;
    result.planning = plan_refresh_for_budget(input,
                                              operation_budget,
                                              target_error,
                                              slots,
                                              security_bits,
                                              optimize_for,
                                              min_chain_remaining_after_compute);
    if (result.planning.status != BootstrapRefreshPlanningStatus::RefreshRequired) {
        result.blocker = result.planning.blocker;
        return result;
    }
    result.refresh = refresh_slots_to_coeffs_first_checked(input, expected, slots, tolerance, evalmod_degree);
    result.refresh_executed = true;
    result.blocker = guarded_refresh_blocker(result.refresh, true);
    return result;
}

BootstrapRefreshPlanningResult Bootstrapper::plan_refresh_for_budget(
    const Cipher& input,
    const CkksOperationBudget& operation_budget,
    double target_error,
    std::size_t slots,
    int security_bits,
    ParameterOptimizeFor optimize_for,
    std::size_t min_chain_remaining_after_compute) const {
    if (!adapter_) {
        throw std::runtime_error("Bootstrapper has no SealAdapter");
    }
    return plan_bootstrap_refresh({
        adapter_->info(input),
        operation_budget,
        target_error,
        slots,
        security_bits,
        optimize_for,
        min_chain_remaining_after_compute
    });
}

const std::vector<BootstrapStage>& Bootstrapper::pipeline() const noexcept {
    return stages_;
}

BootstrapPipelinePlan Bootstrapper::plan(std::size_t slots) const {
    return make_scalable_bootstrap_plan(slots);
}

} // namespace m2424
