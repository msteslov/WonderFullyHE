#include "m2424/bootstrap.hpp"

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

Bootstrapper::Bootstrapper(SealAdapter& adapter) : adapter_(&adapter) {
    stages_ = {
        {"ModRaise", BootstrapStageStatus::PrimitiveReady, {}, {},
         "ciphertext расширяется к первой RNS-базе modulus chain"},
        {"CoeffToSlot", BootstrapStageStatus::PrimitiveReady, {}, {},
         "линейное преобразование выполняется через CKKS-ротации и plaintext-диагонали"},
        {"eval_mod_normalization", BootstrapStageStatus::Ready, {}, {},
         "амплитуда входа EvalMod приводится к рабочему интервалу полинома"},
        {"EvalMod", BootstrapStageStatus::PrimitiveReady, {}, {},
         "модульная редукция приближается полиномом степени 7"},
        {"SlotToCoeff", BootstrapStageStatus::PrimitiveReady, {}, {},
         "обратное линейное преобразование возвращает coefficient-представление"},
        {"post_refresh_mod_raise", BootstrapStageStatus::PrimitiveReady, {}, {},
         "после refresh ciphertext снова поднимается к первой RNS-базе"}
    };
}

std::vector<int> Bootstrapper::refresh_rotation_steps(std::size_t slots) {
    return BootstrapPrototype::required_rotation_steps(slots);
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

const std::vector<BootstrapStage>& Bootstrapper::pipeline() const noexcept {
    return stages_;
}

BootstrapPipelinePlan Bootstrapper::plan(std::size_t slots) const {
    return make_research_bootstrap_plan(slots);
}

} // namespace m2424
