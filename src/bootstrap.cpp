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
        {"ModRaise", BootstrapStageStatus::SpecificationReady, {}, {},
         "подъём ciphertext к расширенной цепочке модулей описан в модели"},
        {"CoeffToSlot", BootstrapStageStatus::PrimitiveReady, {}, {},
         "линейное преобразование опирается на доступные CKKS-ротации"},
        {"EvalMod", BootstrapStageStatus::PrimitiveReady, {}, {},
         "для вычисления доступны multiply, relinearize и rescale"},
        {"SlotToCoeff", BootstrapStageStatus::SpecificationReady, {}, {},
         "обратное линейное преобразование зафиксировано в bootstrapping-модели"}
    };
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
                   BootstrapStageStatus::SpecificationReady,
                   report.depth_boundary,
                   report.depth_boundary,
                   "цель этапа: подготовить ciphertext к расширенной цепочке модулей"),
        make_stage("CoeffToSlot",
                   BootstrapStageStatus::PrimitiveReady,
                   report.depth_boundary,
                   report.depth_boundary,
                   "для линейных преобразований доступны ротации CKKS-слотов"),
        make_stage("EvalMod",
                   BootstrapStageStatus::PrimitiveReady,
                   report.depth_boundary,
                   report.depth_boundary,
                   "для приближённого modular reduction доступны multiply/relinearize/rescale"),
        make_stage("SlotToCoeff",
                   BootstrapStageStatus::SpecificationReady,
                   report.depth_boundary,
                   report.depth_boundary,
                   "этап возвращает данные из slot-представления в coefficient-представление")
    };
    return report;
}

const std::vector<BootstrapStage>& Bootstrapper::pipeline() const noexcept {
    return stages_;
}

} // namespace m2424
