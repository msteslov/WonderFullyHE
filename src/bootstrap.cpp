#include "m2424/bootstrap.hpp"

#include <exception>
#include <stdexcept>

namespace m2424 {

const char* to_string(BootstrapStageStatus status) noexcept {
    switch (status) {
        case BootstrapStageStatus::Done:
            return "done";
        case BootstrapStageStatus::Prototype:
            return "prototype";
        case BootstrapStageStatus::ModelOnly:
            return "model";
        case BootstrapStageStatus::Planned:
            return "planned";
        case BootstrapStageStatus::Failed:
            return "failed";
    }
    return "unknown";
}

Bootstrapper::Bootstrapper(SealAdapter& adapter) : adapter_(&adapter) {
    stages_ = {
        {"Depth diagnostics", BootstrapStageStatus::Done,
         "detects the multiplication depth limit that bootstrap must refresh"},
        {"ModRaise", BootstrapStageStatus::ModelOnly,
         "mathematical stage is fixed; SEAL-level implementation is in progress"},
        {"CoeffToSlot", BootstrapStageStatus::Prototype,
         "rotation primitives required for linear transforms are available"},
        {"EvalMod", BootstrapStageStatus::Prototype,
         "homomorphic multiplication/rescale building block is available"},
        {"SlotToCoeff", BootstrapStageStatus::ModelOnly,
         "inverse linear transform is described in the bootstrap model"},
        {"Refresh ciphertext", BootstrapStageStatus::Planned,
         "full refreshed ciphertext output will replace the diagnostic prototype"}
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
    report.next_exponent = 1;

    for (std::size_t step = 0; step < max_steps; ++step) {
        try {
            current = adapter_->mul_relin_rescale(current, current);
            ++report.successful_multiplications;
            report.next_exponent *= 2;
        } catch (const std::exception& error) {
            report.stop_reason = error.what();
            break;
        }
    }

    if (report.stop_reason.empty()) {
        report.stop_reason = "max_steps reached before depth failure";
    }

    report.stages = stages_;
    return report;
}

const std::vector<BootstrapStage>& Bootstrapper::pipeline() const noexcept {
    return stages_;
}

} // namespace m2424
