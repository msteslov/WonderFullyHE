#pragma once
#include "m2424/evalround_execution.hpp"
namespace m2424 {
struct EvalRoundExecutionPlan::Data {
    BootstrapContractResult certification{BootstrapCertificationStatus::RequiredBoundUnavailable,
        "EvalRoundExecution", "No compiled arithmetic certificate"};
    std::array<std::uint64_t,4> fingerprint{};
    std::vector<EvalRoundExecutionNode> nodes;
    std::vector<Plain> constants;
    EvalRoundPlan mathematicalPlan;
    EvalRoundExecutionDiagnostics diagnostics;
    std::size_t output{};
};
}
