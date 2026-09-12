#pragma once
#include "m2424/bootstrap.hpp"
#include "certified_arithmetic_internal.hpp"
namespace m2424 {
struct BootstrapPlan::Data {
    BootstrapTrace trace;
    BootstrapContractResult readiness{BootstrapCertificationStatus::RequiredBoundUnavailable,"plan","No executable plan"};
    BootstrapRequest request;
    std::optional<CertifiedEvalRoundPlusCoeffToSlot> cts;
    EvalRoundExecutionPlan evalRound;
    PreparedSlotToCoeffPlan stc;
    std::vector<EvalRoundExecutionNode> combination;
    std::vector<Plain> constants;
    std::size_t combinationOutput{};
    std::size_t degree{},ctsDepth{},stcDepth{};
    std::vector<int> keys;
};
}
