#pragma once
#include <cmath>
#include "m2424/sparse_bootstrap.hpp"
#include "m2424/evalround_plus_coeff_to_slot.hpp"
#include "m2424/evalround_execution.hpp"
namespace m2424 {
enum class BootstrapLiftEvidence { Unknown, Analytical, TestFixtureAssumption };
struct BootstrapLiftBound {
    std::optional<std::uint32_t> K;
    BootstrapLiftEvidence evidence{BootstrapLiftEvidence::Unknown};
    std::string provenance;
};
struct BootstrapRequest {
    BootstrapTarget target;
    CoeffToSlotCertificationInput upstream;
    // Legacy test fixture only. Production lift comes from generated sparse keys.
    BootstrapLiftBound lift;
    std::optional<BootstrapBound> sparseEvaluationKeyNoiseSupport;
    std::vector<RlweSecurityEvidence> securityEvidence;
    // Fixed output-side schedule. No scale metadata rewrite is permitted.
    double combinationScale{std::ldexp(1.,49)};
    double minimumArithmeticScale{std::ldexp(1.,40)};
    double minimumConstantScale{std::ldexp(1.,40)};
    std::size_t maxIntegerRescales{3};
    std::size_t maxCleaningRounds{4};
    std::size_t outputChainIndex{2};
};
class BootstrapPlan {
public:
    BootstrapPlan();
    const BootstrapTrace& trace() const;
    const BootstrapContractResult& executionReadiness() const;
    const EvalRoundExecutionPlan& evalRound() const;
    const CertifiedEvalRoundPlusCoeffToSlot* coeffToSlot() const;
    const std::vector<EvalRoundExecutionNode>& combinationNodes() const;
private:
    struct Data; std::shared_ptr<const Data> data_;
    friend class Bootstrapper;
};
struct BootstrapResult { std::optional<Cipher> output; BootstrapTrace trace; };
using BootstrapObserver=std::function<void(const std::string&,std::size_t,std::size_t,const Cipher&,const BootstrapBound&)>;
class Bootstrapper {
public:
    explicit Bootstrapper(std::size_t degree,std::size_t coeffToSlotDepth=4,std::size_t slotToCoeffDepth=4);
    std::vector<int> rotationKeys() const;
    // Optional analysis prepares every executable stage; no secret-key oracle.
    BootstrapPlan prepare(SealAdapter&,const Cipher&,const BootstrapRequest&) const;
    BootstrapContractResult preflight(const SealAdapter&,const Cipher&,const BootstrapPlan&) const;
    BootstrapResult apply(SealAdapter&,const Cipher&,const BootstrapPlan&,const BootstrapObserver& = {}) const;
private:
    EvalRoundPlusCoeffToSlot cts_;
    SlotToCoeffPlan stc_;
};
}
