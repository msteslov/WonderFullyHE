#pragma once

#include "m2424/coeff_to_slot.hpp"
#include "m2424/slot_to_coeff.hpp"
#include "m2424/sparse_bootstrap.hpp"

namespace m2424 {

struct CoeffToSlotFactorTrace {
    std::string operation;
    CipherInfo outputState;
    std::vector<std::uint64_t> activePrimes;
    BootstrapBound kappa;
    BootstrapBound diagonalError;
    BootstrapBound idealMagnitude;
    BootstrapBound propagatedSemanticError;
    BootstrapBound localAddedError;
};

struct CoeffToSlotBranchTrace {
    BootstrapGate gate{BootstrapGate::CoeffToSlotHP};
    CoeffToSlotPrefactor prefactor;
    BootstrapBound inputMagnitude;
    BootstrapBound inputSemanticError;
    std::optional<SlotToCoeffRuntimeStage> firstFactorRestoration;
    std::array<std::vector<CoeffToSlotFactorTrace>, 2> halves;
    std::size_t levelsConsumed{};
    std::size_t rescaleOperations{};
    BootstrapBound outputError;
    BootstrapGateEvidence evidence;
    BootstrapContractResult certificate;
    std::vector<SlotToCoeffFactorTrace> certifiedFactors;
    BootstrapContractResult errorBudget;
};

struct EvalRoundPlusCoeffToSlotResult {
    Cipher hpFirst, hpSecond, lpFirst, lpSecond;
    CoeffToSlotBranchTrace hpTrace, lpTrace;
    std::size_t restorationOperations{};
};


struct CoeffToSlotCertificationInput {
    // Sup norm of the exact raised y=sigma(u/Delta0), including source noise.
    BootstrapBound raisedMagnitude;
    // Coefficient sup norms in unscaled integer units, before ModRaise.
    BootstrapBound messageMagnitude, sourceNoiseMagnitude;
};
struct CoeffToSlotDomainCertificate {
    BootstrapBound rho;
    BootstrapContractResult result;
};
// Optional analysis: unscaled coefficient bounds and the certified LP error.
CoeffToSlotDomainCertificate certifyCoeffToSlotDomain(const BootstrapInputContext&,
    const CoeffToSlotCertificationInput&, const BootstrapBound& lpError);
class CertifiedEvalRoundPlusCoeffToSlot {
public:
    const CoeffToSlotBranchTrace& hp() const;
    const CoeffToSlotBranchTrace& lp() const;
    const CoeffToSlotDomainCertificate& domain() const;
private:
    struct Impl; std::shared_ptr<const Impl> impl_;
    friend class EvalRoundPlusCoeffToSlot;
};
struct RootLinearTransformPlan;

class PreparedEvalRoundPlusCoeffToSlot {
public:
    PreparedEvalRoundPlusCoeffToSlot(PreparedEvalRoundPlusCoeffToSlot&&) noexcept = default;
    PreparedEvalRoundPlusCoeffToSlot& operator=(PreparedEvalRoundPlusCoeffToSlot&&) noexcept = default;
    PreparedEvalRoundPlusCoeffToSlot(const PreparedEvalRoundPlusCoeffToSlot&) = delete;
    PreparedEvalRoundPlusCoeffToSlot& operator=(const PreparedEvalRoundPlusCoeffToSlot&) = delete;
private:
    PreparedEvalRoundPlusCoeffToSlot(PreparedCoeffToSlotPlan, PreparedCoeffToSlotPlan,
        BootstrapInputContext, CoeffToSlotContract, CoeffToSlotContract);
    PreparedCoeffToSlotPlan hp_, lp_;
    BootstrapInputContext source_;
    CoeffToSlotContract hpContract_, lpContract_;
    friend class EvalRoundPlusCoeffToSlot;
};

/// PR-1 branch executor. Never constructs or certifies a full bootstrap plan.
/// Source context must be resolved on the input BEFORE ModRaise.
class EvalRoundPlusCoeffToSlot {
public:
    explicit EvalRoundPlusCoeffToSlot(std::size_t degree, std::size_t depth = 4);
    const CoeffToSlotPlan& plan() const noexcept { return plan_; }
    CoeffToSlotPlanRequirements requirements() const { return plan_.requirements(); }
    BootstrapContractResult preflight(const SealAdapter&, const RaisedCipher&,
        const BootstrapInputContext&, const CoeffToSlotContract& hp,
        const CoeffToSlotContract& lp) const;
    PreparedEvalRoundPlusCoeffToSlot prepare(SealAdapter&, const RaisedCipher&,
        const BootstrapInputContext&, const CoeffToSlotContract& hp,
        const CoeffToSlotContract& lp) const;
    /// Executes independently from the SAME raised state; does not consume input.
    /// inputMagnitude bounds |y|, not decrypted measurements. Unknown is allowed
    /// for diagnostics, but remains Unknown in the trace. nu_b is part of u.
    EvalRoundPlusCoeffToSlotResult apply(SealAdapter&, const RaisedCipher&,
        const PreparedEvalRoundPlusCoeffToSlot&, const BootstrapBound& inputMagnitude = {}) const;
    // Optional analysis preparation; execution remains in the dependency-free core.
    CertifiedEvalRoundPlusCoeffToSlot prepareCertified(SealAdapter&, const RaisedCipher&,
        const BootstrapInputContext&, const CoeffToSlotContract&, const CoeffToSlotContract&,
        const CoeffToSlotCertificationInput&) const;
    BootstrapContractResult preflight(const SealAdapter&, const RaisedCipher&,
        const BootstrapInputContext&, const CertifiedEvalRoundPlusCoeffToSlot&) const;
    EvalRoundPlusCoeffToSlotResult apply(SealAdapter&, const RaisedCipher&,
        const CertifiedEvalRoundPlusCoeffToSlot&,
        const std::function<void(BootstrapGate,const SlotToCoeffRuntimeStage&,const Cipher&)>& observer={}) const;
    // Sparse restoration is owned by this first-factor path and executed once,
    // shared by HP/LP. Preparation only examines sparse state metadata.
    CertifiedEvalRoundPlusCoeffToSlot prepareCertified(SealAdapter&,const SparseRaisedCipher&,
        const SparseBootstrapPlan&,const CoeffToSlotContract&,const CoeffToSlotContract&) const;
    BootstrapContractResult preflight(const SealAdapter&,const SparseRaisedCipher&,
        const CertifiedEvalRoundPlusCoeffToSlot&) const;
    EvalRoundPlusCoeffToSlotResult apply(SealAdapter&,const SparseRaisedCipher&,
        const CertifiedEvalRoundPlusCoeffToSlot&,
        const std::function<void(BootstrapGate,const SlotToCoeffRuntimeStage&,const Cipher&)>& observer={}) const;
    std::vector<ComplexVector> applyPlainTrace(const ComplexVector&, std::size_t half,
        const CoeffToSlotPrefactor& = {}) const;
private:
    CertifiedEvalRoundPlusCoeffToSlot prepareCertifiedImpl(SealAdapter&,const RaisedCipher&,
        const BootstrapInputContext&,const CoeffToSlotContract&,const CoeffToSlotContract&,
        const CoeffToSlotCertificationInput&,const SparseBootstrapPlan*) const;
    BootstrapContractResult preflightCertifiedState(const SealAdapter&,const RaisedCipher&,
        const BootstrapInputContext&,const CertifiedEvalRoundPlusCoeffToSlot&) const;
    EvalRoundPlusCoeffToSlotResult applyCertifiedFactors(SealAdapter&,const RaisedCipher&,
        const CertifiedEvalRoundPlusCoeffToSlot&,
        const std::function<void(BootstrapGate,const SlotToCoeffRuntimeStage&,const Cipher&)>&) const;
    std::vector<std::size_t> certificationBabySteps() const;
    RootLinearTransformPlan certificationLayout(const CoeffToSlotPrefactor&) const;
    CoeffToSlotPlan plan_;
};

} // namespace m2424
