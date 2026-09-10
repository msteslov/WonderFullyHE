#pragma once

#include "m2424/coeff_to_slot.hpp"

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
    std::array<std::vector<CoeffToSlotFactorTrace>, 2> halves;
    std::size_t levelsConsumed{};
    std::size_t rescaleOperations{};
    BootstrapBound outputError;
    BootstrapGateEvidence evidence;
    BootstrapContractResult certificate;
};

struct EvalRoundPlusCoeffToSlotResult {
    Cipher hpFirst, hpSecond, lpFirst, lpSecond;
    CoeffToSlotBranchTrace hpTrace, lpTrace;
};

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
private:
    CoeffToSlotPlan plan_;
};

} // namespace m2424
