#pragma once
#include "m2424/bootstrap_contract.hpp"
namespace m2424 {
struct SparseBootstrapInput {
    BootstrapBound messageMagnitude,sourceNoiseMagnitude;
    // Unknown by default: preparation resolves the configured backend support.
    // A supplied Unknown override deliberately disables certification.
    std::optional<BootstrapBound> evaluationKeyNoiseSupport;
    std::vector<RlweSecurityEvidence> securityEvidence;
    int targetSecurityBits{128};
};
struct SparseBootstrapCertificate {
    BootstrapInputContext input;
    SparseKeyMetadata keys;
    std::uint32_t K{};
    std::string liftProvenance;
    BootstrapBound encapsulationError,restorationError,sourceNoiseMagnitude;
    BootstrapBound raisedMagnitude,restoredMagnitude,rhoBeforeCoeffToSlot;
    BootstrapSecurityReport security;
    BootstrapContractResult result;
};
class SparseBootstrapPlan {
public:
    SparseBootstrapPlan();
    const SparseBootstrapCertificate& certificate() const;
private:
    struct Data; std::shared_ptr<const Data> data_;
    friend SparseBootstrapPlan prepareSparseBootstrap(const SealAdapter&,const Cipher&,const SparseBootstrapInput&);
};
// Optional analysis preparation; bounds use shared exact finite-support helpers.
SparseBootstrapPlan prepareSparseBootstrap(const SealAdapter&,const Cipher&,const SparseBootstrapInput&);
BootstrapContractResult preflightSparseBootstrap(const SealAdapter&,const Cipher&,const SparseBootstrapPlan&);
RaisedCipher executeSparseBootstrap(SealAdapter&,const Cipher&,const SparseBootstrapPlan&);
}
