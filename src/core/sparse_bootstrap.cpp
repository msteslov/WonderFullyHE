#include "sparse_bootstrap_internal.hpp"
#include <cmath>
#include <algorithm>
#include <stdexcept>
namespace m2424 {
BootstrapSecurityReport certifyPublicRlweFamilies(std::vector<PublicRlweFamily> families,
    const std::vector<RlweSecurityEvidence>& evidence,int target) {
    using S=BootstrapCertificationStatus;
    BootstrapSecurityReport out; out.families=std::move(families);
    out.result={S::Certified,"security","All public families have matching external estimator evidence"};
    double minimum=INFINITY; bool allKnown=!out.families.empty();
    if(out.families.empty()||target<=0) out.result={S::SecurityBudgetExceeded,"security","Complete nonempty public-family inventory and positive target required"};
    for(auto& f:out.families) {
        f.evidence.reset(); f.result={S::SecurityBudgetExceeded,"security/"+f.id,"Missing concrete estimator evidence for this family, modulus and known relations"};
        std::size_t matches=0;
        for(const auto& e:evidence) if(e.familyId==f.id) { ++matches; f.evidence=e; }
        bool accepted=false;
        if(matches==1&&f.evidence) {
            const auto& e=*f.evidence;
            const bool valid=e.statement==f.statement&&!f.statement.empty()&&e.lowerSecurityBits&&
                std::isfinite(*e.lowerSecurityBits)&&*e.lowerSecurityBits>=0&&!e.estimator.empty()&&
                !e.version.empty()&&!e.artifact.empty()&&!e.assumptions.empty();
            if(valid&&(!f.searchSpaceCeilingBits||*e.lowerSecurityBits<=*f.searchSpaceCeilingBits)) {
                accepted=true; minimum=std::min(minimum,*e.lowerSecurityBits);
                if(*e.lowerSecurityBits>=target) f.result={S::Certified,"security/"+f.id,"External evidence: "+e.artifact+"; assumptions: "+e.assumptions};
                else f.result.provenance="Concrete family security below requested target";
            } else f.result.provenance="Invalid/stale estimator statement or estimate exceeds signed sparse search-space ceiling";
        }
        allKnown=allKnown&&accepted;
        if(f.result.status!=S::Certified&&out.result.status==S::Certified) out.result=f.result;
    }
    // A failed/Unknown family cannot be hidden by the minimum of known families.
    if(allKnown) out.minimumSecurityBits=minimum;
    return out;
}
SparseBootstrapPlan::SparseBootstrapPlan():data_(std::make_shared<Data>()) {}
const SparseBootstrapCertificate& SparseBootstrapPlan::certificate() const { return data_->certificate; }
BootstrapContractResult preflightSparseBootstrap(const SealAdapter& a,const Cipher& c,const SparseBootstrapPlan& p) {
    using S=BootstrapCertificationStatus; const auto& cert=p.certificate();
    if(!a.hasSparseEncapsulationKey()||!a.hasSparseRestorationKey()) return {S::MissingEvaluationKeys,"sparse.keys","Both directional switch keys required"};
    const auto k=a.sparseKeyMetadata();
    if(k.generation!=cert.keys.generation||k.weight!=cert.keys.weight||k.context!=cert.keys.context||k.degree!=cert.keys.degree)
        return {S::InvalidInput,"sparse.keys","Prepared secret generation/weight/context mismatch"};
    if(k.publicKeySamples!=cert.keys.publicKeySamples||k.relinSamples!=cert.keys.relinSamples||
       k.galoisSamples!=cert.keys.galoisSamples||k.galoisElements!=cert.keys.galoisElements||
       k.encryptionModuli!=cert.keys.encryptionModuli)
        return {S::InvalidInput,"sparse.inventory","Published material changed after security statement preparation"};
    const auto source=resolveBootstrapInput(a,c,cert.input.scaleBinary64Bits);
    if(!source.context) return source.result;
    if(source.context->sourcePrimes!=cert.input.sourcePrimes||source.context->raisedPrimes!=cert.input.raisedPrimes||source.context->contextFingerprint!=cert.input.contextFingerprint||a.info(c).ciphertextSize!=2)
        return {S::InvalidInput,"sparse.context","Prepared active modulus/context/components mismatch"};
    return cert.result;
}
RaisedCipher executeSparseBootstrap(SealAdapter& a,const Cipher& c,const SparseBootstrapPlan& p) {
    const auto gate=preflightSparseBootstrap(a,c,p);
    if(gate.status!=BootstrapCertificationStatus::Certified) throw std::invalid_argument(gate.provenance);
    // Arithmetic readiness is distinct from security. No global certificate is
    // issued by this stage, just as the PR-3/4 executors expose local readiness.
    auto sparse=a.encapsulateSparse(c); auto raised=a.modRaiseSparse(sparse); return a.restoreSparse(raised);
}
}
