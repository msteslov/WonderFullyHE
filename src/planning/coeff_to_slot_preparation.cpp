#include "../core/coeff_to_slot_cert_internal.hpp"
#include "linear_transform_preparation_internal.hpp"
namespace m2424 {
using namespace linear_certificate;
CertifiedEvalRoundPlusCoeffToSlot EvalRoundPlusCoeffToSlot::prepareCertified(SealAdapter& a,const RaisedCipher& x,const BootstrapInputContext& source,const CoeffToSlotContract& hp,const CoeffToSlotContract& lp,const CoeffToSlotCertificationInput& input) const {
    return prepareCertifiedImpl(a,x,source,hp,lp,input,nullptr);
}
CertifiedEvalRoundPlusCoeffToSlot EvalRoundPlusCoeffToSlot::prepareCertified(SealAdapter& a,const SparseRaisedCipher& x,const SparseBootstrapPlan& sparse,const CoeffToSlotContract& hp,const CoeffToSlotContract& lp) const {
    const auto& c=sparse.certificate();
    const auto gate=sparseCoeffToSlotKeyGate(a,x.generation_,sparse);
    if(gate.status!=Status::Certified) {
        CertifiedEvalRoundPlusCoeffToSlot out; auto p=std::make_shared<CertifiedEvalRoundPlusCoeffToSlot::Impl>(); out.impl_=p;
        p->source=c.input; p->sparse=sparse;
        for(auto& t:p->traces) t.certificate=gate;
        p->domain.result=gate; p->domain.rho.provenance=gate.provenance; return out;
    }
    CoeffToSlotCertificationInput input{c.raisedMagnitude,c.messageMagnitude,c.sourceNoiseMagnitude};
    return prepareCertifiedImpl(a,x.cipher_,c.input,hp,lp,input,&sparse);
}
CertifiedEvalRoundPlusCoeffToSlot EvalRoundPlusCoeffToSlot::prepareCertifiedImpl(SealAdapter& a,const RaisedCipher& x,const BootstrapInputContext& source,const CoeffToSlotContract& hp,const CoeffToSlotContract& lp,const CoeffToSlotCertificationInput& input,const SparseBootstrapPlan* sparse) const {
    CertifiedEvalRoundPlusCoeffToSlot result; auto p=std::make_shared<CertifiedEvalRoundPlusCoeffToSlot::Impl>(); result.impl_=p;
    if(sparse) p->sparse=*sparse;
    p->source=source; p->hpContract=hp; p->lpContract=lp; p->factorization=plan_.factorization();
    p->domain.result={Status::RequiredBoundUnavailable,"EvalRoundDomain","Upstream message/noise coefficient bounds required"};
    p->domain.rho.provenance=p->domain.result.provenance;
    try {
        auto pre=preflight(a,x,source,hp,lp); if(pre.status!=Status::Certified) throw Failure{pre.status,pre.provenance};
        if(!known(input.raisedMagnitude)) throw Failure{Status::RequiredBoundUnavailable,"Raised canonical-embedding magnitude bound with provenance required"};
        if(a.info(x).ciphertextSize!=2) throw Failure{Status::InvalidInput,"Certified CtS requires size-two raised input"};
        const auto depth=plan_.depth(),N=plan_.polyModulusDegree();
        BootstrapBound incoming=bound(0,"Exact raised y=sigma(u/Delta0) is the semantic input; source noise belongs to u");
        std::optional<SlotToCoeffRuntimeStage> restoration;
        if(sparse) {
            const auto& sc=sparse->certificate();
            const auto backend=experimental::finiteSupportBackendKeyNoise();
            if(!known(backend)||!known(sc.restorationError)||sc.restorationError.upperBound<=0||
               !known(sc.restorationKeyNoise)||sc.restorationKeyNoise.upperBound<=0||
               !known(sc.restorationModDown)||sc.restorationModDown.upperBound<=0)
                throw Failure{Status::KeySwitchBoundUnavailable,"First-factor restoration requires known nonzero key-switch and ModDown bounds"};
            const mpq_class key=experimental::finiteSupportKeyNoise(N,q(backend.upperBound),source.raisedPrimes,source.specialPrime,q(a.info(x).scale));
            const mpq_class round=experimental::finiteSupportDivideRound(N,1,2,q(a.info(x).scale));
            if(q(sc.restorationError.upperBound)<key+round||q(sc.restorationKeyNoise.upperBound)<key||q(sc.restorationModDown.upperBound)<round)
                throw Failure{Status::KeySwitchBoundUnavailable,"Restoration envelope understates the shared audited backend bounds"};
            incoming=sc.restorationError;
            SlotToCoeffRuntimeStage stage; stage.operation="first-factor sparse restoration (shared HP/LP)";
            stage.activePrimes=source.raisedPrimes; stage.chainIndex=a.info(x).chainIndex;
            stage.inputScale=stage.outputScale=stage.arithmeticScale=scale(q(a.info(x).scale),a.info(x).scale);
            stage.magnitude=input.raisedMagnitude; stage.semanticError=incoming; stage.localError=incoming;
            stage.scaleRepresentationError=bound(0,"Identity switching leaves the exact binary64 scale unchanged");
            stage.centeredHeadroomNumerator=headroom(stage.activePrimes,q(a.info(x).scale),q(input.raisedMagnitude.upperBound)+q(incoming.upperBound),"CtS first-factor restoration");
            stage.headroomProvenance="Actual raised modulus; restoration precedes the first factor and leaves ideal u unchanged";
            restoration=stage;
        }
        for(std::size_t b=0;b<2;++b) {
            auto& trace=p->traces[b]; trace.gate=b?BootstrapGate::CoeffToSlotLP:BootstrapGate::CoeffToSlotHP;
            trace.prefactor=b?CoeffToSlotPrefactor::sourceNormalization(source):CoeffToSlotPrefactor{};
            trace.inputMagnitude=input.raisedMagnitude;
            trace.inputSemanticError=incoming; trace.firstFactorRestoration=restoration;
            auto layout=certificationLayout(trace.prefactor);
            // Shared preparation uses 2^bitlength(actual prime), exactly the legacy CtS schedule.
            SlotToCoeffContract c; c.inputMagnitude[0]=c.inputMagnitude[1]=input.raisedMagnitude;
            c.inputError[0]=c.inputError[1]=incoming;
            auto prepared=prepareRootLinearTransform(layout,a,x.cipher_,x.cipher_,c); p->branches[b]=prepared;
            if(prepared->certificate.result.status!=Status::Certified) throw Failure{prepared->certificate.result.status,prepared->certificate.result.provenance};
            double finalError=0;
            for(std::size_t h=0;h<2;++h) {
                for(std::size_t r=0;r<depth;++r) {
                    auto factor=prepared->certificate.factors[h*depth+r];
                    factor.incomingSemanticError=r?prepared->certificate.factors[h*depth+r-1].propagation.semanticError:incoming;
                    if(sparse&&r==0) {
                        factor.arithmeticTerms.insert(factor.arithmeticTerms.begin(),{
                            {"incoming restoration key noise (not B_local)",sparse->certificate().restorationKeyNoise},
                            {"incoming restoration ModDown (not B_local)",sparse->certificate().restorationModDown}});
                    }
                    CoeffToSlotFactorTrace old; old.operation="certified BSGS factor/rescale"; old.activePrimes=factor.runtime.back().activePrimes;
                    old.outputState=a.info(x);
                    auto stateBits=factor.runtime.back().outputScale.binary64Bits;
                    std::memcpy(&old.outputState.scale,&stateBits,8);
                    old.outputState.chainIndex=factor.runtime.back().chainIndex;
                    old.outputState.coeffModulusSize=old.activePrimes.size(); old.outputState.coeffModulusLog2=0;
                    for(auto prime:old.activePrimes) old.outputState.coeffModulusLog2+=std::log2(static_cast<double>(prime));
                    old.kappa=factor.bounds.kappa; old.diagonalError=factor.bounds.delta; old.localAddedError=factor.bounds.localArithmeticError;
                    old.idealMagnitude=factor.propagation.idealMagnitude; old.propagatedSemanticError=factor.propagation.semanticError;
                    trace.halves[h].push_back(std::move(old)); trace.certifiedFactors.push_back(std::move(factor));
                }
                const auto& last=trace.certifiedFactors.back();
                const auto state=last.runtime.back(); const auto M=last.propagation.idealMagnitude,E=last.propagation.semanticError;
                double outputScale; auto sb=state.outputScale.binary64Bits; std::memcpy(&outputScale,&sb,8);
                const mpq_class local=experimental::finiteSupportKeyNoise(N,q(c.evaluationKeyNoiseSupport.upperBound),state.activePrimes,source.specialPrime,q(outputScale))
                    +experimental::finiteSupportDivideRound(N,1,2,q(outputScale));
                SlotToCoeffFactorTrace projection; projection.branch=h; projection.factor=depth; projection.incomingSemanticError=E;
                projection.bounds={bound(2,"Exact real projection v+conj(v) has real-linear infinity norm two"),bound(0,"Projection has no encoded coefficients"),bound(local,"Shared finite-support conjugation key noise plus two-component ModDown; "+c.evaluationKeyNoiseSupport.provenance)};
                projection.arithmeticTerms={
                    {"conjugation key noise",bound(experimental::finiteSupportKeyNoise(N,q(c.evaluationKeyNoiseSupport.upperBound),state.activePrimes,source.specialPrime,q(outputScale)),"Shared key noise at actual active primes and scale")},
                    {"conjugation ModDown",bound(experimental::finiteSupportDivideRound(N,1,2,q(outputScale)),"Shared two-component divide-round")}
                };
                projection.propagation=propagateLinearTransform(M,E,projection.bounds);
                if(projection.propagation.result.status!=Status::Certified) throw Failure{Status::RequiredBoundUnavailable,"Projection recurrence unavailable"};
                auto stage=state; stage.factor=depth; stage.inputScale=stage.arithmeticScale=stage.outputScale;
                stage.operation="conjugation"; stage.magnitude=M; stage.semanticError=bound(q(E.upperBound)+local,"Incoming error plus actual conjugation arithmetic");
                stage.localError=projection.bounds.localArithmeticError; stage.scaleRepresentationError=bound(0,"Conjugation preserves scale exactly");
                stage.centeredHeadroomNumerator=headroom(stage.activePrimes,q(outputScale),q(M.upperBound)+q(stage.semanticError.upperBound),"conjugation");
                projection.runtime.push_back(stage);
                stage.operation="real projection addition"; stage.magnitude=projection.propagation.idealMagnitude; stage.semanticError=projection.propagation.semanticError;
                stage.localError=bound(0,"Same-scale RNS addition is exact");
                stage.centeredHeadroomNumerator=headroom(stage.activePrimes,q(outputScale),q(stage.magnitude.upperBound)+q(stage.semanticError.upperBound),"projection addition");
                projection.runtime.push_back(stage);
                auto compatibility=trace.halves[h].back(); compatibility.operation="certified conjugation/add projection";
                compatibility.kappa=projection.bounds.kappa; compatibility.diagonalError=projection.bounds.delta;
                compatibility.localAddedError=projection.bounds.localArithmeticError;
                compatibility.idealMagnitude=projection.propagation.idealMagnitude; compatibility.propagatedSemanticError=projection.propagation.semanticError;
                trace.halves[h].push_back(std::move(compatibility)); trace.certifiedFactors.push_back(std::move(projection));
                finalError=std::max(finalError,stage.semanticError.upperBound);
            }
            trace.outputError=bound(finalError,"Maximum of both half certificates, including phased first factor and conjugation/add projection");
            trace.levelsConsumed=depth; trace.rescaleOperations=2*depth;
            const auto& requested=b?lp:hp;
            double finalScale; auto finalBits=prepared->certificate.outputScale.binary64Bits; std::memcpy(&finalScale,&finalBits,8);
            if(std::abs(std::log2(finalScale)-requested.outputScaleLog2)>requested.inputScaleToleranceLog2)
                throw Failure{Status::ScaleScheduleInfeasible,"Prepared output scale does not satisfy branch contract"};
            trace.errorBudget={finalError<=requested.maxAbsError?Status::Certified:Status::ErrorBudgetExceeded,
                b?"CoeffToSlotLP.errorBudget":"CoeffToSlotHP.errorBudget","Computed deterministic error compared with requested branch budget; no tolerance adjustment"};
            trace.evidence={true,"Exact root/rational operator preparation and shared finite-support arithmetic/headroom verified",{trace.outputError}};
            trace.certificate={Status::Certified,b?"CoeffToSlotLP":"CoeffToSlotHP",trace.evidence.provenance};
        }
        p->domain=certifyCoeffToSlotDomain(source,input,p->traces[1].outputError);
    } catch(const Failure& f) {
        for(auto& t:p->traces) { t.certificate={f.status,"CoeffToSlot",f.why}; t.evidence.verified=false; }
    } catch(const std::exception& e) {
        for(auto& t:p->traces) { t.certificate={Status::RequiredBoundUnavailable,"CoeffToSlot",e.what()}; t.evidence.verified=false; }
    }
    return result;
}
CoeffToSlotDomainCertificate certifyCoeffToSlotDomain(const BootstrapInputContext& source,const CoeffToSlotCertificationInput& input,const BootstrapBound& lpError) {
    CoeffToSlotDomainCertificate result;
    try {
        if(!known(input.messageMagnitude)||!known(input.sourceNoiseMagnitude)||!known(lpError))
            throw Failure{Status::RequiredBoundUnavailable,"LP domain requires deterministic message/noise coefficient bounds and certified LP error, each with provenance"};
        if(source.sourcePrimes.empty()) throw Failure{Status::MissingExactModulusContext,"LP domain requires exact source primes"};
        mpz_class Q=1; for(auto prime:source.sourcePrimes) {
            if(prime<2) throw Failure{Status::MissingExactModulusContext,"Invalid source modulus factor"};
            Q*=z(prime);
        }
        result.rho=bound((q(input.messageMagnitude.upperBound)+q(input.sourceNoiseMagnitude.upperBound))/Q+q(lpError.upperBound),
            "(|m|+|nu_b|)/exact product(sourcePrimes) + certified E_LP; "+input.messageMagnitude.provenance+"; "+input.sourceNoiseMagnitude.provenance+"; "+lpError.provenance);
        result.result={result.rho.upperBound<.5?Status::Certified:Status::DomainViolation,"EvalRoundDomain","Computed deterministic rho; EvalRound requires rho < 1/2"};
    } catch(const Failure& f) { result.result={f.status,"EvalRoundDomain",f.why}; result.rho.provenance=f.why; }
    return result;
}

}
