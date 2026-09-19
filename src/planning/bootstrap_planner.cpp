#include "../core/bootstrap_internal.hpp"
#include "../core/coeff_to_slot_cert_internal.hpp"
#include "certified_arithmetic_internal.hpp"
#include "m2424/experimental/evalmod_analysis/evalround_execution.hpp"
#include <cfenv>
#include <sstream>
namespace m2424 {
namespace {
using namespace experimental::arithmetic;
using S=BootstrapCertificationStatus;
struct Rejection { BootstrapContractResult gate; };
void require(bool okay,S status,const std::string& gate,const std::string& why) { if(!okay) throw Rejection{{status,gate,why}}; }
double runtime(std::uint64_t bits) { double x; std::memcpy(&x,&bits,8); return x; }
double down(const mpq_class& x) { double d=x.get_d(); if(q(d)>x) d=std::nextafter(d,0.); return d; }
Cipher stateToken(SealAdapter& a,const Cipher& source,std::size_t chain,double scale) {
    (void)source;
    // Tokens are public encryptions used only for parameter/constant preparation.
    // They never enter the output arithmetic or provide semantic evidence.
    auto top=a.encrypt(a.encodeScalarAtScaleFor(0,scale,a.encrypt(a.encode({0.}))));
    return a.modSwitchToChainIndex(top,chain);
}
BootstrapBound intrinsicStC(const SlotToCoeffCertificate& c,const BootstrapBound& magnitude) {
    BootstrapBound total[2]; const auto depth=c.metrics.depth;
    for(std::size_t h=0;h<2;++h) {
        auto M=magnitude,E=bound(0,"Zero incoming error used only to isolate the affine StC intrinsic term");
        for(std::size_t r=0;r<depth;++r) { auto out=propagateLinearTransform(M,E,c.factors[h*depth+r].bounds); require(out.result.status==S::Certified,S::RequiredBoundUnavailable,"StC","Intrinsic recurrence unavailable"); M=out.idealMagnitude; E=out.semanticError; }
        total[h]=E;
    }
    return bound(q(total[0].upperBound)+q(total[1].upperBound),"StC intrinsic error from the same certified factors at zero incoming error; local bounds reserved for full input envelope");
}
// Uniform reservation, using the same finite-support component rounding and
// linear recurrence as actual arithmetic. All actual scales are checked against
// these minima and the actual compiled local error is checked against this cap.
BootstrapBound combinationReserve(std::size_t N,const BootstrapRequest& request,const mpq_class& gamma,const BootstrapBound& hpError,const mpq_class& hpMagnitude) {
    const mpq_class unit=mpq_class(1)/integer(std::uint64_t(1)<<50);
    const mpq_class Ecap=q(request.target.targetAbsoluteError)/(gamma<1?gamma:mpq_class(1)),K=*request.lift.K;
    const mpq_class R=experimental::finiteSupportDivideRound(N,1,2,q(request.minimumArithmeticScale));
    mpq_class conditioning=0;
    for(std::size_t i=0;i<request.maxIntegerRescales;++i) conditioning+=R+unit*(K+Ecap+conditioning+R);
    const mpq_class delta=mpq_class(1)/(2*q(request.minimumConstantScale));
    const mpq_class finalRound=experimental::finiteSupportDivideRound(N,1,2,q(request.combinationScale));
    auto scalarLocal=[&](const mpq_class& M,const mpq_class& E,const mpq_class& k) {
        mpq_class encoding=delta*(M+E);
        mpq_class representation=unit*(k*(M+E)+encoding);
        return mpq_class(encoding+representation+finalRound+unit*(k*(M+E)+encoding+representation+finalRound));
    };
    return bound(gamma*conditioning+scalarLocal(K,Ecap+conditioning,gamma)+scalarLocal(hpMagnitude,q(hpError.upperBound),1),
        "Uniform combination reservation: finite-support rescale conditioning, exact scalar rounding <=1/(2 minimumConstantScale), normal binary64 ratio <=2^-50, both final rescales; actual DAG must fit reservation");
}
double alignmentScale(double input,double target,std::uint64_t prime,double minimum) {
    double p=target*static_cast<double>(prime)/input;
    for(int direction:{0,-1,1}) {
        double value=p;
        for(int i=0;i<32;++i) {
            if(value>=minimum&&std::isfinite(value)&&value>0&&((input*value)/static_cast<double>(prime))==target) return value;
            value=std::nextafter(value,direction<0?0.:INFINITY);
        }
    }
    throw Rejection{{S::ScaleScheduleInfeasible,"combination","No real multiply/rescale aligns to requested output scale"}};
}
}
BootstrapPlan Bootstrapper::prepare(SealAdapter& a,const Cipher& input,const BootstrapRequest& request) const {
    using namespace experimental::arithmetic;
    BootstrapPlan plan; auto p=std::make_shared<BootstrapPlan::Data>(); plan.data_=p;
    p->request=request; p->degree=cts_.plan().polyModulusDegree(); p->ctsDepth=cts_.plan().depth(); p->stcDepth=stc_.metrics().depth; p->keys=rotationKeys();
    auto& trace=p->trace; trace.plan.id="EvalRoundPlus/deterministic-baseline";
    trace.details["pipeline"]="input -> ModRaise -> certified CtS(HP,LP) -> EvalRoundExecutionPlan(LP) -> certified HP-gamma*I_hat -> PreparedSlotToCoeffPlan -> output";
    for(const auto* name:{"E_HP","E_LP","rho_cert","K","gamma","Gamma_StC","requiredIntegerError","E_I","B_comb","E_pair","E_StC","E_boot"})
        trace.bounds[name].provenance="Not computed; prerequisite gates pending";
    trace.details["headroom"]="Unknown until all executable stages are prepared";
    trace.details["executed"]="false"; trace.details["failureProbabilityLog2"]="-inf (deterministic, no primary probabilistic events)";
    trace.details["securityBits"]=std::to_string(a.securityLevelBits());
    trace.details["K.provenance"]=request.lift.provenance;
    trace.details["K.evidence"]=request.lift.evidence==BootstrapLiftEvidence::TestFixtureAssumption?"TEST-ONLY assumption; never a production certificate":request.lift.evidence==BootstrapLiftEvidence::Analytical?"external analytical evidence":"Unknown";
    auto note=[&](BootstrapContractResult g) { trace.gateResults.push_back(g); if(g.status!=S::Certified&&trace.result.status==S::Certified) trace.result=g; };
    trace.result={S::Certified,"planning","No gate failure encountered yet"};
    try {
        const auto validTarget=validateBootstrapTarget(request.target); require(validTarget.status==S::Certified,validTarget.status,validTarget.gate,validTarget.provenance);
        require(std::fegetround()==FE_TONEAREST,S::ScaleScheduleInfeasible,"schedule","Round-to-nearest required");
        for(double scale:{request.combinationScale,request.minimumArithmeticScale,request.minimumConstantScale}) require(std::isfinite(scale)&&scale>=1&&scale<=std::ldexp(1.,900),S::ScaleScheduleInfeasible,"schedule","Finite normal positive scale limits required");
        require(request.maxIntegerRescales<=8&&request.maxCleaningRounds<=16,S::InvalidInput,"schedule","Bounded baseline search required");
        auto resolved=resolveBootstrapInput(a,input); require(bool(resolved.context),resolved.result.status,resolved.result.gate,resolved.result.provenance);
        trace.plan.input=*resolved.context; const auto& source=*trace.plan.input;
        require(a.slotCount()*2==p->degree&&a.info(input).ciphertextSize==2,S::InvalidInput,"input","Degree and two-component input required");
        if(a.sparseKeyMetadata().generation) {
            SparseBootstrapInput upstream;
            upstream.messageMagnitude=request.upstream.messageMagnitude;
            upstream.sourceNoiseMagnitude=request.upstream.sourceNoiseMagnitude;
            upstream.evaluationKeyNoiseSupport=request.sparseEvaluationKeyNoiseSupport;
            upstream.securityEvidence=request.securityEvidence; upstream.targetSecurityBits=request.target.targetSecurityBits;
            const auto sparse=prepareSparseBootstrap(a,input,upstream);
            const auto& cert=sparse.certificate();
            trace.publicRlwe=cert.security;
            trace.certificate.gates[static_cast<std::size_t>(BootstrapGate::Security)]={
                cert.security.result.status==S::Certified,cert.security.result.provenance,{}};
            if(cert.security.minimumSecurityBits) trace.certificate.minimumSecurityBits=static_cast<int>(std::floor(*cert.security.minimumSecurityBits));
            trace.details["sparse.schedule"]="encapsulation at source chain/scale -> centered ModRaise under s_b -> CtS first-factor shared restoration at raised chain/scale -> factorized HP/LP under s; no standalone restoration";
            trace.details["h_b"]=std::to_string(cert.keys.weight);
            trace.details["sparse.distribution"]=cert.keys.distribution;
            trace.details["sparse.keyGeneration"]=std::to_string(cert.keys.generation);
            trace.details["K.provenance"]=cert.liftProvenance;
            trace.details["K.evidence"]="Derived from generated sparse key certificate; caller lift ignored";
            trace.bounds["encapsulation.KS"]=cert.encapsulationError;
            trace.bounds["restoration.KS"]=cert.restorationError;
            trace.bounds["nu_b"]=cert.sourceNoiseMagnitude;
            trace.bounds["raisedMagnitude"]=cert.raisedMagnitude;
            trace.bounds["restoredMagnitude"]=cert.restoredMagnitude;
            trace.bounds["rho.beforeCtS"]=cert.rhoBeforeCoeffToSlot;
            trace.details["lambda_boot"]=cert.security.minimumSecurityBits?std::to_string(*cert.security.minimumSecurityBits):"Unknown: missing audited per-family concrete estimator evidence";
            require(cert.result.status==S::Certified,cert.result.status,cert.result.gate,cert.result.provenance);
            trace.bounds["K"]=bound(cert.K,cert.liftProvenance);
            trace.certificate.gates[static_cast<std::size_t>(BootstrapGate::SparseSecret)]={true,cert.liftProvenance,{trace.bounds["K"]}};
            trace.certificate.gates[static_cast<std::size_t>(BootstrapGate::KeySwitch)]={true,"Two separate audited identity-key-switch operations",{cert.encapsulationError,cert.restorationError}};
            trace.gateResults.push_back(cert.result);
            // Prepare the owned first-factor path while preserving the K>1
            // execution rejection. Parameter/key shortages remain explicit CtS
            // diagnostics; they cannot authorize a narrower EvalRound domain.
            auto sparseRaised=a.modRaiseSparse(a.encapsulateSparse(input));
            CoeffToSlotContract schedule{"bootstrap/sparse-CtS",a.slotCount(),p->degree,
                std::log2(a.scale(input)),std::log2(a.scale(input)),.25,request.target.targetAbsoluteError};
            p->cts=cts_.prepareCertified(a,sparseRaised,sparse,schedule,schedule);
            trace.bounds["E_HP"]=p->cts->hp().outputError;
            trace.bounds["E_LP"]=p->cts->lp().outputError;
            trace.bounds["rho_cert"]=p->cts->domain().rho;
            trace.gateResults.push_back(p->cts->hp().certificate);
            trace.gateResults.push_back(p->cts->lp().certificate);
            trace.gateResults.push_back(p->cts->domain().result);
            trace.certificate.gates[static_cast<std::size_t>(BootstrapGate::CoeffToSlotHP)]=p->cts->hp().evidence;
            trace.certificate.gates[static_cast<std::size_t>(BootstrapGate::CoeffToSlotLP)]=p->cts->lp().evidence;
            require(cert.K==1,S::UnsupportedEvalRoundDomain,"EvalRound.domain","Derived sparse K="+std::to_string(cert.K)+" exceeds the certified K=1 ciphertext executor; h is never reduced");
            throw Rejection{{S::UnsupportedEvalRoundDomain,"EvalRound.domain","Sparse K=1 path is not enabled"}};
        }
#ifndef M2424_ENABLE_BOOTSTRAP_FIXTURE
        require(request.lift.evidence!=BootstrapLiftEvidence::TestFixtureAssumption,S::TestOnlyAssumption,"lift.fixture","Test fixture execution is disabled in production builds");
#endif
        require(request.lift.evidence!=BootstrapLiftEvidence::Analytical,S::LiftBoundUnavailable,"lift","Caller production K is not trusted; generate sparse bootstrap keys");
        require(a.hasRotationKeys(p->keys)&&a.hasRelinKeys()&&a.hasConjugationKey(),S::MissingEvaluationKeys,"keys","CtS/StC rotations, conjugation and relinearization keys required before execution");
        trace.details["keys"]="relinearization, conjugation, rotations:"; for(auto key:p->keys) trace.details["keys"]+=std::to_string(key)+",";
        require(request.lift.K.has_value()&&request.lift.evidence!=BootstrapLiftEvidence::Unknown&&!request.lift.provenance.empty(),S::LiftBoundUnavailable,"lift","Analytical |I|<=K unavailable; ordinary ModRaise does not imply K=1");
        trace.bounds["K"]=bound(*request.lift.K,request.lift.provenance);
        require(*request.lift.K==1,S::UnsupportedEvalRoundDomain,"lift","The ciphertext executor supports only K=1; supplied K is not substituted");
        note({S::Certified,"lift",request.lift.provenance});
        require(known(request.upstream.messageMagnitude)&&known(request.upstream.sourceNoiseMagnitude)&&known(request.upstream.raisedMagnitude),S::RequiredBoundUnavailable,"upstream","Known message/noise coefficient and raised magnitude bounds required");
        mpz_class Q=1; for(auto prime:source.sourcePrimes) Q*=integer(prime);
        const mpq_class Delta(runtime(source.scaleBinary64Bits)),gamma=mpq_class(Q)/Delta;
        trace.details["gamma.exact"]=gamma.get_str(); trace.bounds["gamma"]=bound(gamma,"Exact product(sourcePrimes) / binary64 dyadic Delta0");
        const mpq_class useful=(q(request.upstream.messageMagnitude.upperBound)+q(request.upstream.sourceNoiseMagnitude.upperBound))/Delta;
        const mpq_class hpMagnitude=useful+gamma*(*request.lift.K);
        trace.bounds["sourceNoise.output"]=bound(p->degree*q(request.upstream.sourceNoiseMagnitude.upperBound)/Delta,"Canonical embedding coefficient-noise sup bound N*M_nu/Delta0");
        const std::size_t top=source.raisedPrimes.size()-1;
        require(request.outputChainIndex<=top&&p->stcDepth<=top-request.outputChainIndex,S::InsufficientLevels,"levels","Requested output/StC schedule exceeds actual chain");
        const std::size_t combinationChain=request.outputChainIndex+p->stcDepth;
        require(combinationChain<top&&p->ctsDepth+2<top-combinationChain,S::InsufficientLevels,"levels","Need CtS, extraction, arithmetic combination, StC and requested output modulus before execution");
        trace.details["ModRaise.chainIndex"]=std::to_string(top); trace.details["ModRaise.scaleBits"]=std::to_string(source.scaleBinary64Bits);
        auto raised=a.modRaiseToTop(input);
        CoeffToSlotContract c{"bootstrap/CtS",a.slotCount(),p->degree,std::log2(a.scale(input)),std::log2(a.scale(input)),.25,request.target.targetAbsoluteError};
        p->cts=cts_.prepareCertified(a,raised,source,c,c,request.upstream);
        require(p->cts->hp().certificate.status==S::Certified,p->cts->hp().certificate.status,"CoeffToSlotHP",p->cts->hp().certificate.provenance);
        require(p->cts->lp().certificate.status==S::Certified,p->cts->lp().certificate.status,"CoeffToSlotLP",p->cts->lp().certificate.provenance);
        trace.bounds["E_HP"]=p->cts->hp().outputError; trace.bounds["E_LP"]=p->cts->lp().outputError;
        note(p->cts->hp().errorBudget); note(p->cts->lp().errorBudget);
        trace.bounds["rho_cert"]=p->cts->domain().rho;
        require(p->cts->domain().result.status==S::Certified,p->cts->domain().result.status,"domain",p->cts->domain().result.provenance);
        note(p->cts->domain().result);
        auto pairState=stateToken(a,input,combinationChain,request.combinationScale);
        SlotToCoeffContract sc; sc.inputMagnitude[0]=sc.inputMagnitude[1]=bound(useful,"Useful coefficient half (m+nu_b)/Delta0 after ideal lift cancellation");
        sc.inputError[0]=sc.inputError[1]=bound(q(request.target.targetAbsoluteError),"Reserved pair error envelope; checked against actual combination certificate");
        auto reserved=stc_.prepare(a,pairState,pairState,sc);
        require(reserved.certificate().result.status==S::Certified,reserved.certificate().result.status,"StC",reserved.certificate().result.provenance);
        trace.bounds["Gamma_StC"]=reserved.certificate().gamma;
        require(known(trace.bounds["Gamma_StC"])&&trace.bounds["Gamma_StC"].upperBound>=1,S::SlotToCoeffGainUnavailable,"StC.gain","Actual prepared gain >=1 required for reservation");
        auto intrinsic=intrinsicStC(reserved.certificate(),sc.inputMagnitude[0]); trace.bounds["StC.intrinsic.reserved"]=intrinsic;
        auto combReserve=combinationReserve(p->degree,request,gamma,trace.bounds["E_HP"],hpMagnitude); trace.bounds["B_comb.reserved"]=combReserve;
        const mpq_class gain=q(trace.bounds["Gamma_StC"].upperBound);
        const mpq_class fixed=gain*(q(trace.bounds["E_HP"].upperBound)+q(combReserve.upperBound))+q(intrinsic.upperBound)+q(trace.bounds["sourceNoise.output"].upperBound);
        trace.bounds["fixedOutputError"]=bound(fixed,"Gamma_StC*(E_HP+B_comb reserve)+StC intrinsic+source noise; no double amplification");
        const mpq_class remainder=q(request.target.targetAbsoluteError)-fixed;
        trace.details["backward.remainder.exact"]=remainder.get_str();
        require(remainder>0,S::ErrorBudgetExceeded,"backwardBudget","Known deterministic output-error terms exhaust global budget; EvalRound not compiled or executed");
        const double required=down(remainder/(gain*gamma));
        require(required>0&&std::isfinite(required),S::ErrorBudgetExceeded,"backwardBudget","Positive representable integer budget required");
        trace.bounds["requiredIntegerError"]={required,BootstrapBoundKind::Deterministic,"Downward rounded residual/(actual Gamma_StC * exact gamma)",{}};
        auto ctsScale=runtime(p->cts->lp().certifiedFactors.back().runtime.back().outputScale.binary64Bits);
        auto ctsState=stateToken(a,input,top-p->ctsDepth,ctsScale);
        EvalRoundProblem problem{1,trace.bounds["rho_cert"].upperBound,required,request.maxCleaningRounds};
        auto candidate=makeEvalRoundReferenceCandidate(problem,EvalRoundRadix::Binary,EvalRoundExtractionMethod::BinaryQuadraticK1);
        auto math=planEvalRoundCandidate(problem,candidate);
        experimental::EvalRoundExecutionOptions options; options.inputSemanticError=trace.bounds["E_LP"];
        p->evalRound=experimental::EvalRoundExecutionCompiler::compile(a,ctsState,math,options);
        require(p->evalRound.certification().status==S::Certified,p->evalRound.certification().status,"EvalRound",p->evalRound.certification().provenance);
        trace.bounds["E_I"]=bound(q(p->evalRound.integerErrorUpper()),p->evalRound.certification().provenance);
        trace.details["EvalRound.extractor"]=p->evalRound.mathematicalPlan().candidateId; trace.details["EvalRound.radix"]="2";
        for(const auto& digit:p->evalRound.mathematicalPlan().digits) trace.details["EvalRound.cleaning"]+=std::to_string(digit.cleaningIterations)+",";
        const auto& end=p->evalRound.nodes()[p->evalRound.outputNode()];
        Builder b(a,ctsState,experimental::finiteSupportBackendKeyNoise(),"Bootstrap combination");
        auto hp=b.input(ctsScale,hpMagnitude,q(trace.bounds["E_HP"].upperBound));
        auto ii=b.input(runtime(end.outputScale.binary64Bits),1,p->evalRound.integerErrorUpper(),a.chainIndex(ctsState)-end.chainIndex);
        std::size_t rescaleCount=0;
        const auto prime=source.raisedPrimes.at(combinationChain+1);
        while((request.combinationScale*static_cast<double>(prime))/b.states[ii].scale<request.minimumConstantScale) {
            require(++rescaleCount<=request.maxIntegerRescales,S::InsufficientLevels,"combination","Integer normalization exceeds reserved rescale count");
            ii=b.add(Op::Rescale,{ii},"integer scale reduction");
            require(b.states[ii].scale>=request.minimumArithmeticScale,S::ScaleScheduleInfeasible,"combination","Actual arithmetic scale below certified reservation");
        }
        const auto level=a.chainIndex(ctsState)-(combinationChain+1);
        require(b.states[ii].level<=level,S::InsufficientLevels,"combination","EvalRound and scale reduction exceed reserved modulus schedule");
        hp=b.switchLevel(hp,level,"HP modulus alignment"); ii=b.switchLevel(ii,level,"integer modulus alignment");
        hp=b.scalar(hp,1,alignmentScale(b.states[hp].scale,request.combinationScale,prime,request.minimumConstantScale),"HP arithmetic scale alignment");
        ii=b.scalar(ii,gamma,alignmentScale(b.states[ii].scale,request.combinationScale,prime,request.minimumConstantScale),"exact gamma multiplication");
        hp=b.add(Op::Rescale,{hp},"HP alignment rescale"); ii=b.add(Op::Rescale,{ii},"gamma rescale");
        auto output=b.add(Op::Subtract,{hp,ii},"HP - gamma*I_hat");
        // Replay the very same arithmetic with zero incoming errors to isolate
        // B_comb; reserve is checked against full-error-minus-propagated bound.
        const mpq_class pairE=b.states[output].E;
        const mpq_class combLocal=pairE-q(trace.bounds["E_HP"].upperBound)-gamma*q(p->evalRound.integerErrorUpper());
        require(combLocal>=0&&combLocal<=q(combReserve.upperBound),S::RequiredBoundUnavailable,"combination","Actual local error exceeds backward reservation");
        trace.bounds["B_comb"]=bound(combLocal,"Shared certified arithmetic DAG total minus exact propagated HP and gamma*E_I terms");
        trace.bounds["E_pair"]=bound(pairE,"E_HP + exact gamma*E_I + B_comb from actual arithmetic schedule");
        require(pairE<=q(request.target.targetAbsoluteError),S::ErrorBudgetExceeded,"combination","Pair error exceeds reserved StC input envelope");
        p->combination=b.nodes; p->combinationOutput=output; p->constants.resize(b.nodes.size());
        for(std::size_t i=0;i<b.nodes.size();++i) if(b.nodes[i].operation==Op::MultiplyPlain||b.nodes[i].operation==Op::AddPlain) {
            std::vector<std::uint64_t> residues; for(auto mod:b.nodes[i].activePrimes) { mpz_class r; mpz_mod(r.get_mpz_t(),b.rounded[i].get_mpz_t(),integer(mod).get_mpz_t()); residues.push_back(std::stoull(r.get_str())); }
            p->constants[i]=a.encodeScalarRnsAtScaleFor(residues,b.nodes[i].constantScale,ctsState,b.states[i].level);
        }
        sc.inputError[0]=sc.inputError[1]=trace.bounds["E_pair"];
        p->stc=stc_.prepare(a,pairState,pairState,sc);
        require(p->stc.certificate().result.status==S::Certified,p->stc.certificate().result.status,"StC",p->stc.certificate().result.provenance);
        require(p->stc.certificate().gamma.upperBound==trace.bounds["Gamma_StC"].upperBound,S::SlotToCoeffGainUnavailable,"StC","Actual gain differs from backward prepared gain");
        trace.bounds["E_StC"]=p->stc.certificate().outputError;
        trace.bounds["E_boot"]=bound(q(trace.bounds["E_StC"].upperBound)+q(trace.bounds["sourceNoise.output"].upperBound),"Actual StC output certificate already includes amplification of E_pair; add only original source noise");
        trace.plan.levelsUsed=top-request.outputChainIndex;
        trace.details["outputScale.exact"]=p->stc.certificate().outputScale.numerator+"/"+p->stc.certificate().outputScale.denominator;
        trace.details["outputScale.bits"]=std::to_string(p->stc.certificate().outputScale.binary64Bits);
        trace.details["headroom"]="All CtS/EvalRound/combination/StC stages checked at actual primes/scales";
        auto evidence=[&](BootstrapGate g,const BootstrapBound& v,const std::string& why) { trace.certificate.gates[static_cast<std::size_t>(g)]={true,why,{v}}; };
        evidence(BootstrapGate::Input,request.upstream.raisedMagnitude,"Resolved exact source and analytical upstream bounds");
        evidence(BootstrapGate::KeySwitch,experimental::finiteSupportBackendKeyNoise(),"Shared backend finite support");
        evidence(BootstrapGate::Domain,trace.bounds["rho_cert"],"Computed LP domain and mandatory lift gate");
        for(auto g:{BootstrapGate::Extraction,BootstrapGate::Reconstruction}) evidence(g,trace.bounds["E_I"],p->evalRound.certification().provenance);
        evidence(BootstrapGate::CoeffToSlotHP,trace.bounds["E_HP"],p->cts->hp().evidence.provenance);
        evidence(BootstrapGate::CoeffToSlotLP,trace.bounds["E_LP"],p->cts->lp().evidence.provenance);
        evidence(BootstrapGate::Combination,trace.bounds["E_pair"],"Actual certified gamma/subtract DAG");
        evidence(BootstrapGate::SlotToCoeff,trace.bounds["E_StC"],p->stc.certificate().result.provenance);
        evidence(BootstrapGate::Arithmetic,trace.bounds["E_boot"],"No double StC amplification");
        mpq_class utilization=0;
        auto headroomUse=[&](const std::vector<std::uint64_t>& primes,const std::string& margin) {
            mpz_class modulus=1; for(auto prime:primes) modulus*=integer(prime);
            const mpz_class exactMargin(margin);
            require(exactMargin>0,S::HeadroomViolation,"headroom","Missing positive compiled centered margin");
            mpq_class used=1-mpq_class(exactMargin)/modulus; if(used>utilization) utilization=used;
        };
        for(const auto* branch:{&p->cts->hp(),&p->cts->lp()}) for(const auto& factor:branch->certifiedFactors)
            for(const auto& state:factor.runtime) headroomUse(state.activePrimes,state.centeredHeadroomNumerator);
        for(const auto& node:p->evalRound.nodes()) headroomUse(node.activePrimes,node.centeredHeadroomNumerator);
        for(const auto& node:p->combination) headroomUse(node.activePrimes,node.centeredHeadroomNumerator);
        for(const auto& factor:p->stc.certificate().factors) for(const auto& state:factor.runtime) headroomUse(state.activePrimes,state.centeredHeadroomNumerator);
        trace.bounds["headroom.utilization"]=bound(utilization,"Maximum 2*ceil(scale*(M+E))/exact Q over compiled output states; internal BSGS QP checks separately verified");
        evidence(BootstrapGate::Headroom,trace.bounds["headroom.utilization"],trace.details["headroom"]);
        for(auto g:{BootstrapGate::ScaleSchedule,BootstrapGate::EvaluationKeys}) trace.certificate.gates[static_cast<std::size_t>(g)]={true,"All compiled actual schedules and required keys checked",{}};
        trace.certificate.minimumSecurityBits=a.securityLevelBits();
        trace.certificate.gates[static_cast<std::size_t>(BootstrapGate::Security)]={a.securityLevelBits()>=request.target.targetSecurityBits,"Actual SEAL context security level",{}};
        trace.certificate.gates[static_cast<std::size_t>(BootstrapGate::SparseSecret)]={false,"Ordinary ModRaise; sparse-secret production certificate is outside PR-5",{}};
        trace.certificate.outputError=trace.bounds["E_boot"];
        note({trace.bounds["E_boot"].upperBound<=request.target.targetAbsoluteError?S::Certified:S::ErrorBudgetExceeded,"accuracy","Actual end-to-end arithmetic bound compared with global target"});
        p->readiness={S::Certified,"execution","All executable stage certificates and schedules valid; global status remains separate"};
        if(request.lift.evidence==BootstrapLiftEvidence::TestFixtureAssumption) note({S::TestOnlyAssumption,"lift.fixture","External K=1 fixture assumption cannot produce production certification"});
        note(validateBootstrapCertificate(request.target,trace.plan,trace.certificate));
    } catch(const Rejection& f) { note(f.gate); p->readiness=f.gate; }
      catch(const experimental::arithmetic::Failure& f) { BootstrapContractResult g{f.status,"arithmetic",f.why}; note(g); p->readiness=g; }
      catch(const std::exception& e) { BootstrapContractResult g{S::InvalidInput,"planning",e.what()}; note(g); p->readiness=g; }
    if(trace.publicRlwe) trace.gateResults.push_back(trace.publicRlwe->result);
    trace.details["firstFailingGate"]=trace.result.gate;
    if(trace.result.status!=S::Certified) for(auto& entry:trace.bounds) if(entry.second.kind==BootstrapBoundKind::Unknown)
        entry.second.provenance="Unavailable after "+trace.result.gate+": "+trace.result.provenance;
    if(trace.certificate.outputError.kind==BootstrapBoundKind::Unknown) trace.certificate.outputError=trace.bounds["E_boot"];
    return plan;
}
}
