#include "../core/bootstrap_internal.hpp"
#include <algorithm>
#include <cfenv>
#include <cstring>
#include <stdexcept>
namespace m2424 {
BootstrapPlan::BootstrapPlan():data_(std::make_shared<Data>()) {}
const BootstrapTrace& BootstrapPlan::trace() const { return data_->trace; }
const BootstrapContractResult& BootstrapPlan::executionReadiness() const { return data_->readiness; }
const EvalRoundExecutionPlan& BootstrapPlan::evalRound() const { return data_->evalRound; }
const std::vector<EvalRoundExecutionNode>& BootstrapPlan::combinationNodes() const { return data_->combination; }
Bootstrapper::Bootstrapper(std::size_t N,std::size_t c,std::size_t s):cts_(N,c),stc_(N,s) {}
std::vector<int> Bootstrapper::rotationKeys() const {
    auto keys=cts_.requirements().rotationSteps,other=stc_.requirements().rotations;
    keys.insert(keys.end(),other.begin(),other.end());
    std::sort(keys.begin(),keys.end()); keys.erase(std::unique(keys.begin(),keys.end()),keys.end()); return keys;
}
BootstrapContractResult Bootstrapper::preflight(const SealAdapter& a,const Cipher& input,const BootstrapPlan& p) const {
    using S=BootstrapCertificationStatus;
    if(p.executionReadiness().status!=S::Certified) return p.executionReadiness();
    const auto& d=*p.data_;
    if(std::fegetround()!=FE_TONEAREST) return {S::ScaleScheduleInfeasible,"schedule","Round-to-nearest required"};
    auto source=resolveBootstrapInput(a,input,d.trace.plan.input->scaleBinary64Bits);
    if(!source.context) return source.result;
    const auto& expected=*d.trace.plan.input; const auto& actual=*source.context;
    if(actual.contextFingerprint!=expected.contextFingerprint||actual.sourcePrimes!=expected.sourcePrimes||actual.raisedPrimes!=expected.raisedPrimes||actual.specialPrime!=expected.specialPrime||d.degree!=cts_.plan().polyModulusDegree()||d.ctsDepth!=cts_.plan().depth()||d.stcDepth!=stc_.metrics().depth)
        return {S::InvalidInput,"context","Prepared source/factorization mismatch"};
    if(a.info(input).ciphertextSize!=2) return {S::InvalidInput,"input","Two-component input required"};
    if(!a.hasRotationKeys(d.keys)||!a.hasConjugationKey()||!a.hasRelinKeys()) return {S::MissingEvaluationKeys,"keys","All CtS/StC rotations, conjugation and relinearization keys required"};
    return d.readiness;
}
BootstrapResult Bootstrapper::apply(SealAdapter& a,const Cipher& input,const BootstrapPlan& p,const BootstrapObserver& observer) const {
    BootstrapResult result; result.trace=p.trace();
    auto gate=preflight(a,input,p); if(gate.status!=BootstrapCertificationStatus::Certified) {
        if(p.executionReadiness().status==BootstrapCertificationStatus::Certified||result.trace.result.status==BootstrapCertificationStatus::Certified) result.trace.result=gate;
        result.trace.gateResults.push_back(gate); return result;
    }
    const auto& d=*p.data_;
    auto audit=[&](const std::string& stage,std::size_t branch,std::size_t index,const Cipher& ct,const BootstrapBound& M,const BootstrapBound& E,const BootstrapBound& local) {
        BootstrapTraceNode node; node.stage=stage+" half "+std::to_string(branch); node.node=index;
        node.state=*result.trace.plan.input; auto info=a.info(ct); std::memcpy(&node.state.scaleBinary64Bits,&info.scale,8);
        node.activePrimes=a.coeffModulusValues(ct); node.chainIndex=info.chainIndex;
        node.valueAbs=M; node.semanticError=E; node.localError=local; result.trace.nodes.push_back(std::move(node));
        if(observer) observer(stage,branch,index,ct,E);
    };
    auto raised=a.modRaiseToTop(input);
    auto halves=cts_.apply(a,raised,*d.cts,[&](BootstrapGate g,const SlotToCoeffRuntimeStage& s,const Cipher& ct) {
        audit(g==BootstrapGate::CoeffToSlotHP?"CtS.HP":"CtS.LP",s.branch,2*s.factor+(s.operation=="rescale"||s.operation=="real projection addition"),ct,s.magnitude,s.semanticError,s.localError);
    });
    Cipher hp[]={halves.hpFirst,halves.hpSecond},lp[]={halves.lpFirst,halves.lpSecond},pair[2];
    for(std::size_t h=0;h<2;++h) {
        auto integer=executeEvalRound(a,lp[h],d.evalRound,[&](std::size_t i,const Cipher& ct) { const auto& n=d.evalRound.nodes()[i]; audit("EvalRound",h,i,ct,n.idealMagnitude,n.semanticError,n.localArithmeticError); });
        pair[h]=executeCertifiedArithmetic(a,{hp[h],integer},d.combination,d.constants,d.combinationOutput,[&](std::size_t i,const Cipher& ct) { const auto& n=d.combination[i]; audit("Combination",h,i,ct,n.idealMagnitude,n.semanticError,n.localArithmeticError); });
    }
    auto output=stc_.apply(a,pair[0],pair[1],d.stc,[&](const SlotToCoeffRuntimeStage& s,const Cipher& ct) { audit("StC",s.branch,2*s.factor+(s.operation=="rescale"),ct,s.magnitude,s.semanticError,s.localError); });
    result.output=std::move(output.ciphertext); result.trace.details["executed"]="true"; return result;
}
}
