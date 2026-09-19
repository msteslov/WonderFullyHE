#include "../core/certified_arithmetic_internal.hpp"
#include "../core/evalround_execution_internal.hpp"
#include <cstring>
#include <stdexcept>
#include <cfenv>
namespace m2424 {
namespace {
std::uint64_t bits(double x) { std::uint64_t b; std::memcpy(&b,&x,sizeof b); return b; }
BootstrapContractResult fail(BootstrapCertificationStatus s, const char* why) { return {s,"EvalRoundExecution",why}; }
}
EvalRoundExecutionPlan::EvalRoundExecutionPlan():data_(std::make_shared<Data>()) {}
const BootstrapContractResult& EvalRoundExecutionPlan::certification() const { return data_->certification; }
const std::vector<EvalRoundExecutionNode>& EvalRoundExecutionPlan::nodes() const { return data_->nodes; }
const EvalRoundPlan& EvalRoundExecutionPlan::mathematicalPlan() const { return data_->mathematicalPlan; }
double EvalRoundExecutionPlan::integerErrorUpper() const { return data_->mathematicalPlan.integerErrorUpper; }
std::size_t EvalRoundExecutionPlan::outputNode() const { return data_->output; }
const EvalRoundExecutionDiagnostics& EvalRoundExecutionPlan::diagnostics() const {
    return data_->diagnostics;
}
BootstrapContractResult preflightEvalRound(const SealAdapter& a, const Cipher& input, const EvalRoundExecutionPlan& p) {
    if(p.certification().status!=BootstrapCertificationStatus::Certified) return p.certification();
    if(std::fegetround()!=FE_TONEAREST) return fail(BootstrapCertificationStatus::ScaleScheduleInfeasible,"Binary64 schedule requires round-to-nearest");
    if(a.contextFingerprint()!=p.data_->fingerprint) return fail(BootstrapCertificationStatus::InvalidInput,"Context fingerprint mismatch");
    const auto info=a.info(input); const auto& first=p.nodes().front();
    if(bits(info.scale)!=first.outputScale.binary64Bits) return fail(BootstrapCertificationStatus::InputScaleMismatch,"Exact input scale bits mismatch");
    if(a.coeffModulusValues(input)!=first.activePrimes || info.chainIndex!=first.chainIndex)
        return fail(BootstrapCertificationStatus::InsufficientLevels,"Input active modulus differs from compiled schedule");
    if(info.ciphertextSize!=2) return fail(BootstrapCertificationStatus::InvalidInput,"Input must have two components");
    for(const auto& n:p.nodes()) {
        if(n.requiredKey==EvalRoundEvaluationKey::Relinearization && !a.hasRelinKeys())
            return fail(BootstrapCertificationStatus::MissingEvaluationKeys,"Missing relinearization key");
        if(n.requiredKey==EvalRoundEvaluationKey::Conjugation && !a.hasConjugationKey())
            return fail(BootstrapCertificationStatus::MissingEvaluationKeys,"Missing conjugation key");
    }
    return p.certification();
}
Cipher executeEvalRound(SealAdapter& a,const Cipher& input,const EvalRoundExecutionPlan& p,
    const std::function<void(std::size_t,const Cipher&)>& observer) {
    const auto gate=preflightEvalRound(a,input,p);
    if(gate.status!=BootstrapCertificationStatus::Certified) throw std::invalid_argument(gate.provenance);
    return executeCertifiedArithmetic(a,{input},p.nodes(),p.data_->constants,p.outputNode(),observer);
}
Cipher executeCertifiedArithmetic(SealAdapter& a,const std::vector<Cipher>& inputs,const std::vector<EvalRoundExecutionNode>& nodes,const std::vector<Plain>& constants,std::size_t outputNode,const std::function<void(std::size_t,const Cipher&)>& observer) {
    std::size_t inputIndex=0;
    std::vector<std::unique_ptr<Cipher>> values(nodes.size());
    std::vector<std::size_t> uses(values.size());
    for(const auto& n:nodes) for(auto i:n.inputs) ++uses[i];
    for(std::size_t i=0;i<values.size();++i) {
        const auto& n=nodes[i];
        auto x=[&](std::size_t j)->const Cipher& { return *values.at(n.inputs.at(j)); };
        Cipher out;
        switch(n.operation) {
        case EvalRoundOperation::Input: out=inputs.at(inputIndex++); break;
        case EvalRoundOperation::Multiply: out=a.multiply(x(0),x(1)); break;
        case EvalRoundOperation::MultiplyPlain: out=a.multiplyPlain(x(0),constants[i]); break;
        case EvalRoundOperation::Add: out=a.add(x(0),x(1)); break;
        case EvalRoundOperation::Subtract: out=a.sub(x(0),x(1)); break;
        case EvalRoundOperation::AddPlain: out=a.addPlain(x(0),constants[i]); break;
        case EvalRoundOperation::Relinearize: out=a.relinearize(x(0)); break;
        case EvalRoundOperation::Rescale: out=a.rescaleToNext(x(0)); break;
        case EvalRoundOperation::ModSwitch: out=n.inputs.size()==1?a.modSwitchToChainIndex(x(0),n.chainIndex):a.modSwitchTo(x(0),x(1)); break;
        case EvalRoundOperation::Conjugate: out=a.conjugate(x(0)); break;
        }
        const auto info=a.info(out);
        if(bits(info.scale)!=n.outputScale.binary64Bits || info.chainIndex!=n.chainIndex || info.ciphertextSize!=n.ciphertextComponents || a.coeffModulusValues(out)!=n.activePrimes)
            throw std::runtime_error("EvalRound runtime departed from certified scale/modulus schedule");
        if(observer) observer(i,out);
        values[i]=std::make_unique<Cipher>(std::move(out));
        for(auto j:n.inputs) if(--uses[j]==0) values[j].reset();
    }
    return std::move(*values.at(outputNode));
}
EvalRoundPairResult executeEvalRoundPair(SealAdapter& a,const Cipher& x,const Cipher& y,const EvalRoundExecutionPlan& p) {
    for(const auto* input:{&x,&y}) {
        const auto gate=preflightEvalRound(a,*input,p);
        if(gate.status!=BootstrapCertificationStatus::Certified) throw std::invalid_argument(gate.provenance);
    }
    return {executeEvalRound(a,x,p),executeEvalRound(a,y,p),p.integerErrorUpper()};
}
}
