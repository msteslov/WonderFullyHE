#include "m2424/m2424.hpp"
#include "m2424/evalround_plus_coeff_to_slot.hpp"
#include "evalround_test_support.hpp"
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <stdexcept>
using namespace m2424;
namespace {
void check(bool v,const char* why) { if(!v) throw std::runtime_error(why); }
double error(const ComplexVector& a,const ComplexVector& b) { double e=0; for(std::size_t i=0;i<a.size();++i) e=std::max(e,std::abs(a[i]-b[i])); return e; }
BootstrapBound bound(const mpq_class& v,const char* why) { double d=v.get_d(); if(mpq_class(d)<v) d=std::nextafter(d,INFINITY); return {d,BootstrapBoundKind::Deterministic,why,{}}; }
void plain() {
    BootstrapInputContext source; source.sourcePrimes={17,97}; double scale=1024; std::memcpy(&source.scaleBinary64Bits,&scale,8);
    auto alpha=CoeffToSlotPrefactor::sourceNormalization(source);
    for(std::size_t N:{8,16,32}) for(std::size_t depth:{1,2,4}) {
        EvalRoundPlusCoeffToSlot plan(N,depth);
        for(std::size_t column=0;column<N;++column) {
            std::vector<double> coefficients(N); coefficients[column]=1;
            auto y=coeffToSlotReference(coefficients);
            for(std::size_t h=0;h<2;++h) {
                auto hp=plan.applyPlainTrace(y,h).back(),lp=plan.applyPlainTrace(y,h,alpha).back();
                for(std::size_t j=0;j<N/2;++j) {
                    check(std::abs(hp[j]-coefficients[h*N/2+j])<1e-12,"Exact-root HP factors vs independent coefficient reference");
                    check(std::abs(lp[j]-coefficients[h*N/2+j]*1024/1649)<1e-12,"Exact-root phased LP factors vs independent reference");
                }
            }
        }
    }
}
void domainTests() {
    BootstrapInputContext source; source.sourcePrimes.assign(20,1152921504606830593ULL);
    CoeffToSlotCertificationInput input; input.messageMagnitude=bound(mpq_class(1e100),"Synthetic exact analytical magnitude");
    input.sourceNoiseMagnitude=bound(0,"Exact zero noise for rational domain unit test");
    auto zero=bound(0,"Exact zero LP error for rational domain unit test");
    auto result=certifyCoeffToSlotDomain(source,input,zero);
    mpz_class Q=1; for(auto prime:source.sourcePrimes) Q*=mpz_class(std::to_string(prime));
    mpq_class expected=mpq_class(input.messageMagnitude.upperBound)/Q;
    check(result.rho.kind==BootstrapBoundKind::Deterministic&&mpq_class(result.rho.upperBound)>=expected,
        "Domain uses exact modulus product beyond binary64 range");
    check(result.rho.upperBound>0&&mpq_class(std::nextafter(result.rho.upperBound,0.))<expected,"Tight outward rational domain bound");
    source.sourcePrimes={17}; input.messageMagnitude=bound(9,"Known excessive message magnitude");
    check(certifyCoeffToSlotDomain(source,input,zero).result.status==BootstrapCertificationStatus::DomainViolation,"Computed rho outside EvalRound domain rejected");
    check(certifyCoeffToSlotDomain(source,input,{}).rho.kind==BootstrapBoundKind::Unknown,"Unknown LP error is not zero");
}
void backend() {
    constexpr std::size_t N=16384,S=N/2;
    CkksProfile profile{N,std::vector<int>(7,60),std::exp2(59.5),S};
    EvalRoundPlusCoeffToSlot plan(N); auto keys=plan.requirements().rotationSteps; keys.push_back(0);
    auto a=SealAdapter::create(profile); a.generateKeys(keys,false);
    auto context=a.encrypt(a.encode({0.})); while(a.chainIndex(context)) context=a.rescaleToNext(context);
    // Exact scalar zero at the actual source modulus. Its encryption noise has
    // an independent finite-support bound, without a measured input certificate.
    auto input=evalround_test::encryptExactScalar(a,context,0,profile.scale);
    auto source=*resolveBootstrapInput(a,input).context; auto raised=a.modRaiseToTop(input);
    const auto alpha=CoeffToSlotPrefactor::sourceNormalization(source);
    mpz_class Q=1; for(auto prime:source.sourcePrimes) Q*=mpz_class(std::to_string(prime));
    CoeffToSlotCertificationInput contract;
    contract.raisedMagnitude=bound(mpq_class(N)*(N+1)*Q/(2*mpq_class(profile.scale)),"Centered source c0,c1 and ternary secret coefficient support one");
    contract.messageMagnitude=bound(0,"Exactly encoded zero message polynomial");
    contract.sourceNoiseMagnitude=bound(mpq_class(21)*(2*N+1)+mpq_class(N+1)/2,"Coefficient support of asymmetric encryption: 21(2N+1)+(N+1)/2, ignoring prime reduction conservatively");
    CoeffToSlotContract schedule{"certified_cts",S,N,59.5,59.5,.25,2e-10};
    auto unknown=contract; unknown.raisedMagnitude={};
    check(plan.prepareCertified(a,raised,source,schedule,schedule,unknown).hp().certificate.status==BootstrapCertificationStatus::RequiredBoundUnavailable,"Unknown raised magnitude is not zero");
    std::puts("Preparing HP/LP certified factors"); std::fflush(stdout);
    auto p=plan.prepareCertified(a,raised,source,schedule,schedule,contract);
    if(p.hp().certificate.status!=BootstrapCertificationStatus::Certified) throw std::runtime_error(p.hp().certificate.provenance);
    check(p.lp().evidence.verified&&p.hp().evidence.verified,"Both gate evidences verified");
    auto changed=source; ++changed.scaleBinary64Bits;
    check(plan.preflight(a,raised,changed,p).status!=BootstrapCertificationStatus::Certified,"Changed exact scale rejected");
    changed=source; --changed.sourcePrimes[0];
    check(plan.preflight(a,raised,changed,p).status!=BootstrapCertificationStatus::Certified,"Changed exact source rejected");
    changed=source; ++changed.contextFingerprint[0];
    check(plan.preflight(a,raised,changed,p).status!=BootstrapCertificationStatus::Certified,"Changed context rejected");
    EvalRoundPlusCoeffToSlot different(N,3);
    check(different.preflight(a,raised,source,p).status!=BootstrapCertificationStatus::Certified,"Changed factorization rejected");
    auto noKeys=SealAdapter::create(profile);
    check(plan.preflight(noKeys,raised,source,p).status==BootstrapCertificationStatus::MissingEvaluationKeys,"Missing rotation/conjugation keys preflight");
    noKeys.generateKeys(plan.requirements().rotationSteps,false);
    check(plan.preflight(noKeys,raised,source,p).status==BootstrapCertificationStatus::MissingEvaluationKeys,"Missing conjugation alone preflight");
    const auto u=a.decryptRaisedCoefficientsAtRaisedModulus(raised),centered=a.decryptRaisedCoefficientsAtSourceModulus(raised);
    const auto y=coeffToSlotReference(u);
    std::vector<ComplexVector> references[2][2];
    for(std::size_t b=0;b<2;++b) for(std::size_t h=0;h<2;++h) references[b][h]=plan.applyPlainTrace(y,h,b?alpha:CoeffToSlotPrefactor{});
    std::size_t nodes=0;
    auto output=plan.apply(a,raised,p,[&](BootstrapGate gate,const SlotToCoeffRuntimeStage& s,const Cipher& ct) {
        std::size_t b=gate==BootstrapGate::CoeffToSlotLP;
        auto expected=references[b][s.branch][s.factor];
        if(s.operation=="conjugation") { expected=references[b][s.branch][s.factor-1]; for(auto& v:expected) v=std::conj(v); }
        const double observed=error(a.decodeComplex(a.decrypt(ct)),expected);
        check(observed<=s.semanticError.upperBound,"Observed factor/stage error exceeds certified error");
        check(!s.semanticError.provenance.empty()&&!s.localError.provenance.empty()&&mpz_class(s.centeredHeadroomNumerator)>0,"Stage bounds/headroom provenance");
        std::printf("%s h=%zu r=%zu %s chain=%zu scale=%.17g observed=%.12g certified=%.12g\n",b?"LP":"HP",s.branch,s.factor,s.operation.c_str(),s.chainIndex,a.scale(ct),observed,s.semanticError.upperBound); ++nodes;
    });
    check(nodes==4*(2*plan.plan().depth()+2),"Every factor kernel/rescale and conjugation/add observed");
    double observed[2]={}; const Cipher* outputs[]={&output.hpFirst,&output.hpSecond,&output.lpFirst,&output.lpSecond};
    bool nonzeroInteger=false;
    for(std::size_t b=0;b<2;++b) for(std::size_t h=0;h<2;++h) {
        auto actual=a.decodeComplex(a.decrypt(*outputs[2*b+h]));
        for(std::size_t j=0;j<S;++j) {
            const auto k=h*S+j; double integer=std::round(alpha.multiplyRounded(u[k]-centered[k])); nonzeroInteger|=integer!=0;
            check(std::abs(centered[k]*profile.scale)<=contract.sourceNoiseMagnitude.upperBound,"Test-only source noise satisfies analytical upstream bound");
            if(b) check(std::abs(actual[j]-integer)<=p.domain().rho.upperBound,"Observed LP domain residual <= computed rho");
            const double target=b?integer+alpha.multiplyRounded(centered[k]):u[k];
            observed[b]=std::max(observed[b],std::abs(actual[j]-target));
        }
    }
    check(nonzeroInteger,"Nontrivial integer lift exercised");
    for(std::size_t b=0;b<2;++b) {
        const auto& t=b?p.lp():p.hp(); check(t.outputError.kind==BootstrapBoundKind::Deterministic,"Final bound known");
        check(observed[b]<=t.outputError.upperBound&&observed[b]<=2e-10,"Independent final oracle, unchanged CtS diagnostic tolerance");
        check(t.levelsConsumed==plan.plan().depth()&&t.rescaleOperations==2*plan.plan().depth(),"No separate alpha rescale");
        check(t.errorBudget.status==BootstrapCertificationStatus::ErrorBudgetExceeded,"Actual large bounds retain feasibility failure");
        for(const auto& f:t.certifiedFactors) {
            for(const auto* v:{&f.bounds.kappa,&f.bounds.delta,&f.bounds.localArithmeticError}) check(v->kind==BootstrapBoundKind::Deterministic&&!v->provenance.empty(),"All factor bounds known with provenance");
            if(f.factor==plan.plan().depth()) check(f.bounds.localArithmeticError.upperBound>0,"Conjugation key-switch error nonzero");
            for(const auto& term:f.arithmeticTerms) std::printf("TERM %s h=%zu r=%zu %s=%.12g\n",b?"LP":"HP",f.branch,f.factor,term.first.c_str(),term.second.upperBound);
            std::printf("%s h=%zu r=%zu kappa=%.17g delta=%.17g B_local=%.12g\n",b?"LP":"HP",f.branch,f.factor,f.bounds.kappa.upperBound,f.bounds.delta.upperBound,f.bounds.localArithmeticError.upperBound);
        }
        std::printf("FINAL %s certified=%.12g observed=%.12g ratio=%.12g levels=%zu rescales=%zu rotations=%zu\n",b?"LP":"HP",t.outputError.upperBound,observed[b],t.outputError.upperBound/observed[b],t.levelsConsumed,t.rescaleOperations,plan.plan().metrics().rotationsPerApply);
    }
    const auto& first=p.lp().certifiedFactors.front().bounds;
    mpq_class exactKappa=mpq_class(p.hp().certifiedFactors.front().bounds.kappa.upperBound)*mpq_class(profile.scale)/Q;
    check(mpq_class(first.kappa.upperBound)>=exactKappa&&mpq_class(std::nextafter(first.kappa.upperBound,0.))<exactKappa,"LP first kappa rounded once from exact dyadic / actual source modulus");
    check(p.domain().rho.kind==BootstrapBoundKind::Deterministic,"Computed rho known");
    mpq_class rho=(mpq_class(contract.messageMagnitude.upperBound)+mpq_class(contract.sourceNoiseMagnitude.upperBound))/Q+mpq_class(p.lp().outputError.upperBound);
    check(mpq_class(p.domain().rho.upperBound)>=rho,"Rho includes real upstream bounds and E_LP");
    std::printf("rho_cert=%.12g factorization=",p.domain().rho.upperBound); for(auto r:plan.plan().factorization().radices) std::printf("%zu,",r); std::puts("");
    // Unknown upstream remains a separate domain blocker while CtS is certified.
    contract.sourceNoiseMagnitude={};
    auto noUpstream=certifyCoeffToSlotDomain(source,contract,p.lp().outputError);
    check(noUpstream.rho.kind==BootstrapBoundKind::Unknown&&noUpstream.result.status==BootstrapCertificationStatus::RequiredBoundUnavailable,"Unknown upstream never deterministic zero");
}
}
int main() { try { plain(); domainTests(); backend(); } catch(const std::exception& e) { std::fprintf(stderr,"%s\n",e.what()); return 1; } }
