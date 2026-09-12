#include "m2424/m2424.hpp"
#include "m2424/experimental/evalmod_analysis/certified_diagonal.hpp"
#include "bootstrap_fixture.hpp"
#include <gmpxx.h>
#include <mpfr.h>
#include <cstdio>
#include <stdexcept>
#include <algorithm>
using namespace m2424;
namespace {
void check(bool yes,const char* why) { if(!yes) throw std::runtime_error(why); }
BootstrapBound bound(double x,const char* why) { return {x,BootstrapBoundKind::Deterministic,why,{}}; }
struct C { mpq_class r,i; };
using V=std::vector<C>;
C plus(C a,C b) { return {a.r+b.r,a.i+b.i}; }
C times(C a,C b) { return {a.r*b.r-a.i*b.i,a.r*b.i+a.i*b.r}; }
double real(const mpq_class& q) { mpfr_t x; mpfr_init2(x,256); mpfr_set_q(x,q.get_mpq_t(),MPFR_RNDN); double d=mpfr_get_d(x,MPFR_RNDN); mpfr_clear(x); return d; }
V evaluate(const EvalRoundExecutionNode& n,const std::vector<V>& values,const V& input) {
    if(n.operation==EvalRoundOperation::Input) return input;
    V out=values.at(n.inputs[0]);
    for(std::size_t j=0;j<out.size();++j) {
        auto& v=out[j];
        switch(n.operation) {
        case EvalRoundOperation::Multiply: v=times(v,values[n.inputs[1]][j]); break;
        case EvalRoundOperation::Add: v=plus(v,values[n.inputs[1]][j]); break;
        case EvalRoundOperation::Subtract: { auto b=values[n.inputs[1]][j]; v=plus(v,{-b.r,-b.i}); break; }
        case EvalRoundOperation::MultiplyPlain: v=times(v,{mpq_class(mpz_class(n.constantNumerator),mpz_class(n.constantDenominator)),0}); break;
        case EvalRoundOperation::AddPlain: v.r+=mpq_class(mpz_class(n.constantNumerator),mpz_class(n.constantDenominator)); break;
        case EvalRoundOperation::Conjugate: v.i=-v.i; break;
        default: break;
        }
    }
    return out;
}
double error(const ComplexVector& actual,const V& expected,const BootstrapBound& certificate) {
    double worst=0; check(certificate.kind==BootstrapBoundKind::Deterministic,"Known stage error required");
    for(std::size_t j=0;j<actual.size();++j) {
        mpq_class r=mpq_class(actual[j].real())-expected[j].r,i=mpq_class(actual[j].imag())-expected[j].i;
        mpq_class square=r*r+i*i,limit(certificate.upperBound);
        check(square<=limit*limit,"Observed stage error exceeds analytical certificate"); worst=std::max(worst,std::sqrt(real(square)));
    } return worst;
}
V rational(const ComplexVector& values) { V out; for(auto v:values) out.push_back({mpq_class(v.real()),mpq_class(v.imag())}); return out; }
void productionRejectsUnprovedLift() {
    auto profile=profiles::fast_demo_ckks(); auto a=SealAdapter::create(profile);
    Bootstrapper bootstrap(profile.polyModulusDegree); auto keys=bootstrap.rotationKeys(); keys.push_back(0); a.generateKeys(keys,true);
    auto input=a.encrypt(a.encode({.125})); BootstrapRequest request;
    auto plan=bootstrap.prepare(a,input,request);
    check(plan.trace().result.status==BootstrapCertificationStatus::LiftBoundUnavailable,"Ordinary production encryption/ModRaise has no implicit K=1 certificate");
    check(plan.evalRound().nodes().empty()&&!bootstrap.apply(a,input,plan).output,"Unproved production lift never executes EvalRound");
    check(plan.trace().bounds.at("E_boot").kind==BootstrapBoundKind::Unknown&&!plan.trace().bounds.at("E_boot").provenance.empty(),"Unreached global errors remain explicit Unknown");
}
void run() {
    const std::size_t N=16,S=N/2; CkksProfile profile{N,std::vector<int>(15,50),std::ldexp(1.,49),S};
    auto a=test::BootstrapFixture::create(profile); Bootstrapper bootstrap(N,2,2); auto requiredKeys=bootstrap.rotationKeys(); requiredKeys.push_back(0); a.generateKeys(requiredKeys,true);
    check(a.securityLevelBits()==0,"Tiny fixture cannot masquerade as secure production parameters");
    auto top=a.encrypt(a.encode({0.})); auto bottom=a.modSwitchToChainIndex(top,0);
    std::vector<double> coefficients(N); for(std::size_t i=0;i<N;++i) coefficients[i]=(int(i%5)-2)*std::ldexp(1.,30);
    std::vector<std::uint64_t> residues;
    auto primes=a.dataModulusValues(); primes.push_back(a.specialKeyModulusValue());
    for(auto prime:primes) for(auto c:coefficients) { mpz_class value(static_cast<long>(c)),mod(std::to_string(prime)),r; mpz_mod(r.get_mpz_t(),value.get_mpz_t(),mod.get_mpz_t()); residues.push_back(std::stoull(r.get_str())); }
    auto encoded=a.encodePolynomialRnsAtKeyScale(residues,profile.scale);
    auto sourcePlain=a.modSwitchPlainTo(encoded,bottom); auto input=test::BootstrapFixture::input(a,sourcePlain);
    BootstrapRequest request; request.target.targetAbsoluteError=1e-6;
    request.upstream.messageMagnitude=bound(std::ldexp(1.,31),"Exact test message coefficient bound");
    request.upstream.sourceNoiseMagnitude=bound(1,"Synthetic test c=(m,1): source noise is ternary s, coefficient support one");
    request.upstream.raisedMagnitude=bound(N*(std::ldexp(1.,31)+1)/profile.scale,"Synthetic fixture has no source wrap: N*(M_m+1)/Delta");
    request.lift={1,BootstrapLiftEvidence::TestFixtureAssumption,"Explicit external K=1 integration fixture assumption; never production evidence"};
    auto missing=request; missing.lift={}; auto blocked=bootstrap.prepare(a,input,missing);
    check(blocked.trace().result.status==BootstrapCertificationStatus::LiftBoundUnavailable&&blocked.evalRound().nodes().empty(),"Unknown K rejects before EvalRound compilation");
    missing=request; missing.lift.K=2; blocked=bootstrap.prepare(a,input,missing);
    check(blocked.trace().result.status==BootstrapCertificationStatus::UnsupportedEvalRoundDomain,"K>1 does not silently become K=1");
    missing=request; missing.upstream.sourceNoiseMagnitude={}; blocked=bootstrap.prepare(a,input,missing);
    check(blocked.trace().result.status==BootstrapCertificationStatus::RequiredBoundUnavailable,"Unknown upstream error rejected");
    missing=request; missing.outputChainIndex=11; blocked=bootstrap.prepare(a,input,missing);
    check(blocked.trace().result.status==BootstrapCertificationStatus::InsufficientLevels,"Insufficient levels rejected before arithmetic");
    missing=request; missing.upstream.messageMagnitude=bound(double(a.coeffModulusValues(input)[0]),"Known excessive domain magnitude"); blocked=bootstrap.prepare(a,input,missing);
    check(blocked.executionReadiness().status==BootstrapCertificationStatus::DomainViolation&&blocked.evalRound().nodes().empty(),"rho>=1/2 blocks EvalRound");
    missing=request; missing.target.targetAbsoluteError=1e-14; blocked=bootstrap.prepare(a,input,missing);
    check(blocked.trace().result.status==BootstrapCertificationStatus::ErrorBudgetExceeded&&blocked.evalRound().nodes().empty(),"Local accuracy failure cannot become global Certified");
    check(blocked.trace().result.gate.find("CoeffToSlot")!=std::string::npos,"First failed CtS branch budget preserved");
    auto rejected=bootstrap.apply(a,input,blocked);
    check(!rejected.output&&rejected.trace.result.gate==blocked.trace().result.gate,"Apply preserves first global failure instead of overwriting with later planning blocker");
    missing=request; missing.target.targetAbsoluteError=1e-10; blocked=bootstrap.prepare(a,input,missing);
    check(blocked.executionReadiness().status==BootstrapCertificationStatus::ErrorBudgetExceeded&&blocked.evalRound().nodes().empty(),"Nonpositive backward budget prevents EvalRound");
    auto plan=bootstrap.prepare(a,input,request);
    if(plan.executionReadiness().status!=BootstrapCertificationStatus::Certified) throw std::runtime_error(plan.executionReadiness().gate+": "+plan.executionReadiness().provenance);
    check(plan.trace().result.status==BootstrapCertificationStatus::TestOnlyAssumption,"Locally executable fixture never globally Certified");
    check(plan.trace().bounds.at("Gamma_StC").upperBound>N,"Actual prepared gain includes diagonal perturbation");
    const auto& trace=plan.trace();
    mpq_class residual(trace.details.at("backward.remainder.exact"));
    mpq_class gamma(trace.details.at("gamma.exact"));
    check(mpq_class(trace.bounds.at("requiredIntegerError").upperBound)<=residual/(mpq_class(trace.bounds.at("Gamma_StC").upperBound)*gamma),"Backward budget uses actual gain and exact gamma, rounded downward");
    auto noKeys=test::BootstrapFixture::create(profile);
    check(bootstrap.preflight(noKeys,input,plan).status==BootstrapCertificationStatus::MissingEvaluationKeys,"Missing keys reject in preflight");
    auto noResult=bootstrap.apply(noKeys,input,plan); check(!noResult.output,"No execution with missing keys");
    auto raised=a.modRaiseToTop(input); auto u=a.decryptRaisedCoefficientsAtRaisedModulus(raised);
    const auto y=coeffToSlotReference(u); EvalRoundPlusCoeffToSlot cts(N,2); auto source=*resolveBootstrapInput(a,input).context; auto alpha=CoeffToSlotPrefactor::sourceNormalization(source);
    std::vector<ComplexVector> ctsRefs[2][2];
    for(std::size_t b=0;b<2;++b) for(std::size_t h=0;h<2;++h) ctsRefs[b][h]=cts.applyPlainTrace(y,h,b?alpha:CoeffToSlotPrefactor{});
    std::vector<V> erRefs[2],combRefs[2]; V hp[2],lp[2]; ComplexVector useful[2];
    for(std::size_t h=0;h<2;++h) { erRefs[h].resize(plan.evalRound().nodes().size()); combRefs[h].resize(plan.combinationNodes().size()); for(std::size_t j=0;j<S;++j) { hp[h].push_back({mpq_class(u[h*S+j]),0}); lp[h].push_back({mpq_class(u[h*S+j])/gamma,0}); useful[h].push_back(u[h*S+j]); } }
    SlotToCoeffPlan stc(N,2); auto stc0=stc.applyPlainTrace(useful[0],0),stc1=stc.applyPlainTrace(useful[1],1); auto combined=stc.applyPlain(useful[0],useful[1]);
    std::size_t visited=0;
    auto output=bootstrap.apply(a,input,plan,[&](const std::string& stage,std::size_t h,std::size_t i,const Cipher& ct,const BootstrapBound& bound) {
        V expected;
        if(stage=="CtS.HP"||stage=="CtS.LP") { auto b=stage=="CtS.LP"; auto values=ctsRefs[b][h][i/2]; if(i/2==2&&i%2==0) { values=ctsRefs[b][h][1]; for(auto& v:values) v=std::conj(v); } expected=rational(values); }
        else if(stage=="EvalRound") { const auto& n=plan.evalRound().nodes()[i]; expected=erRefs[h][i]=evaluate(n,erRefs[h],lp[h]); }
        else if(stage=="Combination") { const auto& n=plan.combinationNodes()[i]; V source=i==0?hp[h]:V(S); expected=combRefs[h][i]=evaluate(n,combRefs[h],source); }
        else { expected=rational(h==2?combined:(h?stc1:stc0)[i/2]); }
        double observed=error(a.decodeComplex(a.decrypt(ct)),expected,bound);
        std::printf("%s h=%zu node=%zu chain=%zu scale=%.17g observed=%.6g certified=%.6g\n",stage.c_str(),h,i,a.chainIndex(ct),a.scale(ct),observed,bound.upperBound); ++visited;
    });
    check(output.output.has_value()&&visited==output.trace.nodes.size(),"All executed stages recorded in unified trace");
    std::vector<double> message=coefficients; for(auto& c:message) c/=profile.scale; auto reference=coeffToSlotReference(message);
    double observed=error(a.decodeComplex(a.decrypt(*output.output)),rational(reference),trace.bounds.at("E_boot"));
    check(observed<=request.target.targetAbsoluteError,"Fixture arithmetic accuracy target");
    check(output.trace.result.status!=BootstrapCertificationStatus::Certified,"Execution does not upgrade global certificate");
    check(a.chainIndex(*output.output)==request.outputChainIndex,"Actual final level");
    for(const auto& item:trace.bounds) std::printf("BOUND %s=%.17g\n",item.first.c_str(),item.second.upperBound);
    std::printf("FINAL levels=%zu output_scale=%.17g observed=%.12g first_global_gate=%s\n",trace.plan.levelsUsed,a.scale(*output.output),observed,trace.result.gate.c_str());
}
}
int main() { try { productionRejectsUnprovedLift(); run(); } catch(const std::exception& e) { std::fprintf(stderr,"%s\n",e.what()); return 1; } }
