#include "m2424/bootstrap.hpp"
#include "m2424/canonical_embedding_reference.hpp"
#include "m2424/experimental/evalmod_analysis/finite_support_arithmetic.hpp"
#include "bootstrap_fixture.hpp"
#include "sparse_bootstrap_oracle.hpp"
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <random>
using namespace m2424;
using S=BootstrapCertificationStatus;
using O=m2424::test::SparseBootstrapOracle;
void check(bool b,const char* why) { if(!b) throw std::runtime_error(why); }
BootstrapBound B(double x) { return {x,BootstrapBoundKind::Deterministic,"Test analytical coefficient envelope",{}}; }
int main() {
 try {
    std::cout<<std::setprecision(14);
    CkksProfile profile{16,std::vector<int>(5,50),std::ldexp(1.,49),8};
    auto a=test::BootstrapFixture::create(profile); a.generateKeys({0,1,-1});
    std::mt19937 random(6021);
    double encObserved=0,restoreObserved=0,liftObserved=0;
    SparseBootstrapPlan last; Cipher input;
    for(std::size_t h:{2u,4u,8u,16u}) for(int repetition=0;repetition<4;++repetition) {
        a.generateSparseBootstrapKeys(h);
        auto secret=O::secret(a);
        check(std::count_if(secret.begin(),secret.end(),[](int x){return x!=0;})==h,"fixed weight");
        for(auto x:secret) check(x>=-1&&x<=1,"alphabet");
        std::vector<double> values(8); for(auto& x:values) x=std::uniform_real_distribution<double>(-.001,.001)(random);
        if(repetition==0) std::fill(values.begin(),values.end(),0.);
        auto top=a.encrypt(a.encode(values)); input=a.modSwitchToChainIndex(top,0);
        SparseBootstrapInput request; request.messageMagnitude=B(std::ldexp(1.,40)); request.sourceNoiseMagnitude=B(10000);
        auto plan=prepareSparseBootstrap(a,input,request); const auto& c=plan.certificate();
        check(c.result.status==S::Certified,"sparse arithmetic certification");
        check(c.K==h&&c.K>1,"analytical K comes from weight");
        check(!c.security.minimumSecurityBits&&c.security.result.status==S::SecurityBudgetExceeded,"no invented security");
        check(c.security.families.size()==9,"full public family inventory");
        for(const auto* id:{"ordinary.public","ordinary.encryption/0","ordinary.relinearization","ordinary.galois","sparse.encapsulation","sparse.restoration","derived.encapsulated","derived.raised","derived.restored"})
            check(std::any_of(c.security.families.begin(),c.security.families.end(),[&](const auto& f){return f.id==id;}),"required public family present");
        if(h==2) check(*c.security.families[4].searchSpaceCeilingBits<128,"low Hamming search space explicit");
        for(const auto& f:c.security.families) check(!f.provenance.empty()&&!f.modulus.empty()&&!f.statement.empty(),"family provenance/modulus");
        auto before=O::original(a,input); auto switched=a.encapsulateSparse(input);
        auto source=O::source(a,switched); auto raised=a.modRaiseSparse(switched); auto up=O::raised(a,raised);
        auto restored=a.restoreSparse(raised); auto final=a.decryptRaisedCoefficientsAtRaisedModulus(restored);
        const double Delta=a.scale(input),q=static_cast<double>(a.coeffModulusValues(input)[0]);
        for(std::size_t j=0;j<before.size();++j) {
            const auto e=std::abs(source[j]-before[j]),r=std::abs(final[j]-up[j]);
            check(e<=c.encapsulationError.upperBound,"encapsulation observed <= certificate");
            check(r<=c.restorationError.upperBound,"restoration observed <= certificate");
            const double I=std::round((up[j]-source[j])*Delta/q);
            check(std::abs(I)<=c.K,"lift <= analytical K");
            check(std::abs((up[j]-source[j])-q/Delta*I)<1e-12,"exact ModRaise semantic equation");
            encObserved=std::max(encObserved,e); restoreObserved=std::max(restoreObserved,r); liftObserved=std::max(liftObserved,std::abs(I));
        }
        auto composed=executeSparseBootstrap(a,input,plan);
        check(a.info(composed).chainIndex==3&&a.info(composed).scale==Delta,"no rescale or scale rewrite");
        auto unknown=request; unknown.evaluationKeyNoiseSupport=BootstrapBound{};
        check(prepareSparseBootstrap(a,input,unknown).certificate().result.status==S::KeySwitchBoundUnavailable,"unknown KS blocks");
        unknown=request; unknown.sourceNoiseMagnitude={};
        check(prepareSparseBootstrap(a,input,unknown).certificate().result.status==S::RequiredBoundUnavailable,"unknown upstream blocks");
        unknown=request; unknown.evaluationKeyNoiseSupport=B(0);
        check(prepareSparseBootstrap(a,input,unknown).certificate().result.status==S::KeySwitchBoundUnavailable,"fake zero KS blocks");
        auto changed=a.normalizeScale(input,Delta*1.001);
        check(preflightSparseBootstrap(a,changed,plan).status==S::InputScaleMismatch,"scale binding");
        last=plan;
    }
    const auto& c=last.certificate();
    std::cout<<"tiny oracle max encapsulation="<<encObserved<<" restoration="<<restoreObserved<<" max observed lift="<<liftObserved<<" K="<<c.K<<" bounds="<<c.encapsulationError.upperBound<<","<<c.restorationError.upperBound<<"\n";
    a.generateSparseBootstrapKeys(4);
    check(preflightSparseBootstrap(a,input,last).status==S::InvalidInput,"key generation invalidates plan");
    SparseBootstrapInput request; request.messageMagnitude=B(std::ldexp(1.,40)); request.sourceNoiseMagnitude=B(10000);
    auto plan=prepareSparseBootstrap(a,input,request);
    O::removeKey(a,true); check(preflightSparseBootstrap(a,input,plan).status==S::MissingEvaluationKeys,"missing restoration preflight");
    a.generateSparseBootstrapKeys(4); plan=prepareSparseBootstrap(a,input,request);
    O::removeKey(a,false); check(preflightSparseBootstrap(a,input,plan).status==S::MissingEvaluationKeys,"missing encapsulation preflight");
    bool refused=false; try { a.generateSparseBootstrapKeys(1); } catch(const std::invalid_argument&) {refused=true;} check(refused,"no h=1 shortcut");

    // Real tc128-validated ordinary context. Sparse h=64 is a candidate, not
    // claimed secure: its concrete all-family security remains Unknown.
    auto production=SealAdapter::create({4096,{36,36,36},std::ldexp(1.,35),2048});
    production.generateKeys({0,1},true); production.generateSparseBootstrapKeys(64);
    auto top=production.encrypt(production.encode({.001})); auto ct=production.modSwitchToChainIndex(top,0);
    SparseBootstrapInput real; real.messageMagnitude=B(std::ldexp(1.,28)); real.sourceNoiseMagnitude=B(10000);
    auto prepared=prepareSparseBootstrap(production,ct,real); const auto& pc=prepared.certificate();
    check(pc.result.status==S::Certified&&pc.K==64,"production analytical sparse lift");
    check(!pc.security.minimumSecurityBits&&production.securityLevelBits()==128,"tc128 not inherited by sparse");
    auto original=O::original(production,ct);
    auto encapsulated=production.encapsulateSparse(ct); auto sparseValues=O::source(production,encapsulated);
    auto sparseRaised=production.modRaiseSparse(encapsulated); auto lifted=O::raised(production,sparseRaised);
    auto result=production.restoreSparse(sparseRaised); auto restored=production.decryptRaisedCoefficientsAtRaisedModulus(result);
    double pe=0,pr=0,pI=0;
    const double qsrc=static_cast<double>(production.coeffModulusValues(ct)[0]),scale=production.scale(ct);
    for(std::size_t i=0;i<original.size();++i) {
        pe=std::max(pe,std::abs(original[i]-sparseValues[i]));
        pr=std::max(pr,std::abs(lifted[i]-restored[i]));
        auto I=std::round((lifted[i]-sparseValues[i])*scale/qsrc); pI=std::max(pI,std::abs(I));
        check(std::abs(lifted[i]-sparseValues[i]-I*qsrc/scale)<1e-12,"production ModRaise identity");
    }
    check(pe<=pc.encapsulationError.upperBound&&pr<=pc.restorationError.upperBound&&pI<=pc.K,"production oracle bounds");
    std::vector<double> encDifference(original.size()),restoreDifference(original.size());
    for(std::size_t i=0;i<original.size();++i) { encDifference[i]=sparseValues[i]-original[i]; restoreDifference[i]=restored[i]-lifted[i]; }
    double canonicalEnc=0,canonicalRestore=0;
    for(auto z:coeffToSlotReference(encDifference)) canonicalEnc=std::max(canonicalEnc,std::abs(z));
    for(auto z:coeffToSlotReference(restoreDifference)) canonicalRestore=std::max(canonicalRestore,std::abs(z));
    check(canonicalEnc<=pc.encapsulationError.upperBound&&canonicalRestore<=pc.restorationError.upperBound,"canonical switching errors <= certificates");
    std::cout<<"production observed canonical enc="<<canonicalEnc<<" restoration="<<canonicalRestore<<"\n";
    std::cout<<"production observed enc="<<pe<<" restoration="<<pr<<" lift="<<pI<<"\n";
    BootstrapRequest req; req.upstream.messageMagnitude=real.messageMagnitude; req.upstream.sourceNoiseMagnitude=real.sourceNoiseMagnitude;
    req.lift={1,BootstrapLiftEvidence::Analytical,"Untrusted caller K=1"};
    auto full=Bootstrapper(4096).prepare(production,ct,req);
    check(full.trace().result.status==S::UnsupportedEvalRoundDomain,"K>1 preflight blocker");
    check(full.trace().bounds.at("K").upperBound==64,"external K ignored");
    check(full.trace().bounds.at("rho_cert").kind==BootstrapBoundKind::Unknown,"missing CtS rho not zero");
    check(full.trace().publicRlwe&&!full.trace().publicRlwe->minimumSecurityBits,"full trace security");
    std::cout<<"production h="<<pc.keys.weight<<" K="<<pc.K<<" enc="<<pc.encapsulationError.upperBound<<" restore="<<pc.restorationError.upperBound<<" nu_b="<<pc.sourceNoiseMagnitude.upperBound<<" rho_before_CtS="<<pc.rhoBeforeCoeffToSlot.upperBound<<"\n";
    for(const auto& f:pc.security.families) { std::cout<<f.id<<" N="<<f.degree<<" samples="<<f.samples<<" components="<<f.components<<" modulus="; for(auto p:f.modulus)std::cout<<p<<","; std::cout<<" lambda=Unknown search-ceiling="<<*f.searchSpaceCeilingBits<<"\n"; }
    check(preflightSparseBootstrap(production,ct,plan).status==S::InvalidInput,"context binding");
    production.encrypt(production.encode({0.}));
    check(preflightSparseBootstrap(production,ct,prepared).status==S::InvalidInput,"new public samples invalidate prepared security statement");
    auto ordinary=SealAdapter::create({4096,{36,36,36},std::ldexp(1.,35),2048});
    ordinary.generateKeys({0}); auto ordinaryInput=ordinary.modSwitchToChainIndex(ordinary.encrypt(ordinary.encode({0.})),0);
    check(Bootstrapper(4096).prepare(ordinary,ordinaryInput,req).trace().result.status==S::LiftBoundUnavailable,"caller analytical K cannot certify ordinary production input");
    std::cout<<"PASS sparse bootstrap\n";
 } catch(const std::exception& e) { std::cerr<<e.what()<<"\n"; return 1; }
}
