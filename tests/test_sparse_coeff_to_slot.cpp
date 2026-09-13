#include "m2424/m2424.hpp"
#include "bootstrap_fixture.hpp"
#include "sparse_bootstrap_oracle.hpp"
#include "../src/core/sparse_bootstrap_internal.hpp"
#include <gmpxx.h>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <cstring>
using namespace m2424;
using S=BootstrapCertificationStatus;
using O=test::SparseBootstrapOracle;
namespace m2424::test {
// Corrupt only a test copy of the immutable certificate to exercise fail-closed
// preparation. No production API accepts a caller restoration error override.
class SparseCoeffToSlotFixture {
public:
    static SparseBootstrapPlan removeRestoration(const SparseBootstrapPlan& original,bool zero) {
        SparseBootstrapPlan p; auto d=std::make_shared<SparseBootstrapPlan::Data>(*original.data_);
        d->certificate.restorationError=zero?BootstrapBound{0,BootstrapBoundKind::Deterministic,"Forged test-only zero restoration",{}}:BootstrapBound{};
        p.data_=d; return p;
    }
};
}
void check(bool okay,const std::string& why) { if(!okay) throw std::runtime_error(why); }
BootstrapBound B(const mpq_class& x,const std::string& why) {
    double d=x.get_d(); if(mpq_class(d)<x) d=std::nextafter(d,INFINITY);
    return {d,BootstrapBoundKind::Deterministic,why,{}};
}
double error(const ComplexVector& x,const ComplexVector& y) {
    double e=0; for(std::size_t i=0;i<x.size();++i)e=std::max(e,std::abs(x[i]-y[i])); return e;
}
int main() {
 try {
    std::cout<<std::setprecision(14);
    constexpr std::size_t N=16384,H=64; const double Delta=std::exp2(59.5);
    EvalRoundPlusCoeffToSlot cts(N); auto keys=cts.requirements().rotationSteps; keys.push_back(0);
    auto a=SealAdapter::create({N,std::vector<int>(7,60),Delta,N/2});
    a.generateKeys(keys,true); a.generateSparseBootstrapKeys(H);
    // Exact zero message; finite-support encryption noise certificate, no
    // observed value participates in the input envelope or K.
    auto top=a.encrypt(a.encode({0.})); auto input=a.modSwitchToChainIndex(top,0);
    SparseBootstrapInput request;
    request.messageMagnitude=B(0,"Exactly encoded zero polynomial");
    request.sourceNoiseMagnitude=B(mpq_class(21)*(2*N+1)+mpq_class(N+1)/2,"Asymmetric encryption finite coefficient support, conservatively ignoring prime reduction");
    auto sparse=prepareSparseBootstrap(a,input,request); const auto& sc=sparse.certificate();
    check(sc.result.status==S::Certified,"Sparse prerequisite: "+sc.result.provenance);
    check(sc.K==H&&sc.keys.weight==H,"K remains exactly h_b");
    auto encapsulated=a.encapsulateSparse(input); auto centered=O::source(a,encapsulated);
    auto raised=a.modRaiseSparse(encapsulated); const auto u=O::raised(a,raised);
    const auto y=coeffToSlotReference(u);
    CoeffToSlotContract schedule{"sparse-first-factor",N/2,N,59.5,59.5,.25,1e-10};
    for(bool zero:{false,true}) {
        auto forged=test::SparseCoeffToSlotFixture::removeRestoration(sparse,zero);
        check(cts.prepareCertified(a,raised,forged,schedule,schedule).hp().certificate.status==S::KeySwitchBoundUnavailable,"Unknown/zero restoration cannot certify");
    }
    const auto before=O::restorationCount(a);
    BootstrapRequest req; req.upstream.messageMagnitude=request.messageMagnitude;
    req.upstream.sourceNoiseMagnitude=request.sourceNoiseMagnitude;
    auto full=Bootstrapper(N).prepare(a,input,req);
    check(full.trace().result.status==S::UnsupportedEvalRoundDomain,"K>1 full production rejection");
    check(full.evalRound().nodes().empty(),"No EvalRound DAG for K>1");
    check(!full.trace().publicRlwe->minimumSecurityBits,"No invented security evidence");
    check(O::restorationCount(a)==before,"Preparation did not restore the sparse-raised input");
    check(full.coeffToSlot()!=nullptr,"Full planner retains the prepared first-factor CtS certificate");
    const auto& p=*full.coeffToSlot();
    check(p.hp().certificate.status==S::Certified,p.hp().certificate.provenance);
    check(p.lp().certificate.status==S::Certified,p.lp().certificate.provenance);
    check(cts.preflight(a,raised,p).status==S::Certified,"Sparse prepared preflight");
    // The two independently prepared sparse certificates describe identical
    // public material; moving restoration ownership must not change inventory.
    const auto& families=full.trace().publicRlwe->families;
    check(families.size()==sc.security.families.size(),"Nine families preserved");
    for(std::size_t i=0;i<families.size();++i) check(families[i].statement==sc.security.families[i].statement,"Security statement unchanged");
    const auto alpha=CoeffToSlotPrefactor::sourceNormalization(sc.input);
    std::vector<ComplexVector> refs[2][2];
    for(std::size_t b=0;b<2;++b)for(std::size_t h=0;h<2;++h) {
        refs[b][h]=cts.applyPlainTrace(y,h,b?alpha:CoeffToSlotPrefactor{});
        for(std::size_t j=0;j<N/2;++j) {
            const auto k=h*N/2+j; const double ideal=b?alpha.multiplyRounded(u[k]):u[k];
            check(std::abs(refs[b][h].back()[j]-ideal)<1e-10,"Plaintext factors versus independent coefficient halves");
        }
    }
    std::size_t nodes=0,restorationNodes=0; double maxRatio=0,observedRest=0;
    auto out=cts.apply(a,raised,p,[&](BootstrapGate gate,const SlotToCoeffRuntimeStage& stage,const Cipher& ct) {
        const bool restoration=stage.operation.find("first-factor sparse restoration")!=std::string::npos;
        const auto b=gate==BootstrapGate::CoeffToSlotLP;
        ComplexVector expected;
        if(restoration) { expected=y; ++restorationNodes; }
        else {
            expected=refs[b][stage.branch][stage.factor];
            if(stage.operation=="conjugation") {expected=refs[b][stage.branch][stage.factor-1];for(auto& v:expected)v=std::conj(v);}
            ++nodes;
        }
        auto e=error(a.decodeComplex(a.decrypt(ct)),expected);
        check(e<=stage.semanticError.upperBound,"Observed stage error exceeds certificate: "+stage.operation);
        check(stage.semanticError.kind==BootstrapBoundKind::Deterministic&&!stage.semanticError.provenance.empty()&&!stage.localError.provenance.empty(),"Known stage provenance");
        check(mpz_class(stage.centeredHeadroomNumerator)>0,"Stage centered headroom");
        maxRatio=std::max(maxRatio,e/stage.semanticError.upperBound); if(restoration)observedRest=e;
        std::cout<<(b?"LP":"HP")<<" h="<<stage.branch<<" r="<<stage.factor<<" "<<stage.operation<<" observed="<<e<<" certified="<<stage.semanticError.upperBound<<"\n";
    });
    check(O::restorationCount(a)==before+1&&restorationNodes==1&&out.restorationOperations==1,"Exactly one shared first-factor restoration");
    check(nodes==4*(2*cts.plan().depth()+2),"Every CtS kernel/rescale/projection stage checked");
    // The input object itself was not replaced by a main-secret RaisedCipher.
    check(O::raised(a,raised)==u,"SparseRaisedCipher stays under s_b and is immutable");
    double finalObserved[2]={}; const Cipher* outputs[]={&out.hpFirst,&out.hpSecond,&out.lpFirst,&out.lpSecond};
    for(std::size_t b=0;b<2;++b)for(std::size_t h=0;h<2;++h) {
        const auto& trace=b?p.lp():p.hp(); auto actual=a.decodeComplex(a.decrypt(*outputs[2*b+h]));
        BootstrapBound M=sc.raisedMagnitude,E=sc.restorationError;
        for(std::size_t r=0;r<=cts.plan().depth();++r) {
            const auto& f=trace.certifiedFactors[h*(cts.plan().depth()+1)+r];
            auto recurrence=propagateLinearTransform(M,E,f.bounds);
            check(recurrence.semanticError.upperBound==f.propagation.semanticError.upperBound,"Restoration propagated by actual recurrence, exactly once");
            if(r==0) {
                check(f.incomingSemanticError.upperBound==sc.restorationError.upperBound,"First-factor incoming error is E_rest");
                check(f.arithmeticTerms[0].first.find("incoming restoration")!=std::string::npos&&f.arithmeticTerms[1].first.find("incoming restoration")!=std::string::npos,"Restoration key-noise and ModDown provenance distinguished from B_local");
            }
            M=recurrence.idealMagnitude; E=recurrence.semanticError;
        }
        check(trace.inputSemanticError.upperBound==sc.restorationError.upperBound,"HP and LP share one incoming restoration envelope");
        check(trace.rescaleOperations==2*cts.plan().depth()&&trace.levelsConsumed==cts.plan().depth(),"Folded alpha without extra rescale");
        for(std::size_t j=0;j<N/2;++j) {
            const auto k=h*N/2+j;
            const double I=std::round(alpha.multiplyRounded(u[k]-centered[k]));
            check(std::abs(I)<=H,"Observed lift respects unchanged analytical K");
            const auto ideal=b?I+alpha.multiplyRounded(centered[k]):u[k];
            finalObserved[b]=std::max(finalObserved[b],std::abs(actual[j]-ideal));
        }
        check(finalObserved[b]<=trace.outputError.upperBound,"Final independent ideal-semantic oracle");
    }
    mpz_class qsrc=1;for(auto prime:sc.input.sourcePrimes)qsrc*=mpz_class(std::to_string(prime));
    const mpq_class rho=(mpq_class(sc.messageMagnitude.upperBound)+mpq_class(sc.sourceNoiseMagnitude.upperBound))/qsrc+mpq_class(p.lp().outputError.upperBound);
    check(p.domain().rho.kind==BootstrapBoundKind::Deterministic,"Known final rho certificate");
    check(mpq_class(p.domain().rho.upperBound)>=rho&&mpq_class(std::nextafter(p.domain().rho.upperBound,0.))<rho,"rho includes E_LP once, no separate E_rest term");
    check(full.trace().bounds.at("nu_b").upperBound==sc.sourceNoiseMagnitude.upperBound,"Restoration does not redefine nu_b");
    check(full.trace().bounds.at("K").upperBound==H,"Full planner retains K=h_b");
    std::cout<<"FINAL E_rest="<<sc.restorationError.upperBound<<" observed_rest="<<observedRest<<" E_HP="<<p.hp().outputError.upperBound<<" observed_HP="<<finalObserved[0]<<" E_LP="<<p.lp().outputError.upperBound<<" observed_LP="<<finalObserved[1]<<" rho_cert="<<p.domain().rho.upperBound<<" max_observed/certified="<<maxRatio<<" h_b=K="<<H<<"\n";
    O::removeKey(a,true); std::size_t callbacks=0;
    check(cts.preflight(a,raised,p).status==S::MissingEvaluationKeys,"Missing restoration key fails preflight");
    bool rejected=false;try {cts.apply(a,raised,p,[&](auto,const auto&,const auto&){++callbacks;});}catch(const std::invalid_argument&){rejected=true;}
    check(rejected&&callbacks==0&&O::restorationCount(a)==before+1,"Missing key executes no first-factor arithmetic");
    a.generateSparseBootstrapKeys(H);
    check(cts.preflight(a,raised,p).status==S::InvalidInput,"New sparse generation invalidates prepared plan");
    auto other=test::BootstrapFixture::create({16,std::vector<int>(5,50),std::ldexp(1.,49),8}); other.generateKeys(std::vector<int>{0});other.generateSparseBootstrapKeys(4);
    check(cts.preflight(other,raised,p).status==S::InvalidInput,"Changed context invalidates prepared plan");
    std::cout<<"PASS sparse first-factor CtS\n";
 }catch(const std::exception& e){std::cerr<<e.what()<<"\n";return 1;}
}
