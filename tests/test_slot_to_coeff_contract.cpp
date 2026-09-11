#include "m2424/slot_to_coeff.hpp"
#include <cmath>
#include <cstdio>
#include <stdexcept>
using namespace m2424;
namespace {
void check(bool x,const char* p) { if(!x) throw std::runtime_error(p); }
BootstrapBound b(double x) { return {x,BootstrapBoundKind::Deterministic,"Analytical test bound",{}}; }
void contracts() {
    LinearTransformFactorBound f{b(2),b(.01),b(.003)};
    auto result=propagateLinearTransform(b(3),b(.1),f);
    check(result.result.status==BootstrapCertificationStatus::Certified&&result.semanticError.upperBound>=.234,"v9 recurrence includes operator perturbation separately");
    check(result.operatorPerturbationError.upperBound>=.031&&result.propagatedError.upperBound>=.2,"Separate propagation terms");
    for(int missing=0;missing<3;++missing) { auto no=f; (missing==0?no.kappa:missing==1?no.delta:no.localArithmeticError)={}; check(propagateLinearTransform(b(3),b(.1),no).result.status==BootstrapCertificationStatus::RequiredBoundUnavailable,"Unknown is never zero"); }
    auto gain=linearTransformGain({f,f}); check(gain.upperBound>=2.01*2.01,"Prepared factor gain");
    f.localArithmeticError={}; check(linearTransformGain({f}).kind!=BootstrapBoundKind::Unknown,"Additive local error not part of operator gain");
    PreparedSlotToCoeffPlan empty; check(empty.certificate().result.status!=BootstrapCertificationStatus::Certified,"Empty prepared plan is not an execution certificate");
}
#ifdef M2424_TEST_STC_BACKEND
void preflight() {
    SlotToCoeffPlan plan(8192,4);
    auto adapter=SealAdapter::create({8192,{40,40,40,40,40},std::ldexp(1.,30),1});
    adapter.generateKeys(plan.requirements().rotations,false);
    auto x=adapter.encrypt(adapter.encode({.001}));
    SlotToCoeffContract c; for(int i=0;i<2;++i) { c.inputMagnitude[i]=b(.001); c.inputError[i]=b(1e-8); }
    auto shallow=plan.prepare(adapter,x,x,c);
    check(shallow.certificate().result.status==BootstrapCertificationStatus::InsufficientLevels,"Insufficient levels rejected before encoding or execution");
    // Rejection contracts are exercised cheaply before any large diagonal preparation.
    SlotToCoeffPlan deeper(8192,3); c.inputError[0]={};
    check(deeper.prepare(adapter,x,x,c).certificate().result.status==BootstrapCertificationStatus::RequiredBoundUnavailable,"Unknown input prevents certification");
    c.inputError[0]=b(1e-8); c.evaluationKeyNoiseSupport={};
    check(deeper.prepare(adapter,x,x,c).certificate().result.status==BootstrapCertificationStatus::RequiredBoundUnavailable,"Unknown key noise prevents certification");
}
#endif
}
int main() { try { contracts();
#ifdef M2424_TEST_STC_BACKEND
preflight();
#endif
std::puts("SlotToCoeff contract tests passed"); } catch(const std::exception& e) { std::fprintf(stderr,"%s\n",e.what()); return 1; } }
