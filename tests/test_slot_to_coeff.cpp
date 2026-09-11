#include "m2424/slot_to_coeff.hpp"
#include "m2424/coeff_to_slot.hpp"
#include <algorithm>
#include <cmath>
#include <cfenv>
#include <cstdio>
#include <stdexcept>
#ifdef M2424_TEST_STC_BACKEND
#include "m2424/experimental/evalmod_analysis/certified_diagonal.hpp"
#include "evalround_test_support.hpp"
#endif
using namespace m2424;
namespace {
void check(bool x,const char* message) { if(!x) throw std::runtime_error(message); }
double error(const ComplexVector& a,const ComplexVector& b) { double e=0; for(std::size_t i=0;i<a.size();++i) e=std::max(e,std::abs(a[i]-b[i])); return e; }
ComplexVector existingForwardOracle(const ComplexVector& x,const ComplexVector& y) {
    std::vector<double> real,imag;
    for(const auto* v:{&x,&y}) for(auto c:*v) { real.push_back(c.real()); imag.push_back(c.imag()); }
    auto out=coeffToSlotReference(real),second=coeffToSlotReference(imag);
    for(std::size_t i=0;i<out.size();++i) out[i]+=Complex(0,1)*second[i]; return out;
}
void plaintext() {
    for(std::size_t N:{8,16,32}) for(std::size_t depth:{1,2,4}) {
        SlotToCoeffPlan plan(N,depth); auto S=N/2;
        for(std::size_t column=0;column<N;++column) {
            ComplexVector x(S),y(S); (column<S?x[column]:y[column-S])=1;
            auto actual=plan.applyPlain(x,y),reference=existingForwardOracle(x,y);
            check(error(actual,reference)<2e-12,"Separate H and D_plus H basis columns");
            auto coefficients=slotToCoeffReference(actual);
            for(std::size_t k=0;k<N;++k) check(std::abs(coefficients[k]-(k==column?1.:0.))<2e-12,"Existing inverse oracle recovers basis coefficient");
        }
        ComplexVector x(S),y(S),slots(S);
        for(std::size_t j=0;j<S;++j) { x[j]={.03*double(j)-.2,.02*double(j)}; y[j]={.1-.05*double(j),-.03}; slots[j]={std::sin(double(j)),std::cos(double(j))}; }
        check(error(plan.applyPlain(x,y),existingForwardOracle(x,y))<2e-12,"Mixed complex halves");
        CoeffToSlotPlan cts(N,depth); const auto halves=cts.applyPlain(slots);
        check(error(plan.applyPlain(halves.first,halves.second),slots)<2e-12,"StC(CtS(a)) roundtrip independent of EvalRound");
    }
}
#ifdef M2424_TEST_STC_BACKEND
void backend() {
    const std::size_t N=16384,S=N/2; SlotToCoeffPlan plan(N,4);
    CkksProfile profile{N,std::vector<int>(7,60),std::ldexp(1.,55),S};
    auto a=SealAdapter::create(profile); a.generateKeys(plan.requirements().rotations,false);
    check(!a.hasRelinKeys()&&!plan.requirements().relinearization,"Linear baseline needs no unused relin key");
    ComplexVector x(S),y(S);
    // Sparse coefficient halves keep the independent O(N^2) oracle meaningful.
    const double amplitude=1./1024;
    x[0]=amplitude; x[3]=-amplitude; y[1]=amplitude; y[7]=-amplitude;
    auto context=a.encrypt(a.encode({0.}));
    experimental::CertifiedRootDiagonalEncoder inputEncoder(N);
    std::vector<int> roots0(S,-1),roots1(S,-1);
    roots0[0]=0; roots0[3]=int(N); roots1[1]=0; roots1[7]=int(N);
    auto encoded0=inputEncoder.encode(a,roots0,profile.scale,mpq_class(amplitude));
    auto encoded1=inputEncoder.encode(a,roots1,profile.scale,mpq_class(amplitude));
    auto plain0=a.modSwitchPlainTo(encoded0.plaintext,context),plain1=a.modSwitchPlainTo(encoded1.plaintext,context);
    check(error(a.decodeComplex(plain0),x)<=encoded0.perturbation.upperBound,"Observed encoded input perturbation <= analytical encoder bound");
    check(error(a.decodeComplex(plain1),y)<=encoded1.perturbation.upperBound,"Second encoded perturbation bound");
    auto first=a.encrypt(plain0),second=a.encrypt(plain1);
    SlotToCoeffContract contract;
    for(std::size_t b=0;b<2;++b) {
        contract.inputMagnitude[b]={amplitude,BootstrapBoundKind::Deterministic,"Exact maximum input coefficient-half magnitude",{}};
        const auto encryption=evalround_test::inputBound(profile);
        contract.inputError[b]={std::nextafter(encryption.upperBound+(b?encoded1:encoded0).perturbation.upperBound,INFINITY),
            BootstrapBoundKind::Deterministic,"Certified root-diagonal input encoding plus analytic ternary/CBD encryption and ModDown bound",{}};
    }
    check(error(a.decodeComplex(a.decrypt(first)),x)<=contract.inputError[0].upperBound,"Input error contract first");
    check(error(a.decodeComplex(a.decrypt(second)),y)<=contract.inputError[1].upperBound,"Input error contract second");
    std::puts("Preparing certified SlotToCoeff diagonals..."); std::fflush(stdout);
    auto prepared=plan.prepare(a,first,second,contract); const auto& cert=prepared.certificate();
    if(cert.result.status!=BootstrapCertificationStatus::Certified) throw std::runtime_error(cert.result.provenance);
    check(plan.preflight(a,first,second,prepared).status==BootstrapCertificationStatus::Certified,"Prepared preflight");
    auto rewritten=a.normalizeScale(first,std::nextafter(a.scale(first),INFINITY));
    check(plan.preflight(a,rewritten,second,prepared).status==BootstrapCertificationStatus::InputScaleMismatch,"First input scale fingerprint");
    check(plan.preflight(a,first,rewritten,prepared).status==BootstrapCertificationStatus::InputScaleMismatch,"Second input scale fingerprint");
    SlotToCoeffPlan changed(N,SlotToCoeffFactorization{{4,5,5,5}});
    check(changed.preflight(a,first,second,prepared).status==BootstrapCertificationStatus::InvalidInput,"Changed factorization invalidates preparation");
    auto changedProfile=profile; changedProfile.coeffModulusBits.back()=59;
    auto otherContext=SealAdapter::create(changedProfile);
    check(plan.preflight(otherContext,first,second,prepared).status==BootstrapCertificationStatus::InvalidInput,"Changed context fingerprint");
    auto noKeys=SealAdapter::create(profile);
    check(plan.preflight(noKeys,first,second,prepared).status==BootstrapCertificationStatus::MissingEvaluationKeys,"Missing rotation preflight");
    bool executed=false;
    try { plan.apply(noKeys,first,second,prepared,[&](const SlotToCoeffRuntimeStage&,const Cipher&){executed=true;}); throw std::runtime_error("Accepted missing keys"); }
    catch(const std::invalid_argument&) {}
    check(!executed,"Missing keys rejected before arithmetic");
    auto lowered=a.rescaleToNext(first);
    auto wrongLevel=a.encrypt(a.encodeScalarAtScaleFor(.001,profile.scale,lowered));
    check(plan.preflight(a,wrongLevel,second,prepared).status==BootstrapCertificationStatus::InsufficientLevels,"Changed input parms/active primes");
    auto needsRelin=a.multiply(first,first);
    auto missingRelin=plan.prepare(a,needsRelin,needsRelin,contract);
    check(missingRelin.certificate().keys.relinearization&&missingRelin.certificate().result.status==BootstrapCertificationStatus::MissingEvaluationKeys,"Missing required input relinearization rejects before encoding/execution");
    auto impossible=contract; impossible.inputMagnitude[0].upperBound=1e100;
    check(plan.prepare(a,first,second,impossible).certificate().result.status==BootstrapCertificationStatus::HeadroomViolation,"Headroom violation rejects before linear execution");
    std::fesetround(FE_UPWARD); auto roundingGate=plan.preflight(a,first,second,prepared); std::fesetround(FE_TONEAREST);
    check(roundingGate.status==BootstrapCertificationStatus::ScaleScheduleInfeasible,"Rounding-mode schedule binding");
    auto references0=plan.applyPlainTrace(x,0),references1=plan.applyPlainTrace(y,1);
    const auto reference=existingForwardOracle(x,y);
    std::size_t observedNodes=0;
    auto output=plan.apply(a,first,second,prepared,[&](const SlotToCoeffRuntimeStage& stage,const Cipher& cipher) {
        const auto& expected=stage.branch==2?reference:(stage.branch?references1:references0)[stage.factor];
        const double e=error(a.decodeComplex(a.decrypt(cipher)),expected);
        check(e<=stage.semanticError.upperBound,"Observed stage error exceeds analytical certificate");
        check(mpz_class(stage.centeredHeadroomNumerator)>0&&!stage.headroomProvenance.empty(),"Positive exact headroom with provenance");
        check(mpq_class(mpz_class(stage.outputScale.numerator),mpz_class(stage.outputScale.denominator))==mpq_class(a.scale(cipher)),"Exact dyadic output scale");
        check(a.coeffModulusValues(cipher)==stage.activePrimes&&a.chainIndex(cipher)==stage.chainIndex,"Runtime exact prime schedule");
        std::printf("branch=%zu factor=%zu %s chain=%zu scale=%.17g observed=%.8g certified=%.8g\n",stage.branch,stage.factor,stage.operation.c_str(),stage.chainIndex,a.scale(cipher),e,stage.semanticError.upperBound);
        ++observedNodes;
    });
    check(observedNodes==4*plan.metrics().depth+1,"All kernel/rescale/add stages observed");
    double observed=error(a.decodeComplex(a.decrypt(output.ciphertext)),reference);
    check(observed<=cert.outputError.upperBound,"Certified final StC error");
    check(std::isfinite(cert.gamma.upperBound)&&cert.gamma.upperBound>N,"Prepared gain includes nonzero encoding perturbation");
    const double expectedKappa[]={1,16,32,16}; mpq_class exactGain=2;
    for(std::size_t r=0;r<cert.metrics.depth;++r) {
        const auto& left=cert.factors[r].bounds; const auto& right=cert.factors[cert.metrics.depth+r].bounds;
        check(left.kappa.upperBound==expectedKappa[r]&&right.kappa.upperBound==expectedKappa[r],"Exact row norms of selected symbolic factors");
        exactGain*=mpq_class(expectedKappa[r])+mpq_class(std::max(left.delta.upperBound,right.delta.upperBound));
    }
    check(mpq_class(cert.gamma.upperBound)>=exactGain,"Outward gamma from actual prepared factors and pair addition");
    for(const auto& f:cert.factors) {
        check(!f.bounds.kappa.provenance.empty()&&!f.bounds.delta.provenance.empty()&&!f.bounds.localArithmeticError.provenance.empty(),"Every factor bound has provenance");
        if(f.branch<2) check(f.bounds.delta.upperBound>0&&f.internalHeadroomProofs.size()==5,"Prepared perturbations and all internal QP/headroom checks");
        check(f.bounds.kappa.kind!=BootstrapBoundKind::Unknown&&f.bounds.delta.kind!=BootstrapBoundKind::Unknown&&f.bounds.localArithmeticError.kind!=BootstrapBoundKind::Unknown,"Required factor bounds known");
        std::printf("factor b=%zu r=%zu kappa=%.17g delta=%.17g B_local=%.12g\n",f.branch,f.factor,f.bounds.kappa.upperBound,f.bounds.delta.upperBound,f.bounds.localArithmeticError.upperBound);
    }
    std::printf("FINAL gamma=%.17g certified=%.12g observed=%.12g rotations=%zu rescales=%zu levels=%zu\n",cert.gamma.upperBound,cert.outputError.upperBound,observed,cert.metrics.rotations,cert.metrics.rescales,cert.metrics.depth);
    // Actual size-three input path: two explicit input relinearizations, certified
    // by the same finite-support formulas, before the linear kernel starts.
    a.generateKeys(plan.requirements().rotations,true);
    auto fresh=a.encrypt(plain0); auto squared=a.multiply(fresh,fresh);
    auto squaredContract=contract;
    const mpq_class inputE(contract.inputError[0].upperBound),amp(amplitude);
    const mpq_class productError=2*amp*inputE+inputE*inputE;
    double productUpper=productError.get_d(); if(mpq_class(productUpper)<productError) productUpper=std::nextafter(productUpper,INFINITY);
    for(std::size_t b=0;b<2;++b) {
        squaredContract.inputMagnitude[b].upperBound=amplitude*amplitude;
        squaredContract.inputError[b]={productUpper,BootstrapBoundKind::Deterministic,"Exact ring square: 2*M*E+E^2 at exact power-of-two product scale",{}};
    }
    auto withRelin=plan.prepare(a,squared,squared,squaredContract);
    if(withRelin.certificate().result.status!=BootstrapCertificationStatus::Certified) throw std::runtime_error(withRelin.certificate().result.provenance);
    check(withRelin.certificate().keys.relinearization&&withRelin.certificate().metrics.relinearizations==2,"Executable input relinearization requirements");
    check(plan.preflight(noKeys,squared,squared,withRelin).status==BootstrapCertificationStatus::MissingEvaluationKeys,"Prepared relin preflight rejects missing key");
    ComplexVector squaredX=x; for(auto& value:squaredX) value*=value;
    auto r0=plan.applyPlainTrace(squaredX,0),r1=plan.applyPlainTrace(squaredX,1);
    auto combined=existingForwardOracle(squaredX,squaredX); std::size_t relinStages=0;
    plan.apply(a,squared,squared,withRelin,[&](const SlotToCoeffRuntimeStage& stage,const Cipher& c) {
        const auto& expected=stage.operation=="input relinearize"?squaredX:(stage.branch==2?combined:(stage.branch?r1:r0)[stage.factor]);
        check(error(a.decodeComplex(a.decrypt(c)),expected)<=stage.semanticError.upperBound,"Relinearized path per-stage error certificate");
        relinStages+=stage.operation=="input relinearize";
    });
    check(relinStages==2,"Both required relin operations observed");
}
#endif
}
int main() { try { plaintext();
#ifdef M2424_TEST_STC_BACKEND
backend();
#endif
std::puts("SlotToCoeff tests passed"); } catch(const std::exception& e) { std::fprintf(stderr,"%s\n",e.what()); return 1; } }
