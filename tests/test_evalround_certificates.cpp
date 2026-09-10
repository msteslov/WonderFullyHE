#include "m2424/experimental/evalmod_analysis/evalround_execution.hpp"
#include <cmath>
#include "evalround_test_support.hpp"
#include <cfenv>
#include <cstring>
#include <cstdio>
#include <stdexcept>
#include <type_traits>
using namespace m2424;
using namespace m2424::experimental;
namespace {
void check(bool x,const char* m) { if(!x) throw std::runtime_error(m); }
template<class F> void rejects(F f) { try { f(); } catch(const std::invalid_argument&) { return; } throw std::runtime_error("Executor accepted uncertified plan"); }
}
int main() {
    try {
        static_assert(!std::is_convertible_v<EvalRoundPlan,EvalRoundExecutionPlan>);
        CkksProfile profile{32768,std::vector<int>(10,60),std::ldexp(1.,59),1};
        auto a=SealAdapter::create(profile); a.generateKeys(std::vector<int>{0},true);
        auto context=a.encrypt(a.encode({0.}));
        auto input=evalround_test::encryptExactScalar(a,context,1./128,profile.scale);
        EvalRoundProblem p{1,1./128,1e-4,4};
        auto candidate=makeEvalRoundReferenceCandidate(p,EvalRoundRadix::Binary,EvalRoundExtractionMethod::BinaryQuadraticK1);
        auto ref=planEvalRoundCandidate(p,candidate);
        EvalRoundExecutionOptions options; options.inputSemanticError=evalround_test::inputBound(profile);
        auto compile=[&](const EvalRoundPlan& r,const EvalRoundExecutionOptions& o) { return EvalRoundExecutionCompiler::compile(a,input,r,o); };
        auto good=compile(ref,options);
        if(good.certification().status!=BootstrapCertificationStatus::Certified) throw std::runtime_error(good.certification().provenance);
        auto looseProblem=p; looseProblem.requiredIntegerError=1e-2;
        auto looseRef=planEvalRoundCandidate(looseProblem,makeEvalRoundReferenceCandidate(looseProblem,EvalRoundRadix::Binary,EvalRoundExtractionMethod::BinaryQuadraticK1));
        auto loose=compile(looseRef,options);
        check(loose.certification().status==BootstrapCertificationStatus::Certified && loose.mathematicalPlan().totalCleaningIterations<good.mathematicalPlan().totalCleaningIterations,"External requiredIntegerError changes minimum count");
        std::printf("required=%.3g rounds=%zu; required=%.3g rounds=%zu\n",looseProblem.requiredIntegerError,loose.mathematicalPlan().totalCleaningIterations,p.requiredIntegerError,good.mathematicalPlan().totalCleaningIterations);
        auto unknown=options; unknown.evaluationKeyNoiseCoefficientSupport={};
        const auto missing=compile(ref,unknown);
        check(missing.certification().status==BootstrapCertificationStatus::RequiredBoundUnavailable,"Unknown key bound must reject execution certification");
        rejects([&]{executeEvalRound(a,input,missing);});
        unknown=options; unknown.inputSemanticError={};
        check(compile(ref,unknown).certification().status==BootstrapCertificationStatus::RequiredBoundUnavailable,"Unknown upstream bound");
        unknown=options; unknown.evaluationKeyNoiseCoefficientSupport.upperBound=0;
        check(compile(ref,unknown).certification().status==BootstrapCertificationStatus::RequiredBoundUnavailable,"Synthetic zero noise support rejected");
        for(auto method:{EvalRoundExtractionMethod::PiecewiseReference,EvalRoundExtractionMethod::TernaryPhaseReferenceK1}) {
            auto r=planEvalRoundCandidate(p,makeEvalRoundReferenceCandidate(p,method==EvalRoundExtractionMethod::PiecewiseReference?EvalRoundRadix::Binary:EvalRoundRadix::BalancedTernary,method));
            check(r.status==EvalRoundPlanStatus::Certified,"Reference-only target is mathematically certified");
            auto no=compile(r,options);
            check(no.certification().status==BootstrapCertificationStatus::ExtractionNotCertified,"Reference-only extraction rejected");
            rejects([&]{executeEvalRound(a,input,no);});
        }
        auto shifted=evalround_test::encryptExactScalar(a,context,1./128,std::nextafter(profile.scale,INFINITY));
        auto shiftedPlan=EvalRoundExecutionCompiler::compile(a,shifted,ref,options);
        check(shiftedPlan.certification().status==BootstrapCertificationStatus::Certified,"Non-power-of-two exact dyadic input scale compiles");
        auto otherProfile=profile; otherProfile.coeffModulusBits.back()=59;
        auto otherContext=SealAdapter::create(otherProfile);
        check(preflightEvalRound(otherContext,input,good).status==BootstrapCertificationStatus::InvalidInput,"Context fingerprint enforced before execution");
        auto rewritten=a.normalizeScale(input,std::nextafter(a.scale(input),INFINITY));
        check(preflightEvalRound(a,rewritten,good).status==BootstrapCertificationStatus::InputScaleMismatch,"Exact scale bits enforced");
        auto down=a.rescaleToNext(input);
        check(preflightEvalRound(a,down,good).status!=BootstrapCertificationStatus::Certified,"Wrong starting chain rejected");
        auto shallow=a.rescaleToNext(a.rescaleToNext(a.rescaleToNext(a.rescaleToNext(input))));
        // Restore input scale by *encoding/encrypting* a fresh plaintext at that level.
        shallow=a.encrypt(a.encodeScalarAtScaleFor(.01,profile.scale,shallow));
        auto noLevels=EvalRoundExecutionCompiler::compile(a,shallow,ref,options);
        check(noLevels.certification().status!=BootstrapCertificationStatus::Certified,"Insufficient levels rejected before execution");
        bool ran=false;
        rejects([&]{executeEvalRound(a,shallow,noLevels,[&](std::size_t,const Cipher&){ran=true;});});
        check(!ran,"No node executed after rejected preflight");
        // Minimum count is recomputed with real local errors, not copied from ref.
        const auto& selected=good.mathematicalPlan();
        std::size_t minimum=100;
        for(std::size_t x=0;x<=selected.digits[0].cleaningIterations;++x) for(std::size_t y=0;y<=selected.digits[1].cleaningIterations;++y) {
            std::vector<double> errors,locals;
            for(std::size_t d=0;d<2;++d) {
                double e=selected.digits[d].extractionError.upperBound;
                for(std::size_t r=0;r<(d?y:x);++r) e=evalRoundCleaningErrorUpper(EvalRoundRadix::Binary,e,selected.digits[d].cleaningLocalErrors[r].upperBound);
                errors.push_back(e); locals.push_back(selected.digits[d].reconstructionLocalError.upperBound);
            }
            if(evalRoundReconstructionErrorUpper(EvalRoundRadix::Binary,errors,locals)<=p.requiredIntegerError) minimum=std::min(minimum,x+y);
        }
        check(selected.totalCleaningIterations==minimum,"Minimum planner cleaning count with backend errors");
        bool representationError=false;
        std::size_t rescaleNodes=0;
        for(const auto& n:good.nodes()) {
            if(n.scaleRepresentationError.upperBound>0) representationError=true;
            if(n.operation==EvalRoundOperation::Rescale) ++rescaleNodes;
            double runtimeScale; auto scaleBits=n.outputScale.binary64Bits;
            std::memcpy(&runtimeScale,&scaleBits,sizeof runtimeScale);
            check(mpq_class(mpz_class(n.outputScale.numerator),mpz_class(n.outputScale.denominator))==mpq_class(runtimeScale),"Exact runtime scale is binary64 dyadic, not nominal power of two");
            const mpq_class exactRuntime(runtimeScale);
            const mpz_class N(std::to_string(profile.polyModulusDegree));
            const mpq_class roundBound=mpq_class(N*(1+N))/2/exactRuntime;
            if(n.operation==EvalRoundOperation::Rescale)
                check(mpq_class(n.localArithmeticError.upperBound)>=roundBound,"Rescale certificate contains hard two-component rounding bound");
            if(n.requiredKey!=EvalRoundEvaluationKey::None) {
                mpz_class sum=0; for(auto prime:n.activePrimes) sum+=mpz_class(std::to_string(prime-1));
                const mpq_class keyBound=mpq_class(N*N*21*sum)/mpz_class(std::to_string(a.specialKeyModulusValue()))/exactRuntime+roundBound;
                check(mpq_class(n.localArithmeticError.upperBound)>=keyBound,"Key-switch certificate contains support noise and ModDown rounding");
            }
            if(n.operation==EvalRoundOperation::Rescale) {
                const auto& source=good.nodes()[n.inputs[0]];
                const mpq_class before(mpz_class(source.outputScale.numerator),mpz_class(source.outputScale.denominator));
                const mpq_class after(mpz_class(n.arithmeticScale.numerator),mpz_class(n.arithmeticScale.denominator));
                check(after==before/mpz_class(std::to_string(source.activePrimes.back())),"Exact scale division uses actual integer prime");
            }
            for(const auto* b:{&n.idealMagnitude,&n.semanticError,&n.localArithmeticError,&n.scaleRepresentationError})
                check(b->kind!=BootstrapBoundKind::Unknown && !b->provenance.empty(),"Every required node bound known with provenance");
            if(n.operation==EvalRoundOperation::Conjugate || n.operation==EvalRoundOperation::Relinearize || n.operation==EvalRoundOperation::Rescale)
                check(n.localArithmeticError.upperBound>0,"Real key switching/rescale local errors cannot be zero");
        }
        check(representationError,"Scale representation error is not silently dropped");
        check(rescaleNodes==2+2*selected.totalCleaningIterations,"Baseline two rescale nodes per cleaner, no extra scale-normalization rescale");
        check(selected.log2FailureProbabilityUpper==-INFINITY,"Deterministic execution certificate has zero failure probability");
        std::fesetround(FE_UPWARD);
        auto roundingGate=preflightEvalRound(a,input,good);
        std::fesetround(FE_TONEAREST);
        check(roundingGate.status==BootstrapCertificationStatus::ScaleScheduleInfeasible,"Rounding mode cannot change after compilation");
        // Regenerate only the available key subsets; the immutable plan must preflight
        // against current key availability before touching ciphertext arithmetic.
        a.generateKeys(std::vector<int>{0},false);
        check(preflightEvalRound(a,input,good).status==BootstrapCertificationStatus::MissingEvaluationKeys,"Missing relin preflight");
        a.generateKeys(std::vector<int>{},true);
        check(preflightEvalRound(a,input,good).status==BootstrapCertificationStatus::MissingEvaluationKeys,"Missing conjugation preflight");
        std::puts("EvalRound execution certificate rejection and minimum-count tests passed");
    } catch(const std::exception& e) { std::fprintf(stderr,"%s\n",e.what()); return 1; }
}
