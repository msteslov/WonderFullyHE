#include "m2424/evalround.hpp"
#include <algorithm>
#include <cmath>
#include <cfenv>
#include <cstdio>
#include <stdexcept>
#include <numeric>

using namespace m2424;
namespace {
void check(bool condition,const char* message) { if (!condition) throw std::runtime_error(message); }
template<class F> void rejects(F f,const char* message) {
    try { f(); } catch (const std::invalid_argument&) { return; }
    throw std::runtime_error(message);
}
void integerTests() {
    for (std::uint32_t K : {0,1,2,3,4,8,16,32}) for (auto radix : {EvalRoundRadix::Binary,EvalRoundRadix::BalancedTernary}) {
        const auto count=evalRoundDigitCount(K,radix);
        std::uint64_t capacity=1; for(std::size_t i=0;i<count;++i) capacity*=static_cast<unsigned>(radix);
        check(capacity>=2ULL*K+1 && (!count || capacity/static_cast<unsigned>(radix)<2ULL*K+1),"minimum digit count");
        for(std::int64_t I=-static_cast<std::int64_t>(K);I<=K;++I) {
            const auto digits=evalRoundIntegerDigits(I,K,radix);
            check(reconstructEvalRoundInteger(digits,K,radix)==I,"exact signed reconstruction");
            if(radix==EvalRoundRadix::Binary) {
                std::uint64_t J=0,weight=1;
                for(int d:digits) { check(d==0||d==1,"binary alphabet"); J+=weight*d; weight*=2; }
                check(J==static_cast<std::uint64_t>(I+K),"binary J=I+K, not twos complement");
            }
        }
    }
    for(auto radix:{EvalRoundRadix::Binary,EvalRoundRadix::BalancedTernary}) {
        const auto K=std::numeric_limits<std::uint32_t>::max();
        for(std::int64_t I:{-static_cast<std::int64_t>(K),std::int64_t(0),static_cast<std::int64_t>(K)})
            check(reconstructEvalRoundInteger(evalRoundIntegerDigits(I,K,radix),K,radix)==I,"large exact K");
    }
    rejects([]{evalRoundDigitCount(1,static_cast<EvalRoundRadix>(4));},"invalid radix");
    rejects([]{evalRoundIntegerDigits(-2,1,EvalRoundRadix::Binary);},"integer outside domain");
    rejects([]{reconstructEvalRoundInteger({-1,1},1,EvalRoundRadix::Binary);},"invalid binary digit");
}
void minimumCounts(const EvalRoundProblem& p,const EvalRoundCandidate& c,const EvalRoundPlan& plan) {
    check(plan.status==EvalRoundPlanStatus::Certified,"representative candidate certified");
    std::size_t minimum=1000;
    // Independent exhaustive enumeration for one trit or two shifted binary digits.
    for(std::size_t a=0;a<=p.maxCleaningRoundsPerDigit;++a) for(std::size_t b=0;b<=(c.digits.size()==2?p.maxCleaningRoundsPerDigit:0);++b) {
        std::vector<double> errors,local;
        bool valid=true;
        for(std::size_t j=0;j<c.digits.size();++j) {
            double e=c.digits[j].extractionError.upperBound;
            for(std::size_t r=0;r<(j==0?a:b);++r) {
                if(e>1) { valid=false; break; }
                e=evalRoundCleaningErrorUpper(c.radix,e,c.digits[j].cleaningLocalErrors[r].upperBound);
            }
            errors.push_back(e); local.push_back(c.digits[j].reconstructionLocalError.upperBound);
        }
        if(valid && evalRoundReconstructionErrorUpper(c.radix,errors,local)<=p.requiredIntegerError) minimum=std::min(minimum,a+b);
    }
    check(plan.totalCleaningIterations==minimum,"globally minimum total cleaning count");
    std::size_t count=0; std::vector<double> errors,locals;
    for(std::size_t j=0;j<plan.digits.size();++j) {
        count+=plan.digits[j].cleaningIterations;
        errors.push_back(plan.digits[j].errorAfterRounds.back()); locals.push_back(c.digits[j].reconstructionLocalError.upperBound);
    }
    check(count==minimum && evalRoundReconstructionErrorUpper(c.radix,errors,locals)==plan.integerErrorUpper,"certificate reconstruction matches counts");
}
void plannerTests() {
    EvalRoundProblem p{1,0.01,1e-8,12};
    auto binary=makeEvalRoundReferenceCandidate(p,EvalRoundRadix::Binary,EvalRoundExtractionMethod::BinaryQuadraticK1,{1,0.1,0.2});
    auto ternary=makeEvalRoundReferenceCandidate(p,EvalRoundRadix::BalancedTernary,EvalRoundExtractionMethod::TernaryPhaseReferenceK1,{10,0.1,0.2});
    auto selection=planEvalRound(p,{binary,ternary});
    check(selection.selected && *selection.selected==0,"K=1 configured cost selects binary over one trit");
    for(std::size_t i=0;i<2;++i) {
        const auto& c=i?ternary:binary; const auto& plan=selection.candidates[i]; minimumCounts(p,c,plan);
        std::printf("[evalround candidate] radix=%u rounds=%zu E_I=%.9e cost=%.3f\n",static_cast<unsigned>(c.radix),plan.totalCleaningIterations,plan.integerErrorUpper,plan.configuredCost);
    }
    binary.cost.extraction=10; ternary.cost.extraction=1;
    selection=planEvalRound(p,{binary,ternary});
    check(selection.selected && *selection.selected==1,"configured cost selects one trit");
    auto tighter=p; tighter.requiredIntegerError=1e-16;
    minimumCounts(tighter,binary,planEvalRoundCandidate(tighter,binary));
    auto withLocal=ternary;
    for(auto& bound:withLocal.digits[0].cleaningLocalErrors) bound.upperBound=1e-12;
    withLocal.digits[0].reconstructionLocalError.upperBound=1e-12;
    minimumCounts(p,withLocal,planEvalRoundCandidate(p,withLocal));
    auto noRounds=p; noRounds.maxCleaningRoundsPerDigit=0;
    check(planEvalRoundCandidate(noRounds,binary).rejection==EvalRoundRejectionReason::ErrorBudgetExceeded,"bounded search does not invent cleaning rounds");
    auto tooMany=p; tooMany.maxCleaningRoundsPerDigit=65;
    check(planEvalRoundCandidate(tooMany,binary).rejection==EvalRoundRejectionReason::InvalidProblem,"bounded search limit validation");
    auto unevenBudget=p; unevenBudget.requiredIntegerError=1.5e-5;
    const auto unevenPlan=planEvalRoundCandidate(unevenBudget,binary);
    minimumCounts(unevenBudget,binary,unevenPlan);
    check(unevenPlan.digits[0].cleaningIterations==3 && unevenPlan.digits[1].cleaningIterations==2,
        "real quadratic bounds need nonuniform per-digit counts");
    auto bad=binary; bad.extraction.verified=false;
    check(planEvalRoundCandidate(p,bad).status==EvalRoundPlanStatus::Diagnostic,"uncertified extraction diagnostic only");
    selection=planEvalRound(p,{bad,ternary});
    check(selection.selected && *selection.selected==1,"cheap diagnostic cannot win");
    bad=binary; bad.digits[0].extractionError={}; bad.digits[0].extractionError.upperBound=0;
    check(planEvalRoundCandidate(p,bad).rejection==EvalRoundRejectionReason::BoundUnavailable,"unknown extraction is not zero");
    bad=binary; bad.digits[0].cleaningLocalErrors[0]={};
    check(planEvalRoundCandidate(p,bad).rejection==EvalRoundRejectionReason::BoundUnavailable,"missing reachable cleaning bound");
    bad=binary; bad.digits[0].reconstructionLocalError={};
    check(planEvalRoundCandidate(p,bad).rejection==EvalRoundRejectionReason::BoundUnavailable,"missing reconstruction bound");
    bad=binary; bad.extraction.certifiedRho=0.02;
    check(planEvalRoundCandidate(p,bad).rejection==EvalRoundRejectionReason::ExtractionNotCertified,"certificate domain binding");
    bad=binary; bad.digits[0].extractionError.upperBound=1.1;
    check(planEvalRoundCandidate(p,bad).rejection==EvalRoundRejectionReason::CleaningDomainViolation,"cleaner a<=1 domain");
    bad=binary; bad.digits[0].reconstructionLocalError.upperBound=0.1;
    check(planEvalRoundCandidate(p,bad).rejection==EvalRoundRejectionReason::ErrorBudgetExceeded,"irreducible reconstruction floor");
    bad=ternary; bad.extraction.method=EvalRoundExtractionMethod::BinaryQuadraticK1;
    check(planEvalRoundCandidate(p,bad).rejection==EvalRoundRejectionReason::InvalidCandidate,"method/radix mismatch");
    bad=binary; bad.extraction.certifiedK=2;
    check(planEvalRoundCandidate(p,bad).rejection==EvalRoundRejectionReason::ExtractionNotCertified,"certified K mismatch");
    bad=binary; bad.cost.extraction=-1;
    check(planEvalRoundCandidate(p,bad).rejection==EvalRoundRejectionReason::InvalidCandidate,"invalid cost");
    bad=binary; bad.digits.pop_back();
    check(planEvalRoundCandidate(p,bad).rejection==EvalRoundRejectionReason::InvalidCandidate,"wrong digit count");
    bad=binary; bad.digits[0].extractionError.provenance.clear();
    check(planEvalRoundCandidate(p,bad).rejection==EvalRoundRejectionReason::BoundUnavailable,"missing bound provenance");
    bad=binary; bad.radix=static_cast<EvalRoundRadix>(4);
    check(planEvalRoundCandidate(p,bad).rejection==EvalRoundRejectionReason::InvalidCandidate,"unknown planner radix");
    bad=binary; bad.extraction.method=static_cast<EvalRoundExtractionMethod>(99);
    check(planEvalRoundCandidate(p,bad).rejection==EvalRoundRejectionReason::InvalidCandidate,"unknown extraction method");
    const int originalRounding=std::fegetround();
    check(std::fesetround(FE_UPWARD)==0,"set unsupported rounding");
    const auto unsupported=planEvalRoundCandidate(p,binary);
    check(std::fesetround(originalRounding)==0,"restore rounding mode");
    check(unsupported.rejection==EvalRoundRejectionReason::InvalidProblem,"unsupported arithmetic mode fails closed");
    auto invalid=p; invalid.requiredIntegerError=0;
    check(planEvalRoundCandidate(invalid,binary).rejection==EvalRoundRejectionReason::InvalidProblem,"budget must be supplied");
    for(double rho:{0.5,0.6,-0.01,std::numeric_limits<double>::quiet_NaN()}) {
        invalid=p; invalid.rho=rho;
        check(planEvalRoundCandidate(invalid,binary).rejection==EvalRoundRejectionReason::DomainViolation,"invalid rho");
    }
    for(double rho:{0.0,0.125,0.2}) {
        auto domain=p; domain.rho=rho;
        auto piecewise=makeEvalRoundReferenceCandidate(domain,EvalRoundRadix::Binary,EvalRoundExtractionMethod::PiecewiseReference);
        check(planEvalRoundCandidate(domain,piecewise).status==EvalRoundPlanStatus::Certified,"DigitExtract restriction is not universal");
        piecewise.extraction.method=EvalRoundExtractionMethod::DigitExtract;
        piecewise.extraction.digitExtractEpsilon=0.25;
        check(planEvalRoundCandidate(domain,piecewise).rejection==EvalRoundRejectionReason::DigitExtractDomainViolation,"DigitExtract specific rho gate");
    }
    bad=binary; bad.extraction.method=EvalRoundExtractionMethod::DigitExtract; bad.extraction.verified=false;
    for(double epsilon:{0.02,0.01,0.26,std::numeric_limits<double>::quiet_NaN()}) {
        bad.extraction.digitExtractEpsilon=epsilon;
        check(planEvalRoundCandidate(p,bad).rejection==EvalRoundRejectionReason::DigitExtractDomainViolation,"DigitExtract epsilon gate");
    }
    bad.extraction.digitExtractEpsilon=0.1;
    check(planEvalRoundCandidate(p,bad).rejection==EvalRoundRejectionReason::ExtractionNotCertified,"valid DigitExtract domain still requires proof");
    bad=binary; bad.failureEvents={{"extraction",-129,"primary tail"}};
    bad.digits[0].extractionError.kind=BootstrapBoundKind::Probabilistic;
    bad.digits[0].extractionError.failureEventIds={"extraction"};
    check(planEvalRoundCandidate(p,bad).status==EvalRoundPlanStatus::Certified,"declared probabilistic certificate");
    bad.failureEvents[0].log2FailureProbability=-127;
    check(planEvalRoundCandidate(p,bad).rejection==EvalRoundRejectionReason::FailureProbabilityExceeded,"failure gate");
    bad.failureEvents.clear();
    check(planEvalRoundCandidate(p,bad).rejection==EvalRoundRejectionReason::BoundUnavailable,"unresolved event");
    check(!planEvalRound(p,{}).selected,"empty search");
    check(planEvalRound(p,{binary,binary}).rejection==EvalRoundRejectionReason::InvalidCandidate,"duplicate candidate IDs");
    for(auto K:{0U,2U,8U}) {
        EvalRoundProblem q{K,0.2,1e-20,12};
        auto b=makeEvalRoundReferenceCandidate(q,EvalRoundRadix::Binary,EvalRoundExtractionMethod::PiecewiseReference,{2,1,1});
        auto t=makeEvalRoundReferenceCandidate(q,EvalRoundRadix::BalancedTernary,EvalRoundExtractionMethod::PiecewiseReference,{1,1,1});
        const auto r=planEvalRound(q,{b,t});
        check(r.selected && r.candidates[*r.selected].integerErrorUpper==0 && r.candidates[*r.selected].totalCleaningIterations==0,"exact target needs no cleaning");
        std::printf("[evalround selection] K=%u rho=0.2 radix=3 rounds=0 E_I=0 reference-only\n",K);
    }
    rejects([]{evalRoundCleaningErrorUpper(EvalRoundRadix::Binary,1.001,0);},"binary cleaner domain");
    rejects([]{evalRoundCleaningErrorUpper(EvalRoundRadix::BalancedTernary,0.1,-1);},"negative local error");
    check(evalRoundCleaningErrorUpper(EvalRoundRadix::Binary,0.125,0)>=5*0.125*0.125,"hard constant 5");
    check(evalRoundCleaningErrorUpper(EvalRoundRadix::BalancedTernary,0.125,0)>=3*0.125*0.125,"hard constant 3");
    check(evalRoundReconstructionErrorUpper(EvalRoundRadix::Binary,{0.125,0.0625},{0,0})==0.25,"binary weighted reconstruction");
    const double trit=evalRoundReconstructionErrorUpper(EvalRoundRadix::BalancedTernary,{0.125},{0});
    check(trit>=0x1.279a74590331dp-3 && trit<=0x1.279a745903320p-3,"ternary 2/sqrt3 outward reconstruction gain");
}
}
int main() {
    try { integerTests(); plannerTests(); std::puts("[test_evalround_planner] PASS"); return 0; }
    catch(const std::exception& e) { std::fprintf(stderr,"FAIL: %s\n",e.what()); return 1; }
}
