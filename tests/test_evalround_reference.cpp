#include "m2424/experimental/evalmod_analysis/evalround_reference.hpp"
#include <cstdio>
#include <stdexcept>

using namespace m2424;
using namespace m2424::experimental;
namespace {
void check(bool value,const char* message) { if(!value) throw std::runtime_error(message); }
template<class F> void rejects(F f,const char* message) {
    try { f(); } catch(const std::invalid_argument&) { return; }
    throw std::runtime_error(message);
}
mpq_class normSquared(const EvalRoundExactComplex& x) { return x.real*x.real+x.imag*x.imag; }
void domains() {
    for(std::uint32_t K:{0,1,2,3,4,8,16}) for(double rho:{0.0,0.0625,0.125,0.25,0.49}) {
        EvalRoundProblem p{K,rho,1e-12,12};
        const mpq_class radius(rho);
        for(auto radix:{EvalRoundRadix::Binary,EvalRoundRadix::BalancedTernary}) {
            for(long I=-static_cast<long>(K);I<=K;++I) {
                for(int side:{-1,0,1}) {
                    const mpq_class z=I+side*radius;
                    const auto digits=extractEvalRoundReferenceDigits(p,radix,z);
                    check(evalRoundReferenceCenter(p,z)==I,"interval centers and both closed endpoints");
                    check(reconstructEvalRoundDigitsExact(digits,K,radix)==I,"exact zero-error reconstruction");
                }
                const mpq_class outside=mpq_class(I)+radius+mpq_class(1,1000000);
                if(rho<0.49) rejects([&]{extractEvalRoundReferenceDigits(p,radix,outside);},"outside interval");
            }
        }
    }
    EvalRoundProblem bad{1,0.5,1e-12,12};
    rejects([&]{evalRoundReferenceCenter(bad,0);},"rho=1/2 reference rejection");
}
void cleaners() {
    for(int bit:{0,1}) for(int k=-256;k<=256;++k) {
        const mpq_class e(k,256), x=bit+e;
        const mpq_class observed=abs(evalRoundBinaryCleanerExact(x)-bit);
        check(observed<=5*e*e,"exact f2 cleaning inequality");
        const mpq_class local(1,1024);
        check(abs(evalRoundBinaryCleanerExact(x)+local-bit)<=5*e*e+local,"binary cleaning plus local error");
    }
    check(evalRoundBinaryCleanerExact(0)==0 && evalRoundBinaryCleanerExact(1)==1,"exact binary fixed points");
    for(int trit:{-1,0,1}) {
        const auto root=evalRoundRootReference(trit,384);
        for(int k=1;k<=256;++k) for(int direction:{-1,1}) {
            const mpq_class a(k,256);
            const EvalRoundExactComplex error{direction*a*mpq_class(3,5),direction*a*mpq_class(4,5)};
            const auto value=evalRoundTernaryCleanerExact({root.real+error.real,root.imag+error.imag});
            const EvalRoundExactComplex difference{value.real-root.real,value.imag-root.imag};
            check(normSquared(difference)<=9*a*a*a*a,"f3 hard inequality at all three roots");
            const mpq_class local(1,1024), upper=3*a*a+local;
            check(normSquared({difference.real+local,difference.imag})<=upper*upper,"ternary cleaning plus local error");
        }
    }
    const auto fixed=evalRoundTernaryCleanerExact({1,0});
    check(fixed.real==1 && fixed.imag==0,"exact ternary root-one fixed point");
}
void extractionAndCleaning() {
    EvalRoundProblem p{1,0.01,1e-10,12};
    for(auto radix:{EvalRoundRadix::Binary,EvalRoundRadix::BalancedTernary}) {
        const auto method=radix==EvalRoundRadix::Binary?EvalRoundExtractionMethod::BinaryQuadraticK1:EvalRoundExtractionMethod::TernaryPhaseReferenceK1;
        const auto candidate=makeEvalRoundReferenceCandidate(p,radix,method);
        const auto plan=planEvalRoundCandidate(p,candidate);
        check(plan.status==EvalRoundPlanStatus::Certified,"reference plan has analytic bound");
        std::vector<std::size_t> counts; for(const auto& digit:plan.digits) counts.push_back(digit.cleaningIterations);
        for(long I=-1;I<=1;++I) for(int k=-32;k<=32;++k) {
            const mpq_class z=I+mpq_class(p.rho)*mpq_class(k,32);
            const auto trace=evaluateEvalRoundReference(p,candidate,z,counts,384);
            for(std::size_t j=0;j<trace.targetDigits.size();++j) {
                const auto target=radix==EvalRoundRadix::Binary?EvalRoundExactComplex{trace.targetDigits[j],0}:evalRoundRootReference(trace.targetDigits[j],384);
                for(std::size_t r=0;r<trace.digitStages[j].size();++r) {
                    const auto& value=trace.digitStages[j][r];
                    const EvalRoundExactComplex difference{value.real-target.real,value.imag-target.imag};
                    const mpq_class bound(plan.digits[j].errorAfterRounds[r]);
                    check(normSquared(difference)<=bound*bound,"observed reference stage <= analytical bound");
                }
            }
            check(abs(trace.reconstructed-I)<=mpq_class(plan.integerErrorUpper),"observed reconstruction <= analytical E_I");
        }
        std::printf("[evalround reference] radix=%u rounds=%zu certified E_I=%.9e all centers/endpoints/grid PASS\n",
            static_cast<unsigned>(radix),plan.totalCleaningIterations,plan.integerErrorUpper);
    }
    // Reference piecewise targets evaluate all digits/roots, including K>1.
    for(auto radix:{EvalRoundRadix::Binary,EvalRoundRadix::BalancedTernary}) {
        EvalRoundProblem q{8,0.25,1e-10,12};
        const auto c=makeEvalRoundReferenceCandidate(q,radix,EvalRoundExtractionMethod::PiecewiseReference);
        for(long I=-8;I<=8;++I) for(int side:{-1,1}) {
            const auto trace=evaluateEvalRoundReference(q,c,mpq_class(I)+mpq_class(side,4),std::vector<std::size_t>(c.digits.size()),384);
            if(radix==EvalRoundRadix::Binary) check(trace.reconstructed==I,"binary piecewise exact reference");
            else check(abs(trace.reconstructed-I)<mpq_class(1,ExactInteger(1)<<350),"ternary MPFR piecewise reference rounding diagnostic");
        }
    }
}
}
int main() {
    try { domains();cleaners();extractionAndCleaning();std::puts("[test_evalround_reference] PASS");return 0; }
    catch(const std::exception& e) {std::fprintf(stderr,"FAIL: %s\n",e.what());return 1;}
}
