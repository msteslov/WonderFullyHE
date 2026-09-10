#include "m2424/experimental/evalmod_analysis/evalround_execution.hpp"
#include "m2424/experimental/evalmod_analysis/evalround_reference.hpp"
#include <mpfr.h>
#include "evalround_test_support.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
using namespace m2424;
using namespace m2424::experimental;
namespace {
void check(bool x,const char* m) { if(!x) throw std::runtime_error(m); }
struct Complex { mpq_class re,im; };
Complex add(Complex a,Complex b) { return {a.re+b.re,a.im+b.im}; }
Complex mul(Complex a,Complex b) { return {a.re*b.re-a.im*b.im,a.re*b.im+a.im*b.re}; }
double asDouble(const mpq_class& x) {
    mpfr_t v; mpfr_init2(v,256); mpfr_set_q(v,x.get_mpq_t(),MPFR_RNDN);
    double result=mpfr_get_d(v,MPFR_RNDN); mpfr_clear(v); return result;
}
}
int main() {
    try {
        CkksProfile profile{32768,std::vector<int>(10,60),std::ldexp(1.,59),9};
        auto a=SealAdapter::create(profile); a.generateKeys(std::vector<int>{0},true);
        EvalRoundProblem problem{1,1./128,1e-4,4};
        const auto candidate=makeEvalRoundReferenceCandidate(problem,EvalRoundRadix::Binary,EvalRoundExtractionMethod::BinaryQuadraticK1);
        const auto reference=planEvalRoundCandidate(problem,candidate);
        auto context=a.encrypt(a.encode({0.}));
        auto input=evalround_test::encryptExactScalar(a,context,0.,profile.scale);
        EvalRoundExecutionOptions options;
        options.inputSemanticError=evalround_test::inputBound(profile);
        auto plan=EvalRoundExecutionCompiler::compile(a,input,reference,options);
        if(plan.certification().status!=BootstrapCertificationStatus::Certified) throw std::runtime_error(plan.certification().provenance);
        check(plan.integerErrorUpper()<=problem.requiredIntegerError,"Required integer budget");
        const auto& math=plan.mathematicalPlan();
        std::printf("N=%zu special_prime=%llu input_primes=",profile.polyModulusDegree,static_cast<unsigned long long>(a.specialKeyModulusValue()));
        for(auto prime:a.coeffModulusValues(input)) std::printf("%llu,",static_cast<unsigned long long>(prime));
        std::puts("");
        std::printf("plan=%s rounds=%zu,%zu required=%.12g certified=%.12g nodes=%zu\n",math.candidateId.c_str(),math.digits[0].cleaningIterations,math.digits[1].cleaningIterations,problem.requiredIntegerError,plan.integerErrorUpper(),plan.nodes().size());
        double worstIntegerError=0;
        for(int center=-1;center<=1;++center) for(double offset:{-problem.rho,0.,problem.rho}) {
            const std::vector<double> z{center+offset};
            input=evalround_test::encryptExactScalar(a,context,z[0],profile.scale);
            std::vector<std::vector<Complex>> values(plan.nodes().size());
            std::size_t seen=0;
            const auto output=executeEvalRound(a,input,plan,[&](std::size_t i,const Cipher& c) {
                const auto& n=plan.nodes()[i]; values[i].resize(z.size());
                const auto decoded=a.decodeComplex(a.decrypt(c));
                double observed=0;
                for(std::size_t j=0;j<z.size();++j) {
                    auto& v=values[i][j];
                    if(n.operation==EvalRoundOperation::Input) v={mpq_class(z[j]),0};
                    else {
                        v=values[n.inputs[0]][j];
                        switch(n.operation) {
                        case EvalRoundOperation::Multiply: v=mul(v,values[n.inputs[1]][j]); break;
                        case EvalRoundOperation::Add: v=add(v,values[n.inputs[1]][j]); break;
                        case EvalRoundOperation::Subtract: { auto b=values[n.inputs[1]][j]; v=add(v,{-b.re,-b.im}); break; }
                        case EvalRoundOperation::MultiplyPlain: {
                            mpq_class k(mpz_class(n.constantNumerator),mpz_class(n.constantDenominator)); v=mul(v,{k,0}); break;
                        }
                        case EvalRoundOperation::AddPlain: v.re+=mpq_class(mpz_class(n.constantNumerator),mpz_class(n.constantDenominator)); break;
                        case EvalRoundOperation::Conjugate: v.im=-v.im; break;
                        default: break; // exact semantic identity: relin/rescale/modswitch
                        }
                    }
                    // Check every CKKS slot against the exact rational polynomial oracle.
                    for(const auto& slot:decoded) {
                    mpq_class re=mpq_class(slot.real())-v.re,im=mpq_class(slot.imag())-v.im;
                    mpq_class squared=re*re+im*im;
                    mpq_class bound(n.semanticError.upperBound);
                    check(squared<=bound*bound,"Per-node observed error exceeds certified semantic error");
                    observed=std::max(observed,std::sqrt(asDouble(squared)));
                    }
                }
                check(a.coeffModulusValues(c)==n.activePrimes,"Exact active primes match runtime");
                check(mpz_class(n.centeredHeadroomNumerator)>0,"Centered headroom");
                ++seen;
                if(center==-1 && offset==-problem.rho) std::printf("node=%zu stage=%s level=%zu scale=%.17g observed=%.4e bound=%.4e local=%.4e\n",i,n.stage.c_str(),n.chainIndex,a.scale(c),observed,n.semanticError.upperBound,n.localArithmeticError.upperBound);
            });
            check(seen==plan.nodes().size(),"All extraction, cleaning and reconstruction nodes observed");
            auto decoded=a.decodeComplex(a.decrypt(output)); double integerError=0;
            for(const auto& slot:decoded) {
                const mpq_class re=mpq_class(slot.real())-center,im(slot.imag()),limit(plan.integerErrorUpper());
                check(re*re+im*im<=limit*limit,"Exact final integer-error comparison");
                integerError=std::max(integerError,std::abs(slot-double(center)));
            }
            worstIntegerError=std::max(worstIntegerError,integerError);
            check(integerError<=plan.integerErrorUpper(),"Observed integer error exceeds reconstruction certificate");
            // Independent PR-2 high precision/exact polynomial trace agrees with DAG.
            const std::vector<std::size_t> counts{math.digits[0].cleaningIterations,math.digits[1].cleaningIterations};
            for(std::size_t j=0;j<z.size();++j) {
                check(evalRoundReferenceCenter(problem,mpq_class(z[j]))==center,"Exact whole-domain center and endpoints");
                mpq_class x(z[j]),square=x*x;
                mpq_class digits[2]={1-square,(square+x)/2};
                for(std::size_t d=0;d<2;++d) for(std::size_t r=0;r<counts[d];++r) digits[d]=evalRoundBinaryCleanerExact(digits[d]);
                check(values[plan.outputNode()][j].re==digits[0]+2*digits[1]-1,"Exact PR-2 cleaner oracle matches compiled DAG");
            }
            std::printf("domain z=%.12g observed_integer_error=%.12g\n",z[0],integerError);
        }
        auto pair=executeEvalRoundPair(a,input,input,plan);
        check(pair.integerErrorUpper==plan.integerErrorUpper(),"Both halves carry same conditional bound");
        for(const auto* c:{&pair.first,&pair.second}) check(a.info(*c).chainIndex==plan.nodes().back().chainIndex,"Pair schedule");
        std::printf("FINAL required=%.12g certified=%.12g observed=%.12g chain=%zu scale=%.17g\n",problem.requiredIntegerError,plan.integerErrorUpper(),worstIntegerError,a.chainIndex(pair.first),a.scale(pair.first));
    } catch(const std::exception& e) { std::fprintf(stderr,"%s\n",e.what()); return 1; }
}
