#include "m2424/experimental/evalmod_analysis/evalround_execution.hpp"
#include "m2424/experimental/evalmod_analysis/evalround_reference.hpp"
#include "m2424/experimental/evalmod_analysis/exact_decimal.hpp"
#include "bootstrap_fixture.hpp"
#include "../src/math/evalround_interval_internal.hpp"
#include "../src/planning/evalround_polynomial_internal.hpp"
#include "../src/core/certified_arithmetic_internal.hpp"
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <cmath>
#include <cstring>
using namespace m2424;
using namespace m2424::experimental;
using S=BootstrapCertificationStatus;
void check(bool b,const char* why){if(!b)throw std::runtime_error(why);}
BootstrapBound B(double x){return {x,BootstrapBoundKind::Deterministic,"Test analytical bound",{}};}
int main(){try {
    auto exactValue=[](const std::vector<mpq_class>& coefficients,const mpq_class& x){
        mpq_class value=0;for(auto it=coefficients.rbegin();it!=coefficients.rend();++it)value=value*x+*it;
        value.canonicalize();return value;
    };
    auto checkShift=[&](const EvalModPolynomial& polynomial,std::int64_t center,const std::vector<mpq_class>& ys){
        std::vector<mpq_class> original;for(const auto& text:polynomial.decimalCoefficients)original.push_back(parseExactDecimal(text));
        const auto shifted=detail::shiftPolynomialToIntegerCenterExact(polynomial,center);
        for(const auto& y:ys)check(exactValue(original,mpq_class(static_cast<long>(center))+y)==exactValue(shifted,y),
            "Exact centered coefficients preserve the same polynomial");
    };
    const std::vector<mpq_class> rationalOffsets{mpq_class(-3,8),mpq_class(0),mpq_class(7,16)};
    for(const auto& polynomial:std::vector<EvalModPolynomial>{
            {PolynomialBasis::Monomial,{"0.125"}},
            {PolynomialBasis::Monomial,{"-0.5","1.25"}},
            {PolynomialBasis::Monomial,{"0.25","-0.75","1.5"}},
            {PolynomialBasis::Monomial,{"-1","0.5","0.125","-0.25"}}}) {
        checkShift(polynomial,-3,rationalOffsets);checkShift(polynomial,2,rationalOffsets);
    }
    EvalModPolynomial sparse{PolynomialBasis::Monomial,std::vector<std::string>(129,"0")};
    sparse.decimalCoefficients[0]="0.25";sparse.decimalCoefficients[3]="-0.5";
    sparse.decimalCoefficients[64]="0.125";sparse.decimalCoefficients[128]="-0.03125";
    checkShift(sparse,-2,rationalOffsets);checkShift(sparse,3,rationalOffsets);

    auto a=test::BootstrapFixture::create({16,std::vector<int>(10,50),std::ldexp(1.,49),8});
    a.generateKeys(std::vector<int>{0},true);
    auto input=a.encrypt(a.encode({0.}));
    EvalRoundExecutionOptions options;options.inputSemanticError=B(1e-8);
    {
        const EvalModPolynomial q0{PolynomialBasis::Chebyshev,{"-2047","0","-2048"}};
        const EvalModPolynomial p0{PolynomialBasis::Monomial,{"1","0","-1"}};
        check(arithmetic::exactPolynomialEqual(
                  convertScaledChebyshevToMonomial(q0,"64"),p0),
              "scaled Chebyshev conversion is exactly equal as a rational polynomial");
        auto chebAdapter=test::BootstrapFixture::create(
            {16,std::vector<int>(48,50),std::ldexp(1.,49),8});
        chebAdapter.generateKeys(std::vector<int>{0},true);
        const auto chebInput=chebAdapter.encrypt(chebAdapter.encode({0.}));
        arithmetic::Builder chebBuilder(chebAdapter,chebInput,
                                        options.evaluationKeyNoiseCoefficientSupport);
        const mpq_class R(65,64),sample(3,4);
        const auto tNode=chebBuilder.input(chebAdapter.scale(chebInput),R,0);
        arithmetic::ScaledChebyshevCompiler chebyshev(chebBuilder,tNode,R);
        const auto t2=chebyshev.compileBasis(2);
        const auto t3=chebyshev.compileBasis(3);
        const auto t8=chebyshev.compileBasis(8);
        const auto t255=chebyshev.compileBasis(255);
        const auto t256=chebyshev.compileBasis(256);
        auto idealAt=[&](std::size_t output,const mpq_class& value) {
            std::vector<mpq_class> ideal(chebBuilder.nodes.size());
            for(std::size_t node=0;node<=output;++node) {
                const auto& operation=chebBuilder.nodes[node];
                if(operation.operation==EvalRoundOperation::Input)ideal[node]=value;
                else {
                    ideal[node]=ideal[operation.inputs[0]];
                    switch(operation.operation) {
                    case EvalRoundOperation::Multiply:
                        ideal[node]*=ideal[operation.inputs[1]];break;
                    case EvalRoundOperation::Add:
                        ideal[node]+=ideal[operation.inputs[1]];break;
                    case EvalRoundOperation::Subtract:
                        ideal[node]-=ideal[operation.inputs[1]];break;
                    case EvalRoundOperation::MultiplyPlain:
                        ideal[node]*=mpq_class(mpz_class(operation.constantNumerator),
                                              mpz_class(operation.constantDenominator));break;
                    case EvalRoundOperation::AddPlain:
                        ideal[node]+=mpq_class(mpz_class(operation.constantNumerator),
                                              mpz_class(operation.constantDenominator));break;
                    default:break;
                    }
                }
            }
            return ideal[output];
        };
        auto referenceT=[](std::size_t degree,const mpq_class& value) {
            if(degree==0)return mpq_class(1);
            mpq_class previous=1,current=value;
            for(std::size_t k=1;k<degree;++k) {
                mpq_class next=2*value*current-previous;
                previous=std::move(current);current=std::move(next);
            }
            return current;
        };
        check(idealAt(t2,sample)==2*sample*sample-1,"T_2 fast-doubling identity is exact");
        check(idealAt(t3,sample)==4*sample*sample*sample-3*sample,
              "T_3 fast-doubling identity is exact");
        check(idealAt(t8,sample)==referenceT(8,sample),"T_8 fast-doubling identity is exact");
        check(idealAt(t255,sample)==referenceT(255,sample)
              &&idealAt(t256,sample)==referenceT(256,sample),
              "T_255 and T_256 DAGs equal the independent exact recurrence");
        std::vector<std::size_t> depth(chebBuilder.nodes.size());
        for(std::size_t node=0;node<chebBuilder.nodes.size();++node) {
            for(const auto in:chebBuilder.nodes[node].inputs)
                depth[node]=std::max(depth[node],depth[in]);
            if(chebBuilder.nodes[node].operation==EvalRoundOperation::Multiply)++depth[node];
        }
        check(depth[t255]<=8&&depth[t256]<=8,
              "fast-doubling T_255/T_256 multiplication depth is logarithmic");
        for(const auto degree:{std::size_t(2),std::size_t(3),std::size_t(8),
                               std::size_t(255),std::size_t(256)}) {
            const auto bound=chebyshev.exactMagnitudeBound(degree);
            check(bound==referenceT(degree,R)
                  &&abs(referenceT(degree,R))<=bound
                  &&abs(referenceT(degree,-R))<=bound,
                  "exact T_k(R) recurrence covers both domain endpoints");
        }
    }
    {
        const std::string repeating="0."+std::string(399,'3')+"7";
        const mpq_class coefficient=parseExactDecimal(repeating);
        check(!arithmetic::projectUp(mpq_class(coefficient.get_den())),
              "test coefficient denominator exceeds finite binary64 scale range");
        arithmetic::Builder fallback(a,input,options.evaluationKeyNoiseCoefficientSupport);
        const auto fallbackInput=fallback.input(a.scale(input),1,options.inputSemanticError.upperBound);
        arithmetic::PolynomialCompiler compiler(fallback,fallbackInput);
        const auto output=compiler.compile({PolynomialBasis::Monomial,{"0",repeating}},
                                           "finite-scale fallback diagnostic");
        const auto& node=fallback.nodes[output];
        check(node.operation==EvalRoundOperation::MultiplyPlain
              &&node.stage.find("finite coefficient 1")!=std::string::npos,
              "denominator above binary64 enters direct finite-scale coefficient fallback");
        const mpq_class scale(mpz_class(node.constantEncodingScale.numerator),
                              mpz_class(node.constantEncodingScale.denominator));
        const mpz_class encoded(node.encodedConstantInteger);
        const mpq_class represented(mpz_class(node.representedConstantNumerator),
                                    mpz_class(node.representedConstantDenominator));
        const mpq_class delta(mpz_class(node.exactConstantEncodingError.numerator),
                              mpz_class(node.exactConstantEncodingError.denominator));
        check(std::isfinite(node.constantScale)&&node.constantScale>0
              &&scale==mpq_class(node.constantScale)
              &&scale<mpq_class(mpz_class(std::to_string(
                    *std::min_element(node.activePrimes.begin(),node.activePrimes.end()))))
              &&encoded==arithmetic::roundq(coefficient*scale),
              "fallback uses actual active primes, a finite exact dyadic scale and n=round(c*S)");
        check(represented==mpq_class(encoded)/scale
              &&delta==abs(represented-coefficient)&&delta>0
              &&mpz_class(node.centeredHeadroomNumerator)>0,
              "fallback trace retains exact n/S and charges exact |n/S-c|");
        const mpq_class incomingError(options.inputSemanticError.upperBound);
        const mpq_class representation(
            mpz_class(node.exactScaleRepresentationError.numerator),
            mpz_class(node.exactScaleRepresentationError.denominator));
        check(fallback.states[output].E
                  ==incomingError*abs(coefficient)+(1+incomingError)*delta+representation,
              "exact coefficient delta enters the MultiplyPlain semantic-error recurrence");

        const std::string tiny="0."+std::string(399,'0')+"1";
        arithmetic::Builder zeroRounded(a,input,options.evaluationKeyNoiseCoefficientSupport);
        const auto zeroInput=zeroRounded.input(a.scale(input),1,0);
        arithmetic::PolynomialCompiler zeroCompiler(zeroRounded,zeroInput);
        const auto zeroOutput=zeroCompiler.compile({PolynomialBasis::Monomial,{"0",tiny}},
                                                   "zero-rounded fallback diagnostic");
        const auto& zeroNode=zeroRounded.nodes[zeroOutput];
        check(zeroNode.encodedConstantInteger=="0"
              &&mpq_class(mpz_class(zeroNode.exactConstantEncodingError.numerator),
                           mpz_class(zeroNode.exactConstantEncodingError.denominator))
                    ==abs(parseExactDecimal(tiny))
              &&zeroRounded.states[zeroOutput].E>0,
              "zero-rounded coefficient remains an explicit node and charges its full error");

        bool excessiveScaleRejected=false;
        try {
            arithmetic::Builder excessive(a,input,options.evaluationKeyNoiseCoefficientSupport);
            const auto excessiveInput=excessive.input(a.scale(input),1,0);
            excessive.scalar(excessiveInput,1,std::numeric_limits<double>::max(),
                             "excessive coefficient scale diagnostic");
        } catch(const arithmetic::Failure& failure) {
            excessiveScaleRejected=failure.status==S::ScaleScheduleInfeasible
                ||failure.status==S::HeadroomViolation;
        }
        check(excessiveScaleRejected,
              "too-large plaintext scale fails a concrete scale/headroom gate");
    }
    EvalRoundProblem p{1,1./128,1e-4,4};
    auto candidate=makeEvalRoundReferenceCandidate(p,EvalRoundRadix::Binary,EvalRoundExtractionMethod::BinaryQuadraticK1);
    const mpq_class exactRho(p.rho);
    for(const auto& polynomial:candidate.extraction.polynomials)
        for(const auto center:{-1,0,1})checkShift(polynomial.polynomial,center,{-exactRho,mpq_class(0),exactRho});
    check(candidate.extraction.polynomials[0].polynomial.decimalCoefficients==std::vector<std::string>({"1","0","-1"}),"b0 exact shared polynomial representation");
    check(candidate.extraction.polynomials[1].polynomial.decimalCoefficients==std::vector<std::string>({"0","0.5","0.5"}),"b1 exact shared polynomial representation");
    auto compile=[&](const EvalRoundCandidate& c,const EvalRoundProblem& domain){return EvalRoundExecutionCompiler::compile(a,input,planEvalRoundCandidate(domain,c),options);};
    auto good=compile(candidate,p);check(good.certification().status==S::Certified,good.certification().provenance.c_str());
    check(good.nodes().size()==54,"K1 exact common-denominator monomial path remains 54 nodes");
    for(const auto& node:good.nodes())
        check(node.stage.find("Chebyshev")==std::string::npos
              &&node.stage.find("normalizeScale")==std::string::npos,
              "K1 does not enter Chebyshev execution or metadata scale normalization");
    {
        auto changedExecution=candidate;
        auto& polynomial=changedExecution.extraction.polynomials[0];
        polynomial.executionRepresentation=EvalRoundPolynomialExecutionRepresentation{
            {PolynomialBasis::Chebyshev,{"-2047","0","-2048"}},"64",
            "test exact rational conversion"};
        polynomial.executionRepresentation->polynomial.decimalCoefficients[0]="-2046";
        const auto rejected=compile(changedExecution,p);
        check(rejected.certification().status==S::ExtractionNotCertified
              &&rejected.nodes().empty()
              &&rejected.certification().provenance.find("not exactly equivalent")!=std::string::npos,
              "a modified Chebyshev coefficient is rejected before ciphertext arithmetic");
        polynomial.executionRepresentation->polynomial.decimalCoefficients[0]="-2047";
        polynomial.executionRepresentation->variableScaleDecimal="32";
        const auto wrongScale=compile(changedExecution,p);
        check(wrongScale.certification().status==S::ExtractionNotCertified
              &&wrongScale.nodes().empty()
              &&wrongScale.certification().provenance.find("not exactly equivalent")!=std::string::npos,
              "a wrong scaled-Chebyshev variable scale is rejected before ciphertext arithmetic");
    }
    for(auto mutation:{0,1,2,3,4,5,6,7,8,9}) {
        auto bad=candidate;auto& poly=bad.extraction.polynomials[0];
        if(mutation==0)poly.polynomial.decimalCoefficients.clear();
        if(mutation==1)poly.polynomial.decimalCoefficients[0]="Unknown";
        if(mutation==2)poly.approximationError={};
        if(mutation==3)poly.certifiedK=64;
        if(mutation==4)poly.certifiedRho=std::nextafter(p.rho,1.);
        if(mutation==5)poly.proof=EvalRoundPolynomialProof::GridDiagnostic;
        if(mutation==6)poly.polynomial.decimalCoefficients[0]="1.01";
        if(mutation==7)poly.digitIndex=1;
        if(mutation==8)poly.radix=EvalRoundRadix::BalancedTernary;
        if(mutation==9)poly.target=EvalRoundDigitTarget::TernaryRoot;
        check(compile(bad,p).certification().status!=S::Certified,"Unknown/stale/grid polynomial evidence rejected");
    }
    EvalRoundProblem general{64,1./128,1e-4,4};
    auto no=makeEvalRoundReferenceCandidate(general,EvalRoundRadix::Binary,EvalRoundExtractionMethod::PiecewiseReference);
    check(compile(no,general).certification().status==S::ExtractionNotCertified,"K64 reference target not executable");
    no.extraction.method=EvalRoundExtractionMethod::DigitExtract;no.extraction.digitExtractEpsilon=.125;
    check(compile(no,general).certification().status!=S::Certified,"DigitExtract without concrete evidence fails");
    no.extraction.method=EvalRoundExtractionMethod::ExternalPolynomial;
    check(compile(no,general).certification().status!=S::Certified,"External K64 coefficients missing");
    no.extraction.polynomials.resize(no.digits.size());
    for(std::size_t j=0;j<no.digits.size();++j){auto& poly=no.extraction.polynomials[j];poly.digitIndex=j;poly.certifiedK=64;poly.certifiedRho=general.rho;poly.polynomial.decimalCoefficients={"0","1"};poly.verified=true;poly.proof=EvalRoundPolynomialProof::OutwardInterval;poly.provenance="Negative test only: no interval error provided";}
    check(compile(no,general).certification().status!=S::Certified,"K64 coefficients with Unknown interval error fail");

    for(int mismatch=0;mismatch<3;++mismatch){auto bad=no;
        for(auto& poly:bad.extraction.polynomials){poly.approximationError=B(0);if(mismatch==0)poly.certifiedK=1;if(mismatch==1)poly.certifiedRho=1./64;if(mismatch==2)poly.proof=EvalRoundPolynomialProof::GridDiagnostic;}
        check(compile(bad,general).certification().status!=S::Certified,"K64 wrong K/rho/grid certificate rejected");
    }
    // Exercise degree-three machinery with the already-existing v9 cleaner.
    // This is arithmetic layer C only, not an extraction certificate for any K.
    auto cleanerInput=a.encrypt(a.encode({.25}));
    arithmetic::Builder builder(a,cleanerInput,options.evaluationKeyNoiseCoefficientSupport);
    auto start=builder.input(a.scale(cleanerInput),1,options.inputSemanticError.upperBound);
    arithmetic::PolynomialCompiler polynomial(builder,start);
    auto cleanerOutput=polynomial.compile({PolynomialBasis::Monomial,{"0","0","3","-2"}},"generic f2 polynomial diagnostic");
    std::vector<Plain> constants(builder.nodes.size());
    for(std::size_t i=0;i<builder.nodes.size();++i){const auto& n=builder.nodes[i];if(n.operation==EvalRoundOperation::MultiplyPlain||n.operation==EvalRoundOperation::AddPlain){
        std::vector<std::uint64_t> residues;for(auto prime:n.activePrimes){mpz_class residue;auto modulus=mpz_class(std::to_string(prime));mpz_mod(residue.get_mpz_t(),builder.rounded[i].get_mpz_t(),modulus.get_mpz_t());residues.push_back(std::stoull(residue.get_str()));}
        constants[i]=a.encodeScalarRnsAtScaleFor(residues,n.constantScale,cleanerInput,builder.states[i].level);
    }}
    auto cleaned=executeCertifiedArithmetic(a,{cleanerInput},builder.nodes,constants,cleanerOutput,{});
    check(std::abs(a.decodeComplex(a.decrypt(cleaned))[0]-0.15625)<=builder.nodes[cleanerOutput].semanticError.upperBound,"Generic cubic evaluates v9 f2 from coefficients");

    auto external=candidate;external.extraction.method=EvalRoundExtractionMethod::ExternalPolynomial;
    for(std::size_t j=0;j<2;++j){
        auto poly=certifyEvalRoundDigitPolynomial(p,j,candidate.extraction.polynomials[j].polynomial,64);
        check(poly.approximationError.upperBound<=poly.directXApproximationError.upperBound
              && poly.centeredApproximationError.upperBound>=0
              && poly.intervalProofPrecisionBits==384,
              "Centered proof cannot weaken the previous direct-x K1 bound");
        external.extraction.polynomials[j]=poly;external.digits[j].extractionError=poly.approximationError;
        for(int I=-1;I<=1;++I)for(double offset:{-p.rho,0.,p.rho}){
            const mpq_class x=mpq_class(I)+mpq_class(offset);
            mpq_class value=0;for(auto it=poly.polynomial.decimalCoefficients.rbegin();it!=poly.polynomial.decimalCoefficients.rend();++it)value=value*x+parseExactDecimal(*it);
            const auto target=evalRoundIntegerDigits(I,1,EvalRoundRadix::Binary)[j];
            check(abs(value-target)<=mpq_class(poly.approximationError.upperBound),"Whole-cell interval proof covers both endpoints and center");
        }
    }
    auto verified=compile(external,p);check(verified.certification().status==S::Certified,verified.certification().provenance.c_str());
    std::cout<<std::hexfloat<<"K1 locked certificate bits="
             <<verified.integerErrorUpper()<<std::defaultfloat<<'\n';
    check(verified.nodes().size()==54
          &&verified.integerErrorUpper()==0x1.9dc9f4187ff64p-17,
          "K1 outward-interval certificate tightening and 54-node DAG are locked");
    auto missing=external;missing.extraction.polynomials.pop_back();
    check(compile(missing,p).certification().status!=S::Certified,"Missing one polynomial rejects the extractor");
    auto grid=external;grid.extraction.polynomials[0].proof=EvalRoundPolynomialProof::GridDiagnostic;
    check(compile(grid,p).certification().status!=S::Certified,"Grid cannot replace interval evidence");
    auto under=external;under.extraction.polynomials[0].approximationError=B(0);
    check(compile(under,p).certification().status!=S::Certified,"Interval bound is recomputed, not trusted");
    auto changed=external;changed.extraction.polynomials[0].polynomial.decimalCoefficients[0]="10";
    check(compile(changed,p).certification().status!=S::Certified,"Modified exact coefficient invalidates or recomputes the certificate");
    auto incomplete=external;incomplete.extraction.polynomials[0].intervalSubdivisions=1;
    check(compile(incomplete,p).certification().status!=S::Certified,"Incomplete interval partition cannot certify");
    auto staleCentered=external;staleCentered.extraction.polynomials[0].centeredApproximationError=B(0);
    check(compile(staleCentered,p).certification().status!=S::Certified,"Stale centered bound cannot certify");
    double worst=0;
    for(int I=-1;I<=1;++I)for(double offset:{-p.rho,0.,p.rho}){
        const double value=I+offset;auto encrypted=a.encrypt(a.encode({value}));
        std::vector<mpq_class> ideal(verified.nodes().size());std::size_t nodes=0;
        auto out=executeEvalRound(a,encrypted,verified,[&](std::size_t i,const Cipher& ct){
            const auto& n=verified.nodes()[i];auto& v=ideal[i];
            if(n.operation==EvalRoundOperation::Input)v=mpq_class(value);
            else {v=ideal[n.inputs[0]];switch(n.operation){
            case EvalRoundOperation::Multiply:v*=ideal[n.inputs[1]];break;
            case EvalRoundOperation::Add:v+=ideal[n.inputs[1]];break;
            case EvalRoundOperation::Subtract:v-=ideal[n.inputs[1]];break;
            case EvalRoundOperation::MultiplyPlain:v*=mpq_class(mpz_class(n.constantNumerator),mpz_class(n.constantDenominator));break;
            case EvalRoundOperation::AddPlain:v+=mpq_class(mpz_class(n.constantNumerator),mpz_class(n.constantDenominator));break;
            default:break;}}
            auto z=a.decodeComplex(a.decrypt(ct))[0];mpq_class re=mpq_class(z.real())-v,im=mpq_class(z.imag());
            check(re*re+im*im<=mpq_class(n.semanticError.upperBound)*mpq_class(n.semanticError.upperBound),"Every generic polynomial/cleaner node observed <= certified");
            std::uint64_t bits;auto scale=a.scale(ct);std::memcpy(&bits,&scale,8);
            check(bits==n.outputScale.binary64Bits&&a.coeffModulusValues(ct)==n.activePrimes&&a.chainIndex(ct)==n.chainIndex,"Actual schedule matches exact certificate");++nodes;
        });
        const auto observed=std::abs(a.decodeComplex(a.decrypt(out))[0]-double(I));
        check(observed<=verified.integerErrorUpper()&&verified.integerErrorUpper()<=p.requiredIntegerError,"Final supplied error budget");worst=std::max(worst,observed);
        check(nodes==verified.nodes().size(),"Every reachable generic node executed");
        auto ref=evaluateEvalRoundReference(p,external,mpq_class(value),{verified.mathematicalPlan().digits[0].cleaningIterations,verified.mathematicalPlan().digits[1].cleaningIterations});
        check(abs(ref.reconstructed-ideal.back())<mpq_class(1,mpz_class(1)<<100),"External polynomial reference agrees with exact DAG semantics");
    }
    std::cout<<std::setprecision(14)<<"generic interval K1 certified="<<verified.integerErrorUpper()<<" observed="<<worst<<" nodes="<<verified.nodes().size()<<"\n";
    std::cout<<"PASS polynomial candidates, interval target proof, generic execution and K64 negative contracts\n";
}catch(const arithmetic::Failure& e){std::cerr<<e.why<<"\n";return 1;}
 catch(const std::exception& e){std::cerr<<e.what()<<"\n";return 1;}}
