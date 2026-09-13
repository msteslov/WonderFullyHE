#include "certified_arithmetic_internal.hpp"
#include "evalround_polynomial_internal.hpp"
#include "m2424/experimental/evalmod_analysis/evalround_execution.hpp"
#include "../core/evalround_execution_internal.hpp"
#include <gmpxx.h>
#include "m2424/experimental/evalmod_analysis/finite_support_arithmetic.hpp"
#include <seal/util/defines.h>
#include <algorithm>
#include <cmath>
#include <cfenv>
#include <cstring>
#include <stdexcept>

namespace m2424::experimental {
using namespace arithmetic;
namespace {
struct DigitPath { std::vector<std::size_t> outputs; std::vector<double> errors; };
}
BootstrapBound evalRoundBackendKeyNoiseSupport() { return finiteSupportBackendKeyNoise(); }

EvalRoundExecutionPlan EvalRoundExecutionCompiler::compile(SealAdapter& adapter,const Cipher& input,
    const EvalRoundPlan& reference,const EvalRoundExecutionOptions& options) {
    EvalRoundExecutionPlan result;
    auto data=std::make_shared<EvalRoundExecutionPlan::Data>();
    result.data_=data;
    try {
        if(std::fegetround()!=FE_TONEAREST) throw Failure{Status::ScaleScheduleInfeasible,"Compilation requires round-to-nearest binary64"};
        if(reference.radix!=EvalRoundRadix::Binary||reference.extraction.method==EvalRoundExtractionMethod::PiecewiseReference||reference.extraction.method==EvalRoundExtractionMethod::TernaryPhaseReferenceK1)
            throw Failure{Status::ExtractionNotCertified,"Reference-only target is not a certified binary polynomial extractor"};
        if(reference.status!=EvalRoundPlanStatus::Certified)
            throw Failure{Status::ExtractionNotCertified,"Mathematical/reference certificate required before compilation"};
        if(!known(options.inputSemanticError)) throw Failure{Status::RequiredBoundUnavailable,"Deterministic input semantic error with provenance is required"};
        const auto backendNoise=evalRoundBackendKeyNoiseSupport();
        if(!known(options.evaluationKeyNoiseCoefficientSupport) || !known(backendNoise) || options.evaluationKeyNoiseCoefficientSupport.upperBound<backendNoise.upperBound)
            throw Failure{Status::RequiredBoundUnavailable,"Unknown/insufficient evaluation-key coefficient support; no zero local bound substitution"};
        if(!adapter.hasRelinKeys() || !adapter.hasConjugationKey()) throw Failure{Status::MissingEvaluationKeys,"Binary real projection requires conjugation and relinearization keys"};
        if(adapter.info(input).ciphertextSize!=2) throw Failure{Status::InvalidInput,"Two-component input required"};
        EvalRoundCandidate candidate;candidate.id=reference.candidateId+"/SEAL-polynomial-baseline";
        candidate.radix=reference.radix;candidate.cost=reference.cost;candidate.extraction=reference.extraction;candidate.failureEvents=reference.failureEvents;
        const auto count=evalRoundDigitCount(reference.problem.K,reference.radix);
        if(!count||candidate.extraction.polynomials.size()!=count)
            throw Failure{Status::ExtractionNotCertified,"One concrete polynomial with whole-domain evidence per digit is required"};
        candidate.digits.resize(count);
        for(std::size_t j=0;j<count;++j) {
            const auto& p=candidate.extraction.polynomials[j];
            if(!p.verified||p.provenance.empty()||p.radix!=reference.radix||p.digitIndex!=j||
               p.certifiedK!=reference.problem.K||bits(p.certifiedRho)!=bits(reference.problem.rho)||
               p.target!=EvalRoundDigitTarget::BinaryOffsetDigit||!known(p.approximationError)||!p.failureEventIds.empty())
                throw Failure{Status::ExtractionNotCertified,"Polynomial/domain/target binding or deterministic approximation evidence missing"};
            BootstrapBound verified;
            if(p.proof==EvalRoundPolynomialProof::BinaryQuadraticIdentity) {
                if(reference.problem.K!=1)throw Failure{Status::ExtractionNotCertified,"Quadratic identity proof covers only K=1"};
                auto exact=makeEvalRoundReferenceCandidate(reference.problem,reference.radix,EvalRoundExtractionMethod::BinaryQuadraticK1);
                const auto& coefficients=exact.extraction.polynomials[j].polynomial;
                if(p.polynomial.basis!=coefficients.basis||p.polynomial.decimalCoefficients.size()!=coefficients.decimalCoefficients.size())
                    throw Failure{Status::ExtractionNotCertified,"Polynomial differs from its analytic identity proof"};
                for(std::size_t k=0;k<coefficients.decimalCoefficients.size();++k)
                    if(parseExactDecimal(p.polynomial.decimalCoefficients[k])!=parseExactDecimal(coefficients.decimalCoefficients[k]))
                        throw Failure{Status::ExtractionNotCertified,"Changed polynomial coefficients invalidate the identity proof"};
                verified=exact.digits[j].extractionError;
            } else if(p.proof==EvalRoundPolynomialProof::OutwardInterval) {
                verified=certifyEvalRoundDigitPolynomial(reference.problem,j,p.polynomial,p.intervalSubdivisions).approximationError;
            } else throw Failure{Status::ExtractionNotCertified,"Unknown/grid-only polynomial evidence is not a whole-domain certificate"};
            if(p.approximationError.upperBound<verified.upperBound)
                throw Failure{Status::ExtractionNotCertified,"Claimed polynomial approximation understates verified whole-domain error"};
            candidate.digits[j].extractionError=p.approximationError;
            candidate.digits[j].cleaningLocalErrors.assign(reference.problem.maxCleaningRoundsPerDigit,bound(0,"Exact reference cleaner only, replaced by backend bounds"));
            candidate.digits[j].reconstructionLocalError=bound(0,"Exact reference reconstruction only, replaced by backend bounds");
        }
        if(planEvalRoundCandidate(reference.problem,candidate).status!=EvalRoundPlanStatus::Certified)
            throw Failure{Status::ExtractionNotCertified,"The actual polynomial candidate does not satisfy the mathematical extraction contract"};
        Builder b(adapter,input,options.evaluationKeyNoiseCoefficientSupport);
        const auto start=b.input(adapter.scale(input),up(integer(reference.problem.K)+q(reference.problem.rho)),options.inputSemanticError.upperBound);
        b.nodes[start].semanticError.provenance=options.inputSemanticError.provenance;
        b.nodes[start].propagatedSemanticError.provenance=options.inputSemanticError.provenance;
        auto conjugate=b.add(Op::Conjugate,{start},"real projection");
        auto sum=b.add(Op::Add,{start,conjugate},"real projection");
        auto x=b.scalar(sum,mpq_class(1,2),2,"real projection");
        PolynomialCompiler polynomials(b,x);
        std::vector<std::size_t> extracted;
        for(std::size_t j=0;j<count;++j)extracted.push_back(polynomials.compile(candidate.extraction.polynomials[j].polynomial,"extraction digit "+std::to_string(j)));
        std::vector<DigitPath> paths(count);
        bool levelLimited=false;
        for(std::size_t digit=0;digit<count;++digit) {
            const auto extraction=extracted[digit];
            double refError=candidate.digits[digit].extractionError.upperBound;
            double error=up(q(refError)+b.states[extraction].E);
            candidate.digits[digit].extractionError=bound(q(error),"Whole-domain digit polynomial approximation plus compiled input/projection/polynomial arithmetic");
            candidate.digits[digit].cleaningLocalErrors.clear();
            paths[digit].outputs.push_back(extraction); paths[digit].errors.push_back(error);
            b.tightenMagnitude(extraction,up(1+q(refError)));
            for(std::size_t round=0;round<reference.problem.maxCleaningRoundsPerDigit && error<=1;++round) {
                const auto saved=b.nodes.size();
                try {
                    const auto previous=paths[digit].outputs.back();
                    auto out=b.cleaner(previous,"cleaning digit "+std::to_string(digit)+" round "+std::to_string(round));
                    // Replay this block with zero *incoming arithmetic* error and
                    // magnitude of the actual incoming digit, not the exact target.
                    auto local=b.localBlock(previous,out,up(1+q(error)));
                    candidate.digits[digit].cleaningLocalErrors.push_back(bound(q(local),"All arithmetic nodes of baseline f2(a)=a^2(3-2a), at |a|<=1+incoming digit error"));
                    error=evalRoundCleaningErrorUpper(EvalRoundRadix::Binary,error,local);
                    refError=evalRoundCleaningErrorUpper(EvalRoundRadix::Binary,refError,0);
                    b.tightenMagnitude(out,up(1+q(refError)));
                    paths[digit].outputs.push_back(out); paths[digit].errors.push_back(error);
                } catch(const Failure& f) {
                    b.nodes.resize(saved); b.states.resize(saved); b.rounded.resize(saved);
                    if(f.status!=Status::InsufficientLevels && f.status!=Status::HeadroomViolation && f.status!=Status::ScaleScheduleInfeasible) throw;
                    levelLimited=true; break;
                }
            }
        }
        // Arithmetic alignment uses plaintext 1 and 2 at the other branch's
        // actual dyadic scale. The products have identical binary64 bits.
        // Bound reconstruction uniformly over all retained round counts so the
        // PR-2 dynamic program can still choose the minimum cleaning count.
        std::vector<mpq_class> rec(count,0);
        if(count==2) {
        for(std::size_t i=0;i<paths[0].outputs.size();++i) for(std::size_t j=0;j<paths[1].outputs.size();++j) {
            const double a=b.states[paths[0].outputs[i]].scale, c=b.states[paths[1].outputs[j]].scale;
            const double out=a*c;
            if(!std::isfinite(out)) throw Failure{Status::ScaleScheduleInfeasible,"Reconstruction scale overflow"};
            const mpq_class ratio=q(a)*q(c)/q(out);
            for(std::size_t d=0;d<2;++d) {
                const double cs=d?a:c; const mpq_class weight=d?2:1;
                const mpq_class delta=absq(mpq_class(roundq(weight*q(cs)))/q(cs)-weight);
                const double e=paths[d].errors[d?j:i];
                mpq_class local=(1+q(e))*(delta+absq(ratio-1)*(weight+delta));
                if(d==0) local+=absq(mpq_class(roundq(-mpq_class(integer(reference.problem.K))*q(out)))/q(out)+integer(reference.problem.K));
                local/=weight;
                if(local>rec[d]) rec[d]=local;
            }
        }
        } else {
            // Uniform arithmetic reservation across reachable cleaning choices.
            // Reuse Builder for exact scalar/scale error, never a grid estimate.
            std::size_t combinations=1;for(const auto& p:paths) { if(combinations>65536/p.outputs.size())throw Failure{Status::ScaleScheduleInfeasible,"Reconstruction search exceeds baseline work limit"};combinations*=p.outputs.size(); }
            for(std::size_t choice=0;choice<combinations;++choice) {
                auto temporary=b;std::vector<std::size_t> roots;auto code=choice;
                for(std::size_t d=0;d<count;++d) {const auto r=code%paths[d].outputs.size();code/=paths[d].outputs.size();auto root=paths[d].outputs[r];roots.push_back(root);temporary.states[root].M=q(up(1+q(paths[d].errors[r])));temporary.states[root].E=0;}
                auto total=roots[0];
                for(std::size_t d=1;d<count;++d) {auto next=roots[d];if(temporary.states[total].level<temporary.states[next].level)total=temporary.alignLevel(total,next,"reconstruction reserve");if(temporary.states[next].level<temporary.states[total].level)next=temporary.alignLevel(next,total,"reconstruction reserve");const auto a=temporary.states[total].scale,c=temporary.states[next].scale;total=temporary.scalar(total,1,c,"reconstruction reserve");next=temporary.scalar(next,mpq_class(mpz_class(1)<<d),a,"reconstruction reserve");total=temporary.add(Op::Add,{total,next},"reconstruction reserve");}
                total=temporary.plus(total,-mpq_class(integer(reference.problem.K)),"reconstruction reserve shift");
                if(temporary.states[total].E>rec[0])rec[0]=temporary.states[total].E;
            }
        }
        for(std::size_t d=0;d<count;++d) candidate.digits[d].reconstructionLocalError=bound(rec[d],
            "Uniform exact scalar alignment/weight encoding, scale representation and constant-shift bound over available cleaning counts");
        auto selected=planEvalRoundCandidate(reference.problem,candidate);
        if(selected.status!=EvalRoundPlanStatus::Certified) throw Failure{levelLimited?Status::InsufficientLevels:Status::ErrorBudgetExceeded,
            "No certified backend cleaning schedule within levels/scales and requiredIntegerError: "+selected.provenance};
        auto output=paths[0].outputs[selected.digits[0].cleaningIterations];
        for(std::size_t d=1;d<count;++d) {
            auto next=paths[d].outputs[selected.digits[d].cleaningIterations];
            if(b.states[output].level<b.states[next].level)output=b.alignLevel(output,next,"reconstruction");
            if(b.states[next].level<b.states[output].level)next=b.alignLevel(next,output,"reconstruction");
            const auto a=b.states[output].scale,c=b.states[next].scale;
            output=b.scalar(output,1,c,"reconstruction accumulated digits");
            next=b.scalar(next,mpq_class(mpz_class(1)<<d),a,"reconstruction weighted digit "+std::to_string(d));
            output=b.add(Op::Add,{output,next},"reconstruction");
        }
        output=b.plus(output,-mpq_class(integer(reference.problem.K)),"reconstruction sum(2^j*b_j)-K");
        // Remove diagnostic search nodes. Only the immutable reachable DAG can execute.
        std::vector<bool> used(b.nodes.size());
        std::function<void(std::size_t)> visit=[&](std::size_t n) { if(used[n]) return; used[n]=true; for(auto j:b.nodes[n].inputs) visit(j); };
        visit(output); std::vector<std::size_t> mapping(b.nodes.size());
        for(std::size_t i=0;i<b.nodes.size();++i) if(used[i]) {
            auto n=b.nodes[i]; for(auto& j:n.inputs) j=mapping[j];
            mapping[i]=data->nodes.size(); data->nodes.push_back(n); data->constants.emplace_back();
            if(n.operation==Op::MultiplyPlain || n.operation==Op::AddPlain) {
                std::vector<std::uint64_t> residues;
                for(auto prime:n.activePrimes) {
                    mpz_class r; mpz_mod(r.get_mpz_t(),b.rounded[i].get_mpz_t(),integer(prime).get_mpz_t());
                    residues.push_back(std::stoull(r.get_str()));
                }
                data->constants.back()=adapter.encodeScalarRnsAtScaleFor(residues,n.constantScale,input,b.states[i].level);
            }
        }
        data->output=mapping[output]; data->fingerprint=adapter.contextFingerprint(); data->mathematicalPlan=std::move(selected);
        data->certification={Status::Certified,"EvalRoundExecution","Generic exact-decimal polynomial DAG; exact SEAL prime/dyadic schedule; deterministic finite-support arithmetic and v9 cleaning/reconstruction"};
    } catch(const Failure& f) { data->certification={f.status,"EvalRoundExecution",f.why}; }
      catch(const std::exception& e) { data->certification={Status::InvalidInput,"EvalRoundExecution",e.what()}; }
    return result;
}
} // namespace m2424::experimental
