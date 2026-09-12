#include "certified_arithmetic_internal.hpp"
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
        if(reference.extraction.method!=EvalRoundExtractionMethod::BinaryQuadraticK1 || reference.radix!=EvalRoundRadix::Binary || reference.problem.K!=1)
            throw Failure{Status::ExtractionNotCertified,"Reference-only or unsupported extractor: executable PR-3 path is BinaryQuadraticK1; exact phase is not a polynomial"};
        if(reference.status!=EvalRoundPlanStatus::Certified)
            throw Failure{Status::ExtractionNotCertified,"Mathematical/reference certificate required before compilation"};
        if(!known(options.inputSemanticError)) throw Failure{Status::RequiredBoundUnavailable,"Deterministic input semantic error with provenance is required"};
        const auto backendNoise=evalRoundBackendKeyNoiseSupport();
        if(!known(options.evaluationKeyNoiseCoefficientSupport) || !known(backendNoise) || options.evaluationKeyNoiseCoefficientSupport.upperBound<backendNoise.upperBound)
            throw Failure{Status::RequiredBoundUnavailable,"Unknown/insufficient evaluation-key coefficient support; no zero local bound substitution"};
        if(!adapter.hasRelinKeys() || !adapter.hasConjugationKey()) throw Failure{Status::MissingEvaluationKeys,"Binary real projection requires conjugation and relinearization keys"};
        if(adapter.info(input).ciphertextSize!=2) throw Failure{Status::InvalidInput,"Two-component input required"};
        auto candidate=makeEvalRoundReferenceCandidate(reference.problem,EvalRoundRadix::Binary,
            EvalRoundExtractionMethod::BinaryQuadraticK1,reference.cost);
        // Validate even if a caller has forged fields of the mutable reference plan.
        if(planEvalRoundCandidate(reference.problem,candidate).status!=EvalRoundPlanStatus::Certified)
            throw Failure{Status::ExtractionNotCertified,"Reference problem does not certify the actual quadratic extractor"};
        candidate.id="BinaryQuadraticK1/SEAL-baseline";
        Builder b(adapter,input,options.evaluationKeyNoiseCoefficientSupport);
        const auto start=b.input(adapter.scale(input),up(1+q(reference.problem.rho)),options.inputSemanticError.upperBound);
        b.nodes[start].semanticError.provenance=options.inputSemanticError.provenance;
        b.nodes[start].propagatedSemanticError.provenance=options.inputSemanticError.provenance;
        auto conjugate=b.add(Op::Conjugate,{start},"real projection");
        auto sum=b.add(Op::Add,{start,conjugate},"real projection");
        auto x=b.scalar(sum,mpq_class(1,2),2,"real projection");
        auto square=b.mul(x,x,"extraction");
        auto squareReduced=b.reduce(square,"extraction b0");
        auto b0=b.plus(b.scalar(squareReduced,-1,1,"extraction b0"),1,"extraction b0");
        auto linear=b.scalar(x,1,b.states[x].scale,"extraction b1");
        auto numerator=b.add(Op::Add,{square,linear},"extraction b1");
        auto b1=b.reduce(b.scalar(numerator,mpq_class(1,2),2,"extraction b1"),"extraction b1");
        DigitPath paths[2];
        bool levelLimited=false;
        for(std::size_t digit=0;digit<2;++digit) {
            const auto extraction=digit?b1:b0;
            double refError=candidate.digits[digit].extractionError.upperBound;
            double error=up(q(refError)+b.states[extraction].E);
            candidate.digits[digit].extractionError=bound(q(error),"Whole-domain quadratic approximation plus compiled input/projection/extraction arithmetic");
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
        mpq_class rec[2]={0,0};
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
                if(d==0) local+=absq(mpq_class(roundq(-q(out)))/q(out)+1);
                local/=weight;
                if(local>rec[d]) rec[d]=local;
            }
        }
        for(std::size_t d=0;d<2;++d) candidate.digits[d].reconstructionLocalError=bound(rec[d],
            "Uniform exact scalar alignment/weight encoding, scale representation and constant-shift bound over available cleaning counts");
        auto selected=planEvalRoundCandidate(reference.problem,candidate);
        if(selected.status!=EvalRoundPlanStatus::Certified) throw Failure{levelLimited?Status::InsufficientLevels:Status::ErrorBudgetExceeded,
            "No certified backend cleaning schedule within levels/scales and requiredIntegerError: "+selected.provenance};
        auto a=paths[0].outputs[selected.digits[0].cleaningIterations];
        auto c=paths[1].outputs[selected.digits[1].cleaningIterations];
        if(b.states[a].level<b.states[c].level) a=b.alignLevel(a,c,"reconstruction");
        if(b.states[c].level<b.states[a].level) c=b.alignLevel(c,a,"reconstruction");
        const auto as=b.states[a].scale,cs=b.states[c].scale;
        a=b.scalar(a,1,cs,"reconstruction b0"); c=b.scalar(c,2,as,"reconstruction 2*b1");
        auto output=b.plus(b.add(Op::Add,{a,c},"reconstruction"),-1,"reconstruction I=b0+2*b1-1");
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
        data->certification={Status::Certified,"EvalRoundExecution","BinaryQuadraticK1 polynomial DAG; exact SEAL prime/dyadic schedule; deterministic finite-support arithmetic and v9 cleaning/reconstruction"};
    } catch(const Failure& f) { data->certification={f.status,"EvalRoundExecution",f.why}; }
      catch(const std::exception& e) { data->certification={Status::InvalidInput,"EvalRoundExecution",e.what()}; }
    return result;
}
} // namespace m2424::experimental
