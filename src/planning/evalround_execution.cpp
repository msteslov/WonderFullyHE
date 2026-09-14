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
struct DigitPath { std::vector<std::size_t> outputs; std::vector<mpq_class> errors; };
double finiteStageBound(const mpq_class& value,BootstrapCertificationStatus failure,
                        const std::string& why) {
    const auto projected=projectUp(value);
    if(!projected) throw Failure{failure,why};
    return *projected;
}
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
                const auto recomputed=certifyEvalRoundDigitPolynomial(
                    reference.problem,j,p.polynomial,p.intervalSubdivisions);
                if(p.intervalProofPrecisionBits!=recomputed.intervalProofPrecisionBits
                   ||p.selectedIntervalProofMethod!=recomputed.selectedIntervalProofMethod
                   ||!known(p.directXApproximationError)||!known(p.centeredApproximationError)
                   ||p.directXApproximationError.upperBound<recomputed.directXApproximationError.upperBound
                   ||p.centeredApproximationError.upperBound<recomputed.centeredApproximationError.upperBound)
                    throw Failure{Status::ExtractionNotCertified,"Stale/incomplete direct-x or centered whole-domain proof metadata"};
                verified=recomputed.approximationError;
            } else throw Failure{Status::ExtractionNotCertified,"Unknown/grid-only polynomial evidence is not a whole-domain certificate"};
            if(p.executionRepresentation) {
                const auto& execution=*p.executionRepresentation;
                if(execution.polynomial.basis!=PolynomialBasis::Chebyshev
                   ||execution.variableScaleDecimal.empty()
                   ||execution.exactEquivalenceProvenance.empty())
                    throw Failure{Status::ExtractionNotCertified,
                        "Incomplete scaled-Chebyshev execution representation metadata"};
                EvalModPolynomial converted;
                try {
                    converted=convertScaledChebyshevToMonomial(
                        execution.polynomial,execution.variableScaleDecimal);
                } catch(const std::exception& error) {
                    throw Failure{Status::ExtractionNotCertified,
                        "Invalid scaled-Chebyshev execution representation: "+std::string(error.what())};
                }
                if(!exactPolynomialEqual(converted,p.polynomial))
                    throw Failure{Status::ExtractionNotCertified,
                        "Scaled-Chebyshev execution representation is not exactly equivalent to canonical certified polynomial"};
            }
            if(p.approximationError.upperBound<verified.upperBound)
                throw Failure{Status::ExtractionNotCertified,"Claimed polynomial approximation understates verified whole-domain error"};
            candidate.digits[j].extractionError=p.approximationError;
            candidate.digits[j].cleaningLocalErrors.assign(reference.problem.maxCleaningRoundsPerDigit,bound(0,"Exact reference cleaner only, replaced by backend bounds"));
            candidate.digits[j].reconstructionLocalError=bound(0,"Exact reference reconstruction only, replaced by backend bounds");
        }
        if(planEvalRoundCandidate(reference.problem,candidate).status!=EvalRoundPlanStatus::Certified)
            throw Failure{Status::ExtractionNotCertified,"The actual polynomial candidate does not satisfy the mathematical extraction contract"};
        Builder b(adapter,input,options.evaluationKeyNoiseCoefficientSupport);
        const auto start=b.input(adapter.scale(input),
            mpq_class(integer(reference.problem.K))+q(reference.problem.rho),
            q(options.inputSemanticError.upperBound));
        b.nodes[start].semanticError.provenance=options.inputSemanticError.provenance;
        b.nodes[start].propagatedSemanticError.provenance=options.inputSemanticError.provenance;
        auto conjugate=b.add(Op::Conjugate,{start},"real projection");
        auto sum=b.add(Op::Add,{start,conjugate},"real projection");
        auto x=b.scalar(sum,mpq_class(1,2),2,"real projection");
        PolynomialCompiler polynomials(b,x);
        std::optional<std::size_t> normalized;
        std::unique_ptr<ScaledChebyshevCompiler> chebyshev;
        if(std::any_of(candidate.extraction.polynomials.begin(),candidate.extraction.polynomials.end(),
            [](const auto& polynomial){return polynomial.executionRepresentation.has_value();})) {
            normalized=b.scalar(x,mpq_class(1,64),64,"Chebyshev normalized t=x/64");
            const mpq_class radius=(mpq_class(integer(reference.problem.K))+q(reference.problem.rho))/64;
            chebyshev=std::make_unique<ScaledChebyshevCompiler>(b,*normalized,radius);
        }
        std::vector<std::size_t> extracted;
        for(std::size_t j=0;j<count;++j) {
            const auto& polynomial=candidate.extraction.polynomials[j];
            const auto stage="extraction digit "+std::to_string(j);
            if(polynomial.executionRepresentation)
                extracted.push_back(chebyshev->compile(
                    polynomial.executionRepresentation->polynomial,stage));
            else extracted.push_back(polynomials.compile(polynomial.polynomial,stage));
        }
        data->diagnostics.digits.resize(count);
        for(std::size_t digit=0;digit<count;++digit) {
            auto& diagnostic=data->diagnostics.digits[digit];
            diagnostic.approximationError=candidate.digits[digit].extractionError;
            std::vector<bool> dependency(b.nodes.size());
            std::function<void(std::size_t)> visit=[&](std::size_t node) {
                if(dependency[node])return;
                dependency[node]=true;
                for(const auto inputNode:b.nodes[node].inputs)visit(inputNode);
            };
            visit(extracted[digit]);
            for(std::size_t node=0;node<dependency.size();++node)
                if(dependency[node]&&b.nodes[node].stage.find("Chebyshev")!=std::string::npos)
                    ++diagnostic.chebyshevNodeCount;
        }
        std::vector<DigitPath> paths(count);
        std::optional<Failure> firstBackendLimit;
        for(std::size_t digit=0;digit<count;++digit) {
            const auto extraction=extracted[digit];
            const bool stableRepresentation=
                candidate.extraction.polynomials[digit].executionRepresentation.has_value();
            mpq_class refError(candidate.digits[digit].extractionError.upperBound);
            b.tightenMagnitude(extraction,1+refError);
            if(stableRepresentation)b.compactSemanticError(extraction,
                "384-bit exact dyadic outward bound for compiled extraction semantic error");
            mpq_class error=refError+b.states[extraction].E;
            auto& diagnostic=data->diagnostics.digits[digit];
            diagnostic.backendExtractionError=exactBound(b.states[extraction].E,
                "Actual scaled-Chebyshev ciphertext DAG semantic error relative to canonical p_j");
            diagnostic.initialCleanerError=exactBound(error,
                "a_0=E_approx+E_backend, with E_approx counted exactly once");
            diagnostic.cleanerErrorAfterRounds.push_back(diagnostic.initialCleanerError);
            if(error>1) throw Failure{Status::ErrorBudgetExceeded,
                "Extraction digit "+std::to_string(digit)+" exact arithmetic error exceeds cleaner domain"
                +(b.firstUnprojectableBound?"; first finite exact bound above binary64: "+*b.firstUnprojectableBound:"")};
            candidate.digits[digit].extractionError={
                finiteStageBound(error,Status::ErrorBudgetExceeded,
                    "Extraction digit error exceeds finite planner representation"),
                BootstrapBoundKind::Deterministic,
                "Whole-domain digit polynomial approximation plus compiled input/projection/polynomial arithmetic",{}};
            candidate.digits[digit].cleaningLocalErrors.clear();
            paths[digit].outputs.push_back(extraction); paths[digit].errors.push_back(error);
            for(std::size_t round=0;round<reference.problem.maxCleaningRoundsPerDigit && error<=1;++round) {
                const auto saved=b.nodes.size();
                try {
                    const auto previous=paths[digit].outputs.back();
                    auto out=b.cleaner(previous,"cleaning digit "+std::to_string(digit)+" round "+std::to_string(round),
                        candidate.extraction.polynomials[digit].executionRepresentation.has_value());
                    // Replay this block with zero *incoming arithmetic* error and
                    // magnitude of the actual incoming digit, not the exact target.
                    auto local=b.localBlock(previous,out,1+error);
                    const mpq_class next=stableRepresentation
                        ?dyadicUpper(5*error*error+local):5*error*error+local;
                    if(!projectUp(local)) throw Failure{Status::ErrorBudgetExceeded,
                        "Cleaning digit "+std::to_string(digit)+" round "+std::to_string(round)+" exact local error exceeds planner budget representation"};
                    candidate.digits[digit].cleaningLocalErrors.push_back(bound(local,"All arithmetic nodes of baseline f2(a)=a^2(3-2a), at |a|<=1+incoming digit error"));
                    diagnostic.cleanerLocalErrors.push_back(exactBound(local,
                        "Actual baseline f2 ciphertext block local arithmetic error"));
                    error=next;
                    refError=stableRepresentation
                        ?dyadicUpper(5*refError*refError):5*refError*refError;
                    b.tightenMagnitude(out,1+refError);
                    if(stableRepresentation)b.compactSemanticError(out,
                        "384-bit exact dyadic outward bound after binary cleaner arithmetic");
                    diagnostic.cleanerErrorAfterRounds.push_back(exactBound(error,
                        "Exact a_next<=5*a^2+B_cln recurrence"));
                    paths[digit].outputs.push_back(out); paths[digit].errors.push_back(error);
                } catch(const Failure& f) {
                    b.nodes.resize(saved); b.states.resize(saved); b.rounded.resize(saved);
                    if(f.status!=Status::InsufficientLevels && f.status!=Status::HeadroomViolation && f.status!=Status::ScaleScheduleInfeasible) throw;
                    if(!firstBackendLimit)firstBackendLimit=f;
                    break;
                }
            }
        }
        {
            auto& diagnostic=data->diagnostics;
            diagnostic.constructedNodes=b.nodes.size();
            std::vector<std::size_t> multiplyDepth(b.nodes.size());
            std::optional<mpz_class> minimumHeadroom;
            for(std::size_t node=0;node<b.nodes.size();++node) {
                const auto& trace=b.nodes[node];
                std::size_t depth=0;
                for(const auto inputNode:trace.inputs)
                    depth=std::max(depth,multiplyDepth[inputNode]);
                if(trace.operation==Op::Multiply)++depth;
                multiplyDepth[node]=depth;
                diagnostic.criticalMultiplicativeDepth=std::max(
                    diagnostic.criticalMultiplicativeDepth,depth);
                diagnostic.criticalPathLevelConsumption=std::max(
                    diagnostic.criticalPathLevelConsumption,b.states[node].level);
                const double scale=b.states[node].scale;
                diagnostic.minimumRuntimeScale=diagnostic.minimumRuntimeScale
                    ?std::min(*diagnostic.minimumRuntimeScale,scale):scale;
                diagnostic.maximumRuntimeScale=diagnostic.maximumRuntimeScale
                    ?std::max(*diagnostic.maximumRuntimeScale,scale):scale;
                if(!trace.centeredHeadroomNumerator.empty()) {
                    const mpz_class headroom(trace.centeredHeadroomNumerator);
                    if(!minimumHeadroom||headroom<*minimumHeadroom)minimumHeadroom=headroom;
                }
                switch(trace.operation) {
                case Op::Multiply:++diagnostic.ciphertextMultiplications;break;
                case Op::Relinearize:++diagnostic.relinearizations;break;
                case Op::Rescale:++diagnostic.rescales;break;
                case Op::ModSwitch:++diagnostic.modSwitches;break;
                case Op::MultiplyPlain:++diagnostic.plaintextMultiplications;break;
                default:break;
                }
            }
            if(minimumHeadroom)
                diagnostic.minimumCenteredHeadroomNumerator=minimumHeadroom->get_str();
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
                const mpq_class e=paths[d].errors[d?j:i];
                mpq_class local=(1+e)*(delta+absq(ratio-1)*(weight+delta));
                if(d==0) local+=absq(mpq_class(roundq(-mpq_class(integer(reference.problem.K))*q(out)))/q(out)+integer(reference.problem.K));
                local/=weight;
                if(local>rec[d]) rec[d]=local;
            }
        }
        } else {
            // Uniform arithmetic reservation across reachable cleaning choices.
            // Replay only the exact scalar/add state transitions.  Copying the
            // entire Builder (including every published GMP certificate) for
            // every Cartesian choice made K64 reconstruction super-linear in
            // memory and time without changing this bound.
            const auto multiplyPlainState=[](State input,const mpq_class& k,
                                              double constantScale) {
                const double outputScale=input.scale*constantScale;
                if(!std::isfinite(outputScale)||outputScale<=0)
                    throw Failure{Status::ScaleScheduleInfeasible,
                        "Reconstruction reservation scale overflow/underflow"};
                const mpz_class encoded=roundq(k*q(constantScale));
                const mpq_class delta=absq(mpq_class(encoded)/q(constantScale)-k);
                const mpq_class magnitude=input.M*absq(k);
                mpq_class error=input.E*absq(k)+(input.M+input.E)*delta;
                error+=absq(q(input.scale)*q(constantScale)/q(outputScale)-1)
                    *(magnitude+error);
                input.scale=outputScale; input.M=magnitude; input.E=error;
                return input;
            };
            std::size_t combinations=1;for(const auto& p:paths) { if(combinations>65536/p.outputs.size())throw Failure{Status::ScaleScheduleInfeasible,"Reconstruction search exceeds baseline work limit"};combinations*=p.outputs.size(); }
            for(std::size_t choice=0;choice<combinations;++choice) {
                std::vector<State> roots; roots.reserve(count); auto code=choice;
                for(std::size_t d=0;d<count;++d) {
                    const auto r=code%paths[d].outputs.size(); code/=paths[d].outputs.size();
                    auto state=b.states[paths[d].outputs[r]];
                    state.M=1+paths[d].errors[r]; state.E=0;
                    roots.push_back(std::move(state));
                }
                auto total=roots[0];
                for(std::size_t d=1;d<count;++d) {
                    auto next=roots[d];
                    total.level=next.level=std::max(total.level,next.level);
                    const double a=total.scale,c=next.scale;
                    total=multiplyPlainState(total,1,c);
                    next=multiplyPlainState(next,mpq_class(mpz_class(1)<<d),a);
                    if(bits(total.scale)!=bits(next.scale))
                        throw Failure{Status::ScaleScheduleInfeasible,
                            "Reconstruction reservation addition scale mismatch"};
                    total.M+=next.M; total.E+=next.E;
                }
                const mpq_class shift=-mpq_class(integer(reference.problem.K));
                const mpq_class delta=absq(
                    mpq_class(roundq(shift*q(total.scale)))/q(total.scale)-shift);
                total.M+=absq(shift); total.E+=delta;
                if(total.E>rec[0])rec[0]=total.E;
            }
        }
        for(std::size_t d=0;d<count;++d) {
            if(!projectUp(rec[d])) throw Failure{Status::ErrorBudgetExceeded,
                "Exact reconstruction arithmetic error exceeds planner budget representation"};
            candidate.digits[d].reconstructionLocalError=bound(rec[d],
                "Uniform exact scalar alignment/weight encoding, scale representation and constant-shift bound over available cleaning counts");
            data->diagnostics.digits[d].reconstructionLocalError=exactBound(rec[d],
                "Uniform exact reconstruction-local bound over available cleaning choices");
        }
        auto selected=planEvalRoundCandidate(reference.problem,candidate);
        if(selected.status!=EvalRoundPlanStatus::Certified) throw Failure{
            firstBackendLimit?firstBackendLimit->status:Status::ErrorBudgetExceeded,
            "No certified backend cleaning schedule within levels/scales and requiredIntegerError: "+selected.provenance
            +(firstBackendLimit?"; first concrete backend limit: "+firstBackendLimit->why:"")};
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
