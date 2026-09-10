#include "m2424/experimental/evalmod_analysis/evalround_execution.hpp"
#include "../core/evalround_execution_internal.hpp"
#include <gmpxx.h>
#include <seal/util/defines.h>
#include <algorithm>
#include <cmath>
#include <cfenv>
#include <cstring>
#include <stdexcept>

namespace m2424::experimental {
namespace {
using Op=EvalRoundOperation;
using Status=BootstrapCertificationStatus;
struct Failure { Status status; std::string why; };
mpq_class q(double x) { return mpq_class(x); }
mpq_class rational(const mpz_class& n,const mpz_class& d) {
    mpq_class result(n,d); result.canonicalize(); return result;
}
mpq_class absq(mpq_class x) { return x<0 ? -x : x; }
mpz_class integer(std::uint64_t x) { return mpz_class(std::to_string(x)); }
mpz_class ceilq(const mpq_class& x) { mpz_class r; mpz_cdiv_q(r.get_mpz_t(),x.get_num_mpz_t(),x.get_den_mpz_t()); return r; }
mpz_class roundq(mpq_class x) {
    const bool neg=x<0; if(neg) x=-x;
    mpz_class r=ceilq(x-mpq_class(1,2)); // ties toward zero is also <= 1/2
    return neg ? -r:r;
}
double up(const mpq_class& x) {
    double d=x.get_d();
    if(!std::isfinite(d)) throw Failure{Status::RequiredBoundUnavailable,"Arithmetic bound overflow"};
    if(q(d)<x) d=std::nextafter(d,INFINITY);
    return d;
}
std::uint64_t bits(double x) { std::uint64_t b; std::memcpy(&b,&x,sizeof b); return b; }
EvalRoundExactScale exact(const mpq_class& x,double runtime) { return {x.get_num().get_str(),x.get_den().get_str(),bits(runtime)}; }
BootstrapBound bound(const mpq_class& x,const std::string& why) { return {up(x),BootstrapBoundKind::Deterministic,why,{}}; }
bool known(const BootstrapBound& b) { return b.kind==BootstrapBoundKind::Deterministic && std::isfinite(b.upperBound) && b.upperBound>=0 && !b.provenance.empty() && b.failureEventIds.empty(); }
struct State { double scale{}; std::size_t level{},components{2}; mpq_class M,E; };
struct Calculation { State state; mpq_class propagated,local,encoding,representation; };
struct Builder {
    std::vector<EvalRoundExecutionNode> nodes;
    std::vector<State> states;
    std::vector<mpz_class> rounded;
    std::vector<std::uint64_t> primes;
    std::uint64_t special;
    std::size_t degree,initialIndex;
    mpq_class keyNoise;
    std::string keyProvenance;
    Builder(const SealAdapter& a,const Cipher& c,const BootstrapBound& noise):
        primes(a.coeffModulusValues(c)),special(a.specialKeyModulusValue()),
        degree(2*a.slotCount()),initialIndex(a.chainIndex(c)),keyNoise(q(noise.upperBound)),
        keyProvenance(noise.provenance+"; N="+std::to_string(degree)+"; secret coefficient support=1; P="+std::to_string(special)+"; B_key="+std::to_string(noise.upperBound)) {}
    mpq_class divideRound(std::size_t components,double scale) const {
        mpz_class power=1,sum=0;
        for(std::size_t i=0;i<components;++i) { sum+=power; power*=degree; }
        return rational(integer(degree)*sum,2)/q(scale);
    }
    mpq_class keySwitch(std::size_t level,double scale) const {
        mpz_class sum=0; for(std::size_t i=0;i<primes.size()-level;++i) sum+=integer(primes[i]-1);
        return rational(integer(degree)*integer(degree)*sum,integer(special))*keyNoise/q(scale)+divideRound(2,scale);
    }
    Calculation calculate(std::size_t index,const std::vector<State>& s) const {
        const auto& n=nodes[index];
        Calculation c; c.state=states[index];
        if(n.operation==Op::Input) { c.propagated=c.state.E; return c; }
        const auto& a=s[n.inputs[0]];
        mpq_class M=a.M,E=a.E,local=0;
        const mpq_class outScale(c.state.scale);
        mpq_class arithmetic(a.scale);
        if(n.operation==Op::Multiply) {
            const auto& b=s[n.inputs[1]];
            M=a.M*b.M; E=a.M*b.E+b.M*a.E+a.E*b.E;
            arithmetic=q(a.scale)*q(b.scale);
        } else if(n.operation==Op::Add || n.operation==Op::Subtract) {
            const auto& b=s[n.inputs[1]]; M=a.M+b.M; E=a.E+b.E;
        } else if(n.operation==Op::MultiplyPlain || n.operation==Op::AddPlain) {
            const mpq_class k=rational(mpz_class(n.constantNumerator),mpz_class(n.constantDenominator));
            const mpq_class delta=absq(mpq_class(rounded[index])/q(n.constantScale)-k);
            c.encoding=delta;
            if(n.operation==Op::MultiplyPlain) {
                M=a.M*absq(k); E=a.E*absq(k); local=(a.M+a.E)*delta;
                arithmetic=q(a.scale)*q(n.constantScale);
            } else { M=a.M+absq(k); local=delta; }
        } else if(n.operation==Op::Relinearize) {
            local=(a.components-2)*keySwitch(a.level,c.state.scale);
        } else if(n.operation==Op::Conjugate) {
            local=keySwitch(a.level,c.state.scale);
        } else if(n.operation==Op::Rescale) {
            arithmetic=q(a.scale)/integer(primes[primes.size()-a.level-1]);
            local=divideRound(a.components,c.state.scale);
        }
        // SEAL stores a binary64 scale. This ratio is an arithmetic error,
        // including prime-to-double rounding during CKKS rescale.
        const mpq_class ratio=arithmetic/outScale;
        c.representation=absq(ratio-1)*(M+E+local);
        c.propagated=E; c.local=local+c.representation;
        c.state.M=q(up(M)); c.state.E=q(up(E+c.local));
        return c;
    }
    void publish(std::size_t i,const Calculation& c) {
        auto& n=nodes[i]; states[i]=c.state; n.ciphertextComponents=c.state.components;
        const std::string prefix="EvalRound "+n.stage+" node "+std::to_string(i)+": ";
        n.idealMagnitude=bound(c.state.M,prefix+"triangle/product norms; v9 digit magnitude at stage boundaries");
        n.propagatedSemanticError=bound(c.propagated,prefix+"semantic error propagation through exact polynomial");
        n.localArithmeticError=bound(c.local,prefix+"exact scalar encoding, finite-support key switching, component divide-round, exact scale ratio");
        if(n.requiredKey!=EvalRoundEvaluationKey::None) n.localArithmeticError.provenance+="; "+keyProvenance;
        if(n.operation==Op::Rescale) n.localArithmeticError.provenance+="; N="+std::to_string(degree)+"; secret coefficient support=1; components="+std::to_string(c.state.components);
        n.semanticError=bound(c.state.E,prefix+"propagated plus local arithmetic error");
        n.constantEncodingError=bound(c.encoding,prefix+"exact rational |roundedInteger/encodingScale-constant|");
        n.scaleRepresentationError=bound(c.representation,prefix+"exact rational scale ratio against binary64 metadata");
        mpz_class Q=1; for(auto prime:n.activePrimes) Q*=integer(prime);
        mpz_class margin=Q-2*ceilq(q(c.state.scale)*(c.state.M+c.state.E));
        if(margin<=0) throw Failure{Status::HeadroomViolation,prefix+"centered no-wrap proof unavailable"};
        n.centeredHeadroomNumerator=margin.get_str();
        n.centeredHeadroomProvenance=prefix+"exact Q/2-ceil(runtimeScale*(M+E)); inverse canonical embedding coefficient norm <= slot sup norm";
    }
    std::size_t input(double scale,double M,double E) {
        EvalRoundExecutionNode n; n.operation=Op::Input; n.stage="input";
        n.chainIndex=initialIndex; n.activePrimes=primes;
        n.outputScale=n.arithmeticScale=exact(q(scale),scale);
        nodes.push_back(n); states.push_back({scale,0,2,q(M),q(E)}); rounded.emplace_back(0);
        publish(0,calculate(0,states)); return 0;
    }
    std::size_t add(Op op,std::vector<std::size_t> inputs,std::string stage,mpq_class constant=0,double constantScale=1) {
        const auto a=states.at(inputs[0]); State s=a;
        EvalRoundExecutionNode n; n.operation=op; n.inputs=inputs; n.stage=std::move(stage);
        mpq_class arithmetic=q(a.scale); mpz_class encoded=0;
        if(op==Op::Multiply || op==Op::Add || op==Op::Subtract) {
            const auto& b=states.at(inputs[1]);
            if(a.level!=b.level) throw Failure{Status::ScaleScheduleInfeasible,"Unaligned modulus"};
            if(op==Op::Multiply) { s.scale=a.scale*b.scale; arithmetic=q(a.scale)*q(b.scale); s.components=a.components+b.components-1; }
            else {
                if(bits(a.scale)!=bits(b.scale)) throw Failure{Status::ScaleScheduleInfeasible,"Addition requires identical scale bits"};
                s.components=std::max(a.components,b.components);
            }
        } else if(op==Op::MultiplyPlain || op==Op::AddPlain) {
            if(op==Op::AddPlain) constantScale=a.scale;
            n.constantNumerator=constant.get_num().get_str(); n.constantDenominator=constant.get_den().get_str(); n.constantScale=constantScale;
            encoded=roundq(constant*q(constantScale));
            if(op==Op::MultiplyPlain) { s.scale=a.scale*constantScale; arithmetic=q(a.scale)*q(constantScale); }
        } else if(op==Op::Rescale) {
            if(a.level+1>=primes.size()) throw Failure{Status::InsufficientLevels,"No prime remains for rescale"};
            const auto p=primes[primes.size()-a.level-1];
            s.level++; s.scale=a.scale/static_cast<double>(p); arithmetic=q(a.scale)/integer(p);
        } else if(op==Op::Relinearize) { s.components=2; n.requiredKey=EvalRoundEvaluationKey::Relinearization; }
        else if(op==Op::Conjugate) { n.requiredKey=EvalRoundEvaluationKey::Conjugation; }
        else if(op==Op::ModSwitch) {
            s.level=states.at(inputs[1]).level;
            if(s.level<a.level) throw Failure{Status::ScaleScheduleInfeasible,"Modulus raising prohibited"};
        }
        if(!std::isfinite(s.scale) || s.scale<=0) throw Failure{Status::ScaleScheduleInfeasible,"Scale overflow/underflow"};
        n.activePrimes.assign(primes.begin(),primes.end()-s.level);
        n.chainIndex=initialIndex-s.level;
        for(auto i:inputs) n.inputScales.push_back(exact(q(states[i].scale),states[i].scale));
        n.outputScale=exact(q(s.scale),s.scale); n.arithmeticScale=exact(arithmetic,s.scale);
        // SEAL requires scale < modulus (also before relinearize/rescale).
        mpz_class Q=1; for(auto p:n.activePrimes) Q*=integer(p);
        if(q(s.scale)>=mpq_class(Q) || absq(mpq_class(encoded))*2>=mpq_class(Q))
            throw Failure{Status::HeadroomViolation,"Scale or plaintext exceeds active modulus"};
        const auto i=nodes.size(); nodes.push_back(n); states.push_back(s); rounded.push_back(encoded);
        publish(i,calculate(i,states)); return i;
    }
    std::size_t mul(std::size_t a,std::size_t b,const std::string& stage) { return add(Op::Multiply,{a,b},stage); }
    std::size_t scalar(std::size_t a,mpq_class k,double scale,const std::string& stage) { return add(Op::MultiplyPlain,{a},stage,k,scale); }
    std::size_t plus(std::size_t a,mpq_class k,const std::string& stage) { return add(Op::AddPlain,{a},stage,k); }
    std::size_t reduce(std::size_t a,const std::string& stage) { return add(Op::Rescale,{add(Op::Relinearize,{a},stage)},stage); }
    std::size_t alignLevel(std::size_t a,std::size_t b,const std::string& stage) {
        return states[a].level==states[b].level ? a : add(Op::ModSwitch,{a,b},stage);
    }
    std::size_t cleaner(std::size_t a,const std::string& stage) {
        // Baseline: a*a -> relinearize -> rescale, then *(3-2a) ->
        // relinearize -> rescale. No deferred/thrifty tensor evaluation.
        auto square=reduce(mul(a,a,stage),stage);
        auto linear=plus(scalar(a,-2,1,stage),3,stage);
        linear=alignLevel(linear,square,stage);
        return reduce(mul(square,linear,stage),stage);
    }
    double localBlock(std::size_t first,std::size_t last,double inputMagnitude) const {
        auto s=states; s[first].M=q(inputMagnitude); s[first].E=0;
        for(auto i=first+1;i<=last;++i) s[i]=calculate(i,s).state;
        return up(s[last].E);
    }
    void tightenMagnitude(std::size_t node,double magnitude) {
        states[node].M=q(magnitude);
        nodes[node].idealMagnitude=bound(q(magnitude),"v9 binary digit recurrence bounds exact polynomial ideal magnitude");
        mpz_class Q=1; for(auto prime:nodes[node].activePrimes) Q*=integer(prime);
        const mpz_class margin=Q-2*ceilq(q(states[node].scale)*(states[node].M+states[node].E));
        if(margin<=0) throw Failure{Status::HeadroomViolation,"Digit boundary headroom"};
        nodes[node].centeredHeadroomNumerator=margin.get_str();
    }
};
struct DigitPath { std::vector<std::size_t> outputs; std::vector<double> errors; };
}
BootstrapBound evalRoundBackendKeyNoiseSupport() {
#ifdef SEAL_USE_GAUSSIAN_NOISE
    return {INFINITY,BootstrapBoundKind::Unknown,"Gaussian evaluation-key finite support has not been certified",{}};
#else
    return {21,BootstrapBoundKind::Deterministic,"SEAL sample_poly_cbd: difference of two 21-bit Hamming weights",{}};
#endif
}

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
