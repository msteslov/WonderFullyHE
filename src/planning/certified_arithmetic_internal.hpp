#pragma once
#include "m2424/evalround_execution.hpp"
#include "m2424/experimental/evalmod_analysis/finite_support_arithmetic.hpp"
#include <gmpxx.h>
#include <cmath>
#include <cstring>
#include <stdexcept>
namespace m2424::experimental::arithmetic {

using Op=EvalRoundOperation;
using Status=BootstrapCertificationStatus;
struct Failure { Status status; std::string why; };
inline mpq_class q(double x) { return mpq_class(x); }
inline mpq_class rational(const mpz_class& n,const mpz_class& d) {
    mpq_class result(n,d); result.canonicalize(); return result;
}
inline mpq_class absq(mpq_class x) { return x<0 ? -x : x; }
inline mpz_class integer(std::uint64_t x) { return mpz_class(std::to_string(x)); }
inline mpz_class ceilq(const mpq_class& x) { mpz_class r; mpz_cdiv_q(r.get_mpz_t(),x.get_num_mpz_t(),x.get_den_mpz_t()); return r; }
inline mpz_class roundq(mpq_class x) {
    const bool neg=x<0; if(neg) x=-x;
    mpz_class r=ceilq(x-mpq_class(1,2)); // ties toward zero is also <= 1/2
    return neg ? -r:r;
}
inline double up(const mpq_class& x) {
    double d=x.get_d();
    if(!std::isfinite(d)) throw Failure{Status::RequiredBoundUnavailable,"Arithmetic bound overflow"};
    if(q(d)<x) d=std::nextafter(d,INFINITY);
    return d;
}
inline std::uint64_t bits(double x) { std::uint64_t b; std::memcpy(&b,&x,sizeof b); return b; }
inline EvalRoundExactScale exact(const mpq_class& x,double runtime) { return {x.get_num().get_str(),x.get_den().get_str(),bits(runtime)}; }
inline BootstrapBound bound(const mpq_class& x,const std::string& why) { return {up(x),BootstrapBoundKind::Deterministic,why,{}}; }
inline bool known(const BootstrapBound& b) { return b.kind==BootstrapBoundKind::Deterministic && std::isfinite(b.upperBound) && b.upperBound>=0 && !b.provenance.empty() && b.failureEventIds.empty(); }
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
    std::string operationLabel;
    std::string keyProvenance;
    Builder(const SealAdapter& a,const Cipher& c,const BootstrapBound& noise,std::string label="EvalRound"):
        primes(a.coeffModulusValues(c)),special(a.specialKeyModulusValue()),
        degree(2*a.slotCount()),initialIndex(a.chainIndex(c)),keyNoise(q(noise.upperBound)),
        operationLabel(std::move(label)),keyProvenance(noise.provenance+"; N="+std::to_string(degree)+"; secret coefficient support=1; P="+std::to_string(special)+"; B_key="+std::to_string(noise.upperBound)) {}
    mpq_class divideRound(std::size_t components,double scale) const {
        return finiteSupportDivideRound(degree,1,components,q(scale));
    }
    mpq_class keySwitch(std::size_t level,double scale) const {
        const std::vector<std::uint64_t> active(primes.begin(),primes.end()-level);
        return finiteSupportKeyNoise(degree,keyNoise,active,special,q(scale))+divideRound(2,scale);
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
        const std::string prefix=operationLabel+" "+n.stage+" node "+std::to_string(i)+": ";
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
    std::size_t input(double scale,double M,double E,std::size_t level=0) {
        EvalRoundExecutionNode n; n.operation=Op::Input; n.stage="input";
        n.chainIndex=initialIndex-level; n.activePrimes.assign(primes.begin(),primes.end()-level);
        n.outputScale=n.arithmeticScale=exact(q(scale),scale);
        const auto i=nodes.size(); nodes.push_back(n); states.push_back({scale,level,2,q(M),q(E)}); rounded.emplace_back(0);
        publish(i,calculate(i,states)); return i;
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
    std::size_t switchLevel(std::size_t a,std::size_t level,const std::string& stage) {
        if(level<states[a].level||level>=primes.size()) throw Failure{Status::InsufficientLevels,"Invalid target modulus level"};
        auto n=nodes[a]; n.operation=Op::ModSwitch; n.inputs={a}; n.stage=stage;
        n.chainIndex=initialIndex-level; n.activePrimes.assign(primes.begin(),primes.end()-level);
        n.inputScales={nodes[a].outputScale}; n.arithmeticScale=n.outputScale;
        auto state=states[a]; state.level=level;
        auto i=nodes.size(); nodes.push_back(n); states.push_back(state); rounded.emplace_back(0);
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

}
