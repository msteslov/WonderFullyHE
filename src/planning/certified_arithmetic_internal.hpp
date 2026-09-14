#pragma once
#include "m2424/evalround_execution.hpp"
#include "m2424/experimental/evalmod_analysis/finite_support_arithmetic.hpp"
#include <gmpxx.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <optional>
#include <stdexcept>
#include <utility>
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
inline mpq_class dyadicUpper(const mpq_class& x,std::size_t precisionBits=384) {
    if(x<0)throw Failure{Status::InvalidInput,"Negative arithmetic bound"};
    const mpz_class denominator=mpz_class(1)<<precisionBits;
    return mpq_class(ceilq(x*denominator),denominator);
}
inline std::optional<double> projectUp(const mpq_class& x) {
    if(x<0) throw Failure{Status::InvalidInput,"Negative arithmetic bound"};
    double d=x.get_d();
    if(!std::isfinite(d)) return std::nullopt;
    if(q(d)<x) d=std::nextafter(d,INFINITY);
    if(!std::isfinite(d)) return std::nullopt;
    return d;
}
inline std::uint64_t bits(double x) { std::uint64_t b; std::memcpy(&b,&x,sizeof b); return b; }
inline EvalRoundExactScale exact(const mpq_class& x,double runtime) { return {x.get_num().get_str(),x.get_den().get_str(),bits(runtime)}; }
inline EvalRoundExactBound exactBound(const mpq_class& x,const std::string& why) {
    if(x<0) throw Failure{Status::InvalidInput,"Negative exact arithmetic bound"};
    return {x.get_num().get_str(),x.get_den().get_str(),
        BootstrapBoundKind::Deterministic,why,projectUp(x)};
}
inline BootstrapBound bound(const mpq_class& x,const std::string& why) {
    const auto projected=projectUp(x);
    return {projected.value_or(std::numeric_limits<double>::infinity()),
        BootstrapBoundKind::Deterministic,why,{}};
}
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
    std::optional<std::string> firstUnprojectableBound;
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
    mpz_class activeModulus(std::size_t level) const {
        mpz_class Q=1;
        for(auto i=primes.begin();i!=primes.end()-level;++i) Q*=integer(*i);
        return Q;
    }
    double baselineDyadicScale(std::size_t level) const {
        const auto minimumPrime=*std::min_element(primes.begin(),primes.end()-level);
        int exponent=0;
        std::frexp(static_cast<double>(minimumPrime),&exponent);
        --exponent;
        while(exponent>=std::numeric_limits<double>::min_exponent-1
              &&q(std::ldexp(1.,exponent))>=mpq_class(integer(minimumPrime))) --exponent;
        if(exponent<std::numeric_limits<double>::min_exponent-1)
            throw Failure{Status::ScaleScheduleInfeasible,"No positive normal dyadic scale below active primes"};
        return std::ldexp(1.,exponent);
    }
    bool scalarScaleFeasible(std::size_t a,const mpq_class& k,double S) const {
        const auto& input=states.at(a);
        const auto Q=activeModulus(input.level);
        const double outputScale=input.scale*S;
        if(!std::isfinite(S)||S<=0||!std::isfinite(outputScale)||outputScale<=0
           ||q(outputScale)>=mpq_class(Q))return false;
        const mpz_class encoded=roundq(k*q(S));
        if(absq(mpq_class(encoded))*2>=mpq_class(Q))return false;
        const mpq_class delta=absq(mpq_class(encoded)/q(S)-k);
        const mpq_class M=input.M*absq(k);
        mpq_class E=input.E*absq(k)+(input.M+input.E)*delta;
        E+=absq(q(input.scale)*q(S)/q(outputScale)-1)*(M+E);
        return Q-2*ceilq(q(outputScale)*(M+E))>0;
    }
    double scheduleScalarScaleForOutput(std::size_t a,const mpq_class& k,double outputScale,
                                        const std::string& stage) const {
        const double S=outputScale/states.at(a).scale;
        if(std::isfinite(S)&&S>0&&states.at(a).scale*S==outputScale
           &&scalarScaleFeasible(a,k,S))return S;
        return scheduleScalarScale(a,k,stage);
    }
    double scheduleScalarScale(std::size_t a,const mpq_class& k,const std::string& stage) const {
        const auto& input=states.at(a);
        const auto Q=activeModulus(input.level);
        int exponent=std::ilogb(baselineDyadicScale(input.level));
        Status lastStatus=Status::ScaleScheduleInfeasible;
        std::string lastWhy="no finite positive candidate";
        constexpr std::size_t maximumCandidates=2048;
        std::size_t attempted=0;
        for(;exponent>=std::numeric_limits<double>::min_exponent-1
              &&attempted<maximumCandidates;--exponent,++attempted) {
            const double S=std::ldexp(1.,exponent);
            const double outputScale=input.scale*S;
            if(!std::isfinite(S)||S<=0||!std::isfinite(outputScale)||outputScale<=0) {
                lastStatus=Status::ScaleScheduleInfeasible;
                lastWhy="non-finite coefficient or output scale";
                continue;
            }
            if(q(outputScale)>=mpq_class(Q)) {
                lastStatus=Status::ScaleScheduleInfeasible;
                lastWhy="output scale is not below the exact active modulus";
                continue;
            }
            const mpz_class encoded=roundq(k*q(S));
            if(absq(mpq_class(encoded))*2>=mpq_class(Q)) {
                lastStatus=Status::HeadroomViolation;
                lastWhy="encoded coefficient integer does not fit the exact active modulus";
                continue;
            }
            const mpq_class delta=absq(mpq_class(encoded)/q(S)-k);
            const mpq_class M=input.M*absq(k);
            mpq_class E=input.E*absq(k)+(input.M+input.E)*delta;
            const mpq_class representation=absq(q(input.scale)*q(S)/q(outputScale)-1)*(M+E);
            E+=representation;
            if(Q-2*ceilq(q(outputScale)*(M+E))<=0) {
                lastStatus=Status::HeadroomViolation;
                lastWhy="non-positive exact centered headroom after coefficient multiplication";
                continue;
            }
            return S;
        }
        throw Failure{lastStatus,stage+": no finite dyadic coefficient scale after bounded search ("
            +std::to_string(attempted)+" candidates); "+lastWhy};
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
        c.state.M=M; c.state.E=E+c.local;
        return c;
    }
    void publish(std::size_t i,const Calculation& c) {
        auto& n=nodes[i]; states[i]=c.state; n.ciphertextComponents=c.state.components;
        const std::string prefix=operationLabel+" "+n.stage+" node "+std::to_string(i)+": ";
        const std::string magnitudeWhy=prefix+"triangle/product norms; v9 digit magnitude at stage boundaries";
        const std::string propagatedWhy=prefix+"semantic error propagation through exact polynomial";
        std::string localWhy=prefix+"exact scalar encoding, finite-support key switching, component divide-round, exact scale ratio";
        if(n.requiredKey!=EvalRoundEvaluationKey::None) localWhy+="; "+keyProvenance;
        if(n.operation==Op::Rescale) localWhy+="; N="+std::to_string(degree)+"; secret coefficient support=1; components="+std::to_string(c.state.components);
        const std::string semanticWhy=prefix+"propagated plus local arithmetic error";
        const std::string encodingWhy=prefix+"exact rational |roundedInteger/encodingScale-constant|";
        const std::string representationWhy=prefix+"exact rational scale ratio against binary64 metadata";
        n.exactIdealMagnitude=exactBound(c.state.M,magnitudeWhy);
        n.exactPropagatedSemanticError=exactBound(c.propagated,propagatedWhy);
        n.exactLocalArithmeticError=exactBound(c.local,localWhy);
        n.exactSemanticError=exactBound(c.state.E,semanticWhy);
        n.exactConstantEncodingError=exactBound(c.encoding,encodingWhy);
        n.exactScaleRepresentationError=exactBound(c.representation,representationWhy);
        if(!firstUnprojectableBound) {
            if(!n.exactIdealMagnitude.outwardBinary64)
                firstUnprojectableBound=prefix+"ideal magnitude";
            else if(!n.exactPropagatedSemanticError.outwardBinary64)
                firstUnprojectableBound=prefix+"propagated semantic error";
            else if(!n.exactLocalArithmeticError.outwardBinary64)
                firstUnprojectableBound=prefix+"local arithmetic error";
            else if(!n.exactSemanticError.outwardBinary64)
                firstUnprojectableBound=prefix+"semantic error";
            else if(!n.exactConstantEncodingError.outwardBinary64)
                firstUnprojectableBound=prefix+"constant encoding error";
            else if(!n.exactScaleRepresentationError.outwardBinary64)
                firstUnprojectableBound=prefix+"scale representation error";
        }
        n.idealMagnitude=bound(c.state.M,magnitudeWhy);
        n.propagatedSemanticError=bound(c.propagated,propagatedWhy);
        n.localArithmeticError=bound(c.local,localWhy);
        n.semanticError=bound(c.state.E,semanticWhy);
        n.constantEncodingError=bound(c.encoding,encodingWhy);
        n.scaleRepresentationError=bound(c.representation,representationWhy);
        mpz_class Q=1; for(auto prime:n.activePrimes) Q*=integer(prime);
        mpz_class margin=Q-2*ceilq(q(c.state.scale)*(c.state.M+c.state.E));
        if(margin<=0) throw Failure{Status::HeadroomViolation,prefix+"centered no-wrap proof unavailable"};
        n.centeredHeadroomNumerator=margin.get_str();
        n.centeredHeadroomProvenance=prefix+"exact Q/2-ceil(runtimeScale*(M+E)); inverse canonical embedding coefficient norm <= slot sup norm";
    }
    std::size_t input(double scale,mpq_class M,mpq_class E,std::size_t level=0) {
        EvalRoundExecutionNode n; n.operation=Op::Input; n.stage="input";
        n.chainIndex=initialIndex-level; n.activePrimes.assign(primes.begin(),primes.end()-level);
        n.outputScale=n.arithmeticScale=exact(q(scale),scale);
        const auto i=nodes.size(); nodes.push_back(n); states.push_back({scale,level,2,std::move(M),std::move(E)}); rounded.emplace_back(0);
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
            const mpq_class represented=mpq_class(encoded)/q(constantScale);
            n.constantEncodingScale=exact(q(constantScale),constantScale);
            n.encodedConstantInteger=encoded.get_str();
            n.representedConstantNumerator=represented.get_num().get_str();
            n.representedConstantDenominator=represented.get_den().get_str();
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
        if(!std::isfinite(s.scale) || s.scale<=0)
            throw Failure{Status::ScaleScheduleInfeasible,n.stage+": output scale overflow/underflow"};
        n.activePrimes.assign(primes.begin(),primes.end()-s.level);
        n.chainIndex=initialIndex-s.level;
        for(auto i:inputs) n.inputScales.push_back(exact(q(states[i].scale),states[i].scale));
        n.outputScale=exact(q(s.scale),s.scale); n.arithmeticScale=exact(arithmetic,s.scale);
        // SEAL requires scale < modulus (also before relinearize/rescale).
        const mpz_class Q=activeModulus(s.level);
        if(q(s.scale)>=mpq_class(Q) || absq(mpq_class(encoded))*2>=mpq_class(Q))
            throw Failure{Status::HeadroomViolation,n.stage+": scale or plaintext exceeds active modulus"};
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
    std::pair<std::size_t,std::size_t> alignForAdd(std::size_t a,std::size_t c,
                                                   const std::string& stage) {
        if(states[a].level<states[c].level)a=alignLevel(a,c,stage);
        if(states[c].level<states[a].level)c=alignLevel(c,a,stage);
        const double sa=states[a].scale,sc=states[c].scale;
        if(bits(sa)!=bits(sc)) {
            const double multiplier=sc*sc==sa?sc:sa/sc;
            if(std::isfinite(multiplier)&&multiplier>0&&sc*multiplier==sa)
                c=scalar(c,1,multiplier,stage);
            else {
                a=scalar(a,1,sc,stage);
                c=scalar(c,1,sa,stage);
            }
        }
        return {a,c};
    }
    std::pair<std::size_t,std::size_t> alignForAddByProduct(std::size_t a,std::size_t c,
                                                            const std::string& stage) {
        if(states[a].level<states[c].level)a=alignLevel(a,c,stage);
        if(states[c].level<states[a].level)c=alignLevel(c,a,stage);
        if(bits(states[a].scale)!=bits(states[c].scale)) {
            const double sa=states[a].scale,sc=states[c].scale;
            a=scalar(a,1,sc,stage);
            c=scalar(c,1,sa,stage);
        }
        return {a,c};
    }
    std::size_t stabilizeScale(std::size_t node,const std::string& stage) {
        double target=baselineDyadicScale(states[node].level);
        while(states[node].scale>=target*target) {
            node=add(Op::Rescale,{node},stage+" scale stabilization");
            target=baselineDyadicScale(states[node].level);
        }
        if(states[node].scale<target) {
            const double ratio=target/states[node].scale;
            const int exponent=static_cast<int>(std::ceil(std::log2(ratio)));
            if(exponent>0)
                node=scalar(node,1,std::ldexp(1.,exponent),stage+" scale stabilization");
        }
        return node;
    }
    std::size_t cleaner(std::size_t a,const std::string& stage,bool stabilizeInput=false) {
        // Baseline: a*a -> relinearize -> rescale, then *(3-2a) ->
        // relinearize -> rescale. No deferred/thrifty tensor evaluation.
        if(stabilizeInput)a=stabilizeScale(a,stage+" input");
        auto square=reduce(mul(a,a,stage),stage);
        auto linear=plus(scalar(a,-2,1,stage),3,stage);
        linear=alignLevel(linear,square,stage);
        return reduce(mul(square,linear,stage),stage);
    }
    mpq_class localBlock(std::size_t first,std::size_t last,const mpq_class& inputMagnitude) const {
        // Cleaner blocks only refer to their input or to an earlier node in the
        // same block.  Do not copy every previously published exact state: K64
        // carries large rational certificates and that quadratic copying cost
        // otherwise dominates compilation.
        std::vector<State> s(last+1);
        s[first]=states[first]; s[first].M=inputMagnitude; s[first].E=0;
        for(auto i=first+1;i<=last;++i) s[i]=calculate(i,s).state;
        return s[last].E;
    }
    void tightenMagnitude(std::size_t node,const mpq_class& magnitude,
                          const std::string& provenance="binary digit recurrence bounds exact polynomial ideal magnitude") {
        states[node].M=magnitude;
        const std::string& why=provenance;
        nodes[node].exactIdealMagnitude=exactBound(magnitude,why);
        nodes[node].idealMagnitude=bound(magnitude,why);
        mpz_class Q=1; for(auto prime:nodes[node].activePrimes) Q*=integer(prime);
        const mpz_class margin=Q-2*ceilq(q(states[node].scale)*(states[node].M+states[node].E));
        if(margin<=0) throw Failure{Status::HeadroomViolation,"Digit boundary headroom"};
        nodes[node].centeredHeadroomNumerator=margin.get_str();
    }
    void compactSemanticError(std::size_t node,const std::string& provenance) {
        const auto compact=dyadicUpper(states[node].E);
        states[node].E=compact;
        nodes[node].exactSemanticError=exactBound(compact,provenance);
        nodes[node].semanticError=bound(compact,provenance);
        const mpz_class margin=activeModulus(states[node].level)
            -2*ceilq(q(states[node].scale)*(states[node].M+compact));
        if(margin<=0)throw Failure{Status::HeadroomViolation,
            nodes[node].stage+": compact exact semantic bound has no centered headroom"};
        nodes[node].centeredHeadroomNumerator=margin.get_str();
        nodes[node].centeredHeadroomProvenance=provenance
            +"; exact Q/2-ceil(runtimeScale*(M+E))";
    }
};

}
