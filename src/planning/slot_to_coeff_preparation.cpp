#include "../core/slot_to_coeff_internal.hpp"
#include "m2424/experimental/evalmod_analysis/certified_diagonal.hpp"
#include "m2424/experimental/evalmod_analysis/finite_support_arithmetic.hpp"
#include <cfenv>
#include <cmath>
#include <cstring>
#include <set>
#include <stdexcept>
namespace m2424 {
namespace {
using Status=BootstrapCertificationStatus;
struct Failure { Status status; std::string why; };
mpq_class q(double x) { return mpq_class(x); }
mpq_class absq(mpq_class x) { return x<0?-x:x; }
mpz_class z(std::uint64_t x) { return mpz_class(std::to_string(x)); }
double up(const mpq_class& x) { double d=x.get_d(); if(!std::isfinite(d)) throw Failure{Status::RequiredBoundUnavailable,"Nonfinite linear-transform bound"}; if(q(d)<x) d=std::nextafter(d,INFINITY); return d; }
BootstrapBound bound(const mpq_class& x,std::string p) { return {up(x),BootstrapBoundKind::Deterministic,std::move(p),{}}; }
bool known(const BootstrapBound& x) { return x.kind==BootstrapBoundKind::Deterministic&&x.failureEventIds.empty()&&std::isfinite(x.upperBound)&&x.upperBound>=0&&!x.provenance.empty(); }
std::uint64_t bits(double x) { std::uint64_t b; std::memcpy(&b,&x,8); return b; }
LinearTransformScale scale(const mpq_class& x,double d) { return {x.get_num().get_str(),x.get_den().get_str(),bits(d)}; }
std::string headroom(const std::vector<std::uint64_t>& primes,const mpq_class& S,const mpq_class& magnitude,const std::string& label) {
    mpz_class Q=1; for(auto p:primes) Q*=z(p);
    const mpq_class scaled=S*magnitude; mpz_class ceiling; mpz_cdiv_q(ceiling.get_mpz_t(),scaled.get_num_mpz_t(),scaled.get_den_mpz_t());
    const mpz_class margin=Q-2*ceiling;
    if(margin<=0) throw Failure{Status::HeadroomViolation,label+": no centered headroom"};
    return margin.get_str();
}
}
PreparedSlotToCoeffPlan SlotToCoeffPlan::prepare(SealAdapter& adapter,const Cipher& x,const Cipher& y,const SlotToCoeffContract& contract) const {
    PreparedSlotToCoeffPlan result; auto p=std::make_shared<PreparedSlotToCoeffPlan::Impl>(); result.impl_=p;
    auto& cert=p->certificate;
    try {
        if(std::fegetround()!=FE_TONEAREST) throw Failure{Status::ScaleScheduleInfeasible,"Round-to-nearest required"};
        const auto N=impl_->degree,S=N/2,depth=impl_->metrics.depth;
        if(adapter.slotCount()!=S) throw Failure{Status::InvalidInput,"Degree/context mismatch"};
        const auto info=adapter.info(x),other=adapter.info(y);
        if(info.ciphertextSize<2||info.ciphertextSize>3||other.ciphertextSize<2||other.ciphertextSize>3) throw Failure{Status::InvalidInput,"Baseline accepts size-two or size-three ciphertexts"};
        p->inputComponents[0]=info.ciphertextSize; p->inputComponents[1]=other.ciphertextSize;
        cert.keys=impl_->keys; cert.keys.relinearization=info.ciphertextSize==3||other.ciphertextSize==3;
        if(bits(info.scale)!=bits(other.scale)) throw Failure{Status::InputScaleMismatch,"Baseline requires equal input scales; no metadata alignment"};
        p->primes=adapter.coeffModulusValues(x); p->chainIndex=info.chainIndex;
        if(p->primes!=adapter.coeffModulusValues(y)) throw Failure{Status::InvalidInput,"Input active primes differ"};
        if(depth>=p->primes.size()) throw Failure{Status::InsufficientLevels,"One data prime per factor and an output prime required"};
        if(!known(contract.evaluationKeyNoiseSupport)||!known(experimental::finiteSupportBackendKeyNoise())||contract.evaluationKeyNoiseSupport.upperBound<experimental::finiteSupportBackendKeyNoise().upperBound)
            throw Failure{Status::RequiredBoundUnavailable,"Evaluation-key support bound unknown or below actual backend support"};
        for(std::size_t b=0;b<2;++b) if(!known(contract.inputMagnitude[b])||!known(contract.inputError[b])) throw Failure{Status::RequiredBoundUnavailable,"Input magnitude/error bound with provenance required"};
        if(cert.keys.relinearization&&!adapter.hasRelinKeys()) throw Failure{Status::MissingEvaluationKeys,"Input relinearization key required before execution"};
        if(!adapter.hasRotationKeys(impl_->keys.rotations)) throw Failure{Status::MissingEvaluationKeys,"Required SlotToCoeff rotations unavailable"};
        p->degree=N; p->factorization=impl_->factorization; p->fingerprint=adapter.contextFingerprint();
        cert.metrics=impl_->metrics;
        const auto special=adapter.specialKeyModulusValue();
        experimental::CertifiedRootDiagonalEncoder encoder(N);
        BootstrapBound finalM[2],finalE[2]; double finalScale=0;
        std::vector<LinearTransformFactorBound> gainFactors(depth+1);
        for(std::size_t b=0;b<2;++b) {
            auto M=contract.inputMagnitude[b],E=contract.inputError[b]; double current=info.scale;
            cert.inputScales[b]=scale(q(current),current);
            headroom(p->primes,q(current),q(M.upperBound)+q(E.upperBound),"input");
            if(p->inputComponents[b]==3) {
                const mpq_class local=experimental::finiteSupportKeyNoise(N,q(contract.evaluationKeyNoiseSupport.upperBound),p->primes,special,q(current))
                    +experimental::finiteSupportDivideRound(N,1,2,q(current));
                E=bound(q(E.upperBound)+local,"Incoming semantic error plus PR-3 finite-support input relinearization");
                SlotToCoeffRuntimeStage stage; stage.operation="input relinearize"; stage.branch=b; stage.chainIndex=info.chainIndex; stage.activePrimes=p->primes;
                stage.inputScale=stage.outputScale=stage.arithmeticScale=scale(q(current),current);
                stage.magnitude=M; stage.semanticError=E;
                stage.localError=bound(local,"Shared PR-3 key noise plus two-component ModDown for one eliminated component");
                stage.scaleRepresentationError=bound(0,"Relinearization preserves scale bits exactly");
                stage.centeredHeadroomNumerator=headroom(p->primes,q(current),q(M.upperBound)+q(E.upperBound),stage.operation);
                stage.headroomProvenance="Exact centered modulus after relinearization";
                cert.inputStages.push_back(std::move(stage)); ++cert.metrics.relinearizations;
            }
            for(std::size_t r=0;r<depth;++r) {
                std::vector<std::uint64_t> active(p->primes.begin(),p->primes.end()-r);
                // Fixed dyadic plaintext scale selected from actual last-prime bit length.
                auto prime=active.back(); int primeBits=0; for(auto v=prime;v;v>>=1) ++primeBits;
                const double T=std::ldexp(1.,primeBits),product=current*T,next=product/static_cast<double>(prime);
                if(!std::isfinite(product)||!std::isfinite(next)||next<=0) throw Failure{Status::ScaleScheduleInfeasible,"Scale overflow"};
                if(b==1 && r+1<depth) p->factors[b].push_back(p->factors[0][r]);
                else {
                    std::map<int,StCPreparedGroup> groups;
                    for(const auto& d:impl_->factors[b][r]) {
                        const auto baby=d.first%impl_->babySteps[r],giant=d.first-baby;
                        std::vector<int> shifted(S,-1); for(std::size_t row=0;row<S;++row) shifted[(row+giant)%S]=d.second[row];
                        auto encoded=encoder.encode(adapter,shifted,T);
                        auto& group=groups[int(giant)]; group.giant=int(giant);
                        group.terms.push_back({int(baby),std::move(encoded.plaintext),encoded.perturbation});
                    }
                    p->factors[b].emplace_back(); for(auto& g:groups) p->factors[b].back().push_back(std::move(g.second));
                }
                std::size_t kappa=0,babyKappa=0;
                for(std::size_t row=0;row<S;++row) {
                    std::size_t count=0,babyCount=0;
                    for(const auto& d:impl_->factors[b][r]) if(d.second[row]>=0) { ++count; babyCount+=d.first%impl_->babySteps[r]!=0; }
                    kappa=std::max(kappa,count); babyKappa=std::max(babyKappa,babyCount);
                }
                mpq_class delta=0,babyDelta=0;
                std::size_t giants=0; const auto& groups=p->factors[b][r];
                const mpq_class babyNoise=experimental::finiteSupportKeyNoise(N,q(contract.evaluationKeyNoiseSupport.upperBound),active,special,q(current));
                const mpq_class giantNoise=experimental::finiteSupportKeyNoise(N,q(contract.evaluationKeyNoiseSupport.upperBound),active,special,q(product));
                const mpq_class rounding=experimental::finiteSupportDivideRound(N,1,2,q(product));
                for(const auto& g:groups) {
                    giants+=g.giant!=0;
                    for(const auto& term:g.terms) {
                        if(!known(term.perturbation)) throw Failure{Status::RequiredBoundUnavailable,"Unknown encoded diagonal perturbation"};
                        delta+=q(term.perturbation.upperBound);
                        if(term.baby) babyDelta+=q(term.perturbation.upperBound);
                    }
                }
                // At each output row, the sum of baby-noise multipliers is
                // bounded by the row support of baby!=0 diagonals plus their
                // encoding perturbation. Giant permutations restore original rows.
                const mpq_class weightedBabyNoise=(babyKappa+babyDelta)*babyNoise;
                SlotToCoeffFactorTrace trace; trace.branch=b; trace.factor=r;
                trace.bounds.kappa=bound(mpq_class(kappa),"Exact symbolic roots: maximum nonzero row count, every entry has modulus one; unique FFT paths");
                trace.bounds.delta=bound(delta,"Sum of certified encoded diagonal sup errors bounds operator infinity norm perturbation");
                const mpq_class m=q(M.upperBound),e=q(E.upperBound),ideal=mpq_class(kappa)*m;
                const mpq_class propagated=mpq_class(kappa)*e+delta*(m+e);
                const mpq_class productRatio=q(current)*q(T)/q(product);
                const mpq_class productRepresentation=absq(productRatio-1)*((kappa+delta)*(m+e)+weightedBabyNoise);
                const mpq_class localDouble=weightedBabyNoise+(groups.size()+1)*rounding+giants*giantNoise+productRepresentation;
                const mpq_class doubleError=propagated+localDouble;
                const mpq_class rescaleRatio=q(product)/z(prime)/q(next);
                const mpq_class rescaleRepresentation=absq(rescaleRatio-1)*(ideal+doubleError);
                const mpq_class rescaleRound=experimental::finiteSupportDivideRound(N,1,2,q(next));
                trace.bounds.localArithmeticError=bound(localDouble+rescaleRepresentation+rescaleRound,
                    "Shared PR-3 finite support: weighted baby key noise + group inner ModDown + giant key noise + final ModDown + rescale + exact dyadic scale ratios; "+contract.evaluationKeyNoiseSupport.provenance+"; N="+std::to_string(N)+"; P="+std::to_string(special));
                trace.propagation=propagateLinearTransform(M,E,trace.bounds);
                if(trace.propagation.result.status!=Status::Certified) throw Failure{Status::RequiredBoundUnavailable,"Linear recurrence failed"};
                // These conservative envelopes hold for every internal baby and group,
                // including pre-ModDown QP accumulators (divide both QP and scale by P).
                auto record=[&](const std::string& operation,const mpq_class& s,const mpq_class& norm) {
                    trace.internalHeadroomProofs.push_back(operation+": twice margin="+headroom(active,s,norm,operation)+"; QP stages use the equivalent Q/P-normalized no-wrap inequality");
                };
                record("baby extension",q(current),m+e+babyNoise);
                const mpq_class groupEnvelope=(kappa+delta)*(m+e)+weightedBabyNoise;
                record("each extended group accumulator",q(current)*q(T),groupEnvelope);
                const mpq_class innerEnvelope=productRatio*groupEnvelope+rounding;
                record("each inner ModDown",q(product),innerEnvelope);
                record("each giant extension",q(product),innerEnvelope+giantNoise);
                record("outer extended accumulator",q(product),groups.size()*innerEnvelope+giants*giantNoise);
                auto stage=[&](std::string operation,std::vector<std::uint64_t> primes,double before,double after,const mpq_class& arithmetic,const mpq_class& error,const mpq_class& local,const mpq_class& representation) {
                    SlotToCoeffRuntimeStage n; n.operation=std::move(operation); n.branch=b; n.factor=r; n.chainIndex=primes.size()-1; n.activePrimes=std::move(primes);
                    n.inputScale=scale(q(before),before); n.outputScale=scale(q(after),after); n.arithmeticScale=scale(arithmetic,after);
                    n.magnitude=bound(ideal,"kappa*M for the exact ideal factor"); n.semanticError=bound(error,"Separate operator perturbation and finite-support runtime errors");
                    n.localError=bound(local,"Operation arithmetic only, excluding diagonal encoding"); n.scaleRepresentationError=bound(representation,"Exact arithmetic scale / binary64 runtime scale ratio");
                    n.centeredHeadroomNumerator=headroom(n.activePrimes,q(after),q(n.magnitude.upperBound)+q(n.semanticError.upperBound),n.operation);
                    n.headroomProvenance="Exact Q - 2*ceil(scale*(M+E)); positive means centered coefficient headroom";
                    trace.runtime.push_back(std::move(n));
                };
                stage("double-hoisted BSGS",active,current,product,q(current)*q(T),doubleError,localDouble,productRepresentation);
                active.pop_back();
                stage("rescale",active,product,next,q(product)/z(prime),q(trace.propagation.semanticError.upperBound),rescaleRound+rescaleRepresentation,rescaleRepresentation);
                if(b==0) gainFactors[r]=trace.bounds;
                else if(trace.bounds.delta.upperBound>gainFactors[r].delta.upperBound) gainFactors[r].delta=trace.bounds.delta;
                M=trace.propagation.idealMagnitude; E=trace.propagation.semanticError;
                cert.factors.push_back(std::move(trace)); current=next;
            }
            finalM[b]=M; finalE[b]=E; finalScale=current;
        }
        // Pair sup norm -> sum has kappa=2 and no encoding/runtime error.
        gainFactors.back()={bound(2,"Exact addition on the pair sup norm"),bound(0,"No encoded coefficient in pair addition"),bound(0,"Exact RNS addition at identical scales")};
        cert.gamma=linearTransformGain(gainFactors);
        if(!known(cert.gamma)) throw Failure{Status::RequiredBoundUnavailable,"Prepared gain unavailable"};
        cert.outputError=bound(q(finalE[0].upperBound)+q(finalE[1].upperBound),"Sum of the two certified branch errors; exact same-scale addition");
        cert.outputScale=scale(q(finalScale),finalScale);
        SlotToCoeffFactorTrace combine; combine.branch=2; combine.factor=depth; combine.bounds=gainFactors.back();
        auto M=bound(std::max(finalM[0].upperBound,finalM[1].upperBound),"Pair maximum magnitude"),E=bound(std::max(finalE[0].upperBound,finalE[1].upperBound),"Pair maximum error");
        combine.propagation=propagateLinearTransform(M,E,combine.bounds);
        SlotToCoeffRuntimeStage sum; sum.operation="pair addition"; sum.branch=2; sum.factor=depth;
        sum.activePrimes.assign(p->primes.begin(),p->primes.end()-depth); sum.chainIndex=sum.activePrimes.size()-1;
        sum.inputScale=sum.outputScale=sum.arithmeticScale=cert.outputScale;
        sum.magnitude=bound(q(finalM[0].upperBound)+q(finalM[1].upperBound),"Sum of branch magnitudes"); sum.semanticError=cert.outputError;
        sum.localError=sum.scaleRepresentationError=bound(0,"Exact same-scale RNS addition");
        sum.centeredHeadroomNumerator=headroom(sum.activePrimes,q(finalScale),q(sum.magnitude.upperBound)+q(cert.outputError.upperBound),"pair addition"); sum.headroomProvenance="Exact centered headroom after branch addition";
        combine.runtime.push_back(sum); cert.factors.push_back(std::move(combine));
        cert.result={Status::Certified,"SlotToCoeff","Standalone radix-2 double-hoisted baseline with prepared-operator and runtime certificates"};
    } catch(const Failure& f) { cert.result={f.status,"SlotToCoeff",f.why}; }
      catch(const std::exception& e) { cert.result={Status::RequiredBoundUnavailable,"SlotToCoeff",e.what()}; }
    return result;
}
}
