#include "../core/coeff_to_slot_cert_internal.hpp"
#include <algorithm>
#include <cmath>
#include <cfenv>
#include <cstring>
#include <stdexcept>
namespace m2424 {
namespace {
std::uint64_t bits(double x) { std::uint64_t v; std::memcpy(&v,&x,8); return v; }
}
const CoeffToSlotBranchTrace& CertifiedEvalRoundPlusCoeffToSlot::hp() const { if(!impl_) throw std::invalid_argument("Empty CtS certificate"); return impl_->traces[0]; }
const CoeffToSlotBranchTrace& CertifiedEvalRoundPlusCoeffToSlot::lp() const { if(!impl_) throw std::invalid_argument("Empty CtS certificate"); return impl_->traces[1]; }
const CoeffToSlotDomainCertificate& CertifiedEvalRoundPlusCoeffToSlot::domain() const { if(!impl_) throw std::invalid_argument("Empty CtS certificate"); return impl_->domain; }
RootLinearTransformPlan EvalRoundPlusCoeffToSlot::certificationLayout(const CoeffToSlotPrefactor& scalar) const {
    auto radices=plan_.factorization().radices; std::reverse(radices.begin(),radices.end());
    // Reverse adjoints of the existing exact forward factors, preserving the
    // original CtS grouping. No floating-point recovery of root exponents.
    SlotToCoeffPlan forward(plan_.polyModulusDegree(),SlotToCoeffFactorization{radices});
    RootLinearTransformPlan out; out.degree=plan_.polyModulusDegree();
    out.factorization.radices=plan_.factorization().radices; out.babySteps=certificationBabySteps();
    out.prefactor=scalar; out.inverse=true; out.combine=false;
    const auto S=out.degree/2;
    for(auto it=forward.impl_->factors[0].rbegin();it!=forward.impl_->factors[0].rend();++it) {
        RootDiagonalMap adjoint;
        for(const auto& d:*it) {
            auto& values=adjoint[(S-d.first)%S]; values.assign(S,-1);
            for(std::size_t r=0;r<S;++r) if(d.second[r]>=0) values[(r+d.first)%S]=(4*S-d.second[r])%(4*S);
        }
        out.factors[0].push_back(std::move(adjoint));
    }
    out.factors[1]=out.factors[0]; const auto powers=canonicalEmbeddingRootExponents(out.degree);
    for(auto& d:out.factors[1].front()) for(std::size_t r=0;r<S;++r) if(d.second[r]>=0)
        d.second[r]=(d.second[r]+4*S-(powers[(r+d.first)%S]%4)*S)%(4*S);
    out.keys.rotations=requirements().rotationSteps;
    const auto m=plan_.metrics(); out.metrics.depth=plan_.depth(); out.metrics.rotations=m.rotationsPerApply;
    out.metrics.rescales=m.rescalesPerApply; out.metrics.innerModDowns=m.innerModDownsPerApply;
    out.metrics.finalModDowns=m.finalModDownsPerApply; out.metrics.plaintexts=m.plaintextMultiplicationsPerApply;
    return out;
}
std::vector<ComplexVector> EvalRoundPlusCoeffToSlot::applyPlainTrace(const ComplexVector& input,std::size_t half,const CoeffToSlotPrefactor& scalar) const {
    auto layout=certificationLayout(scalar); const auto S=layout.degree/2;
    if(input.size()!=S||half>1) throw std::invalid_argument("CtS diagnostic shape mismatch");
    auto x=input; std::vector<ComplexVector> result;
    for(std::size_t r=0;r<layout.factors[half].size();++r) {
        ComplexVector y(S);
        for(const auto& d:layout.factors[half][r]) for(std::size_t row=0;row<S;++row) if(d.second[row]>=0)
            y[row]+=std::polar(1.,2*std::acos(-1.)*d.second[row]/(4*S))*x[(row+d.first)%S];
        if(r==0) for(auto& v:y) v={scalar.multiplyRounded(v.real()/layout.degree),scalar.multiplyRounded(v.imag()/layout.degree)};
        result.push_back(y); x=std::move(y);
    }
    for(auto& v:x) v+=std::conj(v); result.push_back(std::move(x)); return result;
}
BootstrapContractResult EvalRoundPlusCoeffToSlot::preflight(const SealAdapter& a,const RaisedCipher& x,const BootstrapInputContext& source,const CertifiedEvalRoundPlusCoeffToSlot& p) const {
    using Status=BootstrapCertificationStatus;
    if(!p.impl_) return {Status::InvalidInput,"CoeffToSlot","Empty prepared certificate"};
    if(a.info(x).ciphertextSize!=2) return {Status::InvalidInput,"CoeffToSlot","Prepared execution requires size-two raised input"};
    const auto& c=*p.impl_;
    auto gate=preflight(a,x,source,c.hpContract,c.lpContract); if(gate.status!=Status::Certified) return gate;
    if(std::fegetround()!=FE_TONEAREST) return {Status::ScaleScheduleInfeasible,"CoeffToSlot","Round-to-nearest required"};
    if(source.sourcePrimes!=c.source.sourcePrimes||source.scaleBinary64Bits!=c.source.scaleBinary64Bits||source.contextFingerprint!=c.source.contextFingerprint||source.raisedPrimes!=c.source.raisedPrimes||source.specialPrime!=c.source.specialPrime||source.chainIndex!=c.source.chainIndex||plan_.factorization().radices!=c.factorization.radices)
        return {Status::InvalidInput,"CoeffToSlot","Prepared source/scale/factorization mismatch"};
    for(const auto& trace:c.traces) if(trace.certificate.status!=Status::Certified) return trace.certificate;
    return {Status::Certified,"CoeffToSlot","Prepared HP/LP factors, exact source and evaluation keys verified"};
}
EvalRoundPlusCoeffToSlotResult EvalRoundPlusCoeffToSlot::apply(SealAdapter& a,const RaisedCipher& x,const CertifiedEvalRoundPlusCoeffToSlot& p,const std::function<void(BootstrapGate,const SlotToCoeffRuntimeStage&,const Cipher&)>& observer) const {
    if(!p.impl_) throw std::invalid_argument("Empty CtS certificate");
    auto gate=preflight(a,x,p.impl_->source,p); if(gate.status!=BootstrapCertificationStatus::Certified) throw std::invalid_argument(gate.provenance);
    Cipher outputs[4]; const auto depth=plan_.depth();
    for(std::size_t b=0;b<2;++b) {
        const auto& prepared=*p.impl_->branches[b]; const auto& trace=p.impl_->traces[b];
        auto audit=[&](const SlotToCoeffRuntimeStage& s,const Cipher& ct) {
            const auto info=a.info(ct);
            if(info.chainIndex!=s.chainIndex||bits(info.scale)!=s.outputScale.binary64Bits||a.coeffModulusValues(ct)!=s.activePrimes) throw std::runtime_error("CtS certified schedule mismatch");
            if(observer) observer(trace.gate,s,ct);
        };
        for(std::size_t h=0;h<2;++h) {
            auto ct=x.cipher_;
            for(std::size_t r=0;r<depth;++r) {
                std::vector<HoistedBsgsGroup> groups;
                for(const auto& g:prepared.factors[h][r]) { HoistedBsgsGroup group; group.giantRotation=g.giant; for(const auto& t:g.terms) group.terms.push_back({t.baby,&t.diagonal}); groups.push_back(std::move(group)); }
                const auto& f=trace.certifiedFactors[h*(depth+1)+r];
                ct=a.applyBsgsDoubleHoisted(ct,groups); audit(f.runtime[0],ct);
                ct=a.rescaleToNext(ct); audit(f.runtime[1],ct);
            }
            const auto& projection=trace.certifiedFactors[h*(depth+1)+depth];
            auto conjugate=a.conjugate(ct); audit(projection.runtime[0],conjugate);
            ct=a.add(ct,conjugate); audit(projection.runtime[1],ct); outputs[2*b+h]=std::move(ct);
        }
    }
    return {std::move(outputs[0]),std::move(outputs[1]),std::move(outputs[2]),std::move(outputs[3]),p.hp(),p.lp()};
}
}
