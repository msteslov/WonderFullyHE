#include "../core/slot_to_coeff_internal.hpp"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <cfenv>
#include <numeric>
#include <set>
#include <stdexcept>
namespace m2424 {
namespace {
using Map=RootDiagonalMap;
std::uint64_t bits(double v) { std::uint64_t b; std::memcpy(&b,&v,8); return b; }
void entry(Map& m,std::size_t S,std::size_t row,std::size_t col,int e) {
    auto [it,created]=m.emplace((col+S-row)%S,std::vector<int>{});
    if(created) it->second.assign(S,-1);
    if(it->second[row]!=-1) throw std::logic_error("Nonunique FFT path");
    it->second[row]=e;
}
Map compose(const Map& outer,const Map& inner,std::size_t S) {
    Map result;
    for(const auto& a:outer) for(const auto& b:inner) for(std::size_t r=0;r<S;++r)
        if(a.second[r]>=0 && b.second[(r+a.first)%S]>=0)
            entry(result,S,r,(r+a.first+b.first)%S,(a.second[r]+b.second[(r+a.first)%S])%(4*S));
    return result;
}
ComplexVector applyFactor(const Map& m,const ComplexVector& x) {
    ComplexVector y(x.size()); const double pi=std::acos(-1.);
    for(const auto& d:m) for(std::size_t r=0;r<x.size();++r) if(d.second[r]>=0)
        y[r]+=std::polar(1.,2*pi*d.second[r]/(4*x.size()))*x[(r+d.first)%x.size()];
    return y;
}
std::set<int> rotations(const Map& m,std::size_t baby) {
    std::set<int> out;
    for(const auto& d:m) { auto b=d.first%baby,g=d.first-b; if(b) out.insert(int(b)); if(g) out.insert(int(g)); }
    return out;
}
BootstrapContractResult fail(BootstrapCertificationStatus s,const char* p) { return {s,"SlotToCoeff",p}; }
}
SlotToCoeffPlan::SlotToCoeffPlan(std::size_t N,std::size_t depth):SlotToCoeffPlan(N,SlotToCoeffFactorization{}) {
    if(!depth) throw std::invalid_argument("SlotToCoeff depth must be positive");
    std::size_t raw=0; for(std::size_t S=N/2;S>1;S/=2) ++raw;
    raw+=raw/2; depth=std::min(depth,raw);
    std::vector<std::size_t> r(depth,raw/depth); for(std::size_t i=0;i<raw%depth;++i) ++r[i];
    *this=SlotToCoeffPlan(N,SlotToCoeffFactorization{r});
}
SlotToCoeffPlan::SlotToCoeffPlan(std::size_t N,SlotToCoeffFactorization f) {
    if(N<4||(N&(N-1))) throw std::invalid_argument("SlotToCoeff degree must be a power of two >=4");
    auto p=std::make_shared<Impl>(); p->degree=N;
    const auto S=N/2; std::size_t log=0; for(auto s=S;s>1;s/=2) ++log;
    std::vector<Map> raw;
    // Reverse adjoints of the unnormalized canonical inverse FFT: bit swaps
    // first, then conjugate-transposed butterflies from block 2 to S.
    for(std::size_t count=log/2;count>0;--count) {
        auto low=count-1,high=log-1-low; Map m;
        for(std::size_t r=0;r<S;++r) {
            auto c=r; if(((r>>low)&1)!=((r>>high)&1)) c^=(std::size_t(1)<<low)|(std::size_t(1)<<high);
            entry(m,S,r,c,0);
        } raw.push_back(std::move(m));
    }
    for(std::size_t block=2;block<=S;block*=2) {
        Map m; auto half=block/2;
        for(std::size_t base=0;base<S;base+=block) {
            if(block==2) {
                entry(m,S,base,base,0); entry(m,S,base+1,base,0);
                entry(m,S,base,base+1,int(S/2)); entry(m,S,base+1,base+1,int(3*S/2));
            } else {
                std::size_t power=1;
                for(std::size_t j=0;j<half;++j) {
                    int e=int(power*S/block);
                    entry(m,S,base+j,base+j,0); entry(m,S,base+half+j,base+j,0);
                    entry(m,S,base+j,base+half+j,e); entry(m,S,base+half+j,base+half+j,(e+2*S)%(4*S));
                    power=(power*3)%(4*block);
                }
            }
        } raw.push_back(std::move(m));
    }
    if(f.radices.empty()) f.radices.assign(raw.size(),1);
    if(std::any_of(f.radices.begin(),f.radices.end(),[](auto r){return !r;}) || std::accumulate(f.radices.begin(),f.radices.end(),std::size_t{})!=raw.size()) throw std::invalid_argument("Invalid SlotToCoeff factorization");
    p->factorization=f; std::size_t cursor=0;
    for(auto count:f.radices) { auto m=raw[cursor++]; for(std::size_t k=1;k<count;++k) m=compose(raw[cursor++],m,S); p->factors[0].push_back(std::move(m)); }
    p->factors[1]=p->factors[0]; const auto powers=canonicalEmbeddingRootExponents(N);
    for(auto& d:p->factors[1].back()) for(std::size_t r=0;r<S;++r) if(d.second[r]>=0) d.second[r]=(d.second[r]+powers[r]*S)%(4*S);
    std::set<int> keys;
    for(const auto& m:p->factors[0]) {
        std::size_t best=1,cost=m.size()+1;
        for(std::size_t step=1;step<=S;step*=2) { const auto c=rotations(m,step).size(); if(c<cost) { cost=c; best=step; } }
        p->babySteps.push_back(best); auto rs=rotations(m,best); keys.insert(rs.begin(),rs.end());
        std::set<std::size_t> groups; for(const auto& d:m) groups.insert(d.first-d.first%best);
        p->metrics.rotations+=2*rs.size(); p->metrics.innerModDowns+=2*groups.size(); p->metrics.plaintexts+=2*m.size();
    }
    p->keys.rotations.assign(keys.begin(),keys.end()); p->metrics.depth=f.radices.size();
    p->metrics.rescales=p->metrics.finalModDowns=2*f.radices.size(); impl_=p;
}
std::size_t SlotToCoeffPlan::polyModulusDegree() const { return impl_->degree; }
const SlotToCoeffFactorization& SlotToCoeffPlan::factorization() const { return impl_->factorization; }
SlotToCoeffRequirements SlotToCoeffPlan::requirements() const { return impl_->keys; }
SlotToCoeffMetrics SlotToCoeffPlan::metrics() const { return impl_->metrics; }
std::vector<ComplexVector> SlotToCoeffPlan::applyPlainTrace(const ComplexVector& x,std::size_t branch) const {
    if(x.size()!=impl_->degree/2||branch>1) throw std::invalid_argument("SlotToCoeff plaintext shape mismatch");
    std::vector<ComplexVector> out; auto y=x;
    for(const auto& m:impl_->factors[branch]) { y=applyFactor(m,y); out.push_back(y); } return out;
}
ComplexVector SlotToCoeffPlan::applyPlain(const ComplexVector& x,const ComplexVector& y) const {
    auto a=applyPlainTrace(x,0).back(),b=applyPlainTrace(y,1).back(); for(std::size_t i=0;i<a.size();++i) a[i]+=b[i]; return a;
}
PreparedSlotToCoeffPlan::PreparedSlotToCoeffPlan():impl_(std::make_shared<Impl>()) {}
const SlotToCoeffCertificate& PreparedSlotToCoeffPlan::certificate() const { return impl_->certificate; }
BootstrapContractResult SlotToCoeffPlan::preflight(const SealAdapter& a,const Cipher& x,const Cipher& y,const PreparedSlotToCoeffPlan& p) const {
    const auto& c=p.certificate(); if(c.result.status!=BootstrapCertificationStatus::Certified) return c.result;
    if(std::fegetround()!=FE_TONEAREST) return fail(BootstrapCertificationStatus::ScaleScheduleInfeasible,"Binary64 round-to-nearest required");
    if(impl_->degree!=p.impl_->degree||impl_->factorization.radices!=p.impl_->factorization.radices || a.contextFingerprint()!=p.impl_->fingerprint) return fail(BootstrapCertificationStatus::InvalidInput,"Prepared plan context/factorization mismatch");
    const Cipher* inputs[]={&x,&y};
    for(std::size_t b=0;b<2;++b) {
        auto info=a.info(*inputs[b]);
        if(info.ciphertextSize!=p.impl_->inputComponents[b]) return fail(BootstrapCertificationStatus::InvalidInput,"Prepared input component count mismatch");
        if(bits(info.scale)!=c.inputScales[b].binary64Bits) return fail(BootstrapCertificationStatus::InputScaleMismatch,"Prepared input scale bits mismatch");
        if(info.chainIndex!=p.impl_->chainIndex||a.coeffModulusValues(*inputs[b])!=p.impl_->primes) return fail(BootstrapCertificationStatus::InsufficientLevels,"Prepared active modulus mismatch");
    }
    if(c.keys.relinearization&&!a.hasRelinKeys()) return fail(BootstrapCertificationStatus::MissingEvaluationKeys,"Missing input relinearization keys");
    if(!a.hasRotationKeys(c.keys.rotations)) return fail(BootstrapCertificationStatus::MissingEvaluationKeys,"Missing SlotToCoeff rotations");
    return c.result;
}
SlotToCoeffResult SlotToCoeffPlan::apply(SealAdapter& a,const Cipher& x,const Cipher& y,const PreparedSlotToCoeffPlan& p,const std::function<void(const SlotToCoeffRuntimeStage&,const Cipher&)>& observer) const {
    auto gate=preflight(a,x,y,p); if(gate.status!=BootstrapCertificationStatus::Certified) throw std::invalid_argument(gate.provenance);
    Cipher out[2]={x,y};
    auto audit=[&](const SlotToCoeffRuntimeStage& s,const Cipher& ct) {
        const auto info=a.info(ct);
        if(bits(info.scale)!=s.outputScale.binary64Bits||info.chainIndex!=s.chainIndex||a.coeffModulusValues(ct)!=s.activePrimes) throw std::runtime_error("SlotToCoeff runtime schedule mismatch");
        if(observer) observer(s,ct);
    };
    for(const auto& stage:p.certificate().inputStages) { out[stage.branch]=a.relinearize(out[stage.branch]); audit(stage,out[stage.branch]); }
    for(std::size_t b=0;b<2;++b) for(std::size_t r=0;r<impl_->metrics.depth;++r) {
        std::vector<HoistedBsgsGroup> groups;
        for(const auto& g:p.impl_->factors[b][r]) { HoistedBsgsGroup h; h.giantRotation=g.giant; for(const auto& t:g.terms) h.terms.push_back({t.baby,&t.diagonal}); groups.push_back(std::move(h)); }
        out[b]=a.applyBsgsDoubleHoisted(out[b],groups);
        const auto& trace=p.certificate().factors[b*impl_->metrics.depth+r]; audit(trace.runtime[0],out[b]);
        out[b]=a.rescaleToNext(out[b]); audit(trace.runtime[1],out[b]);
    }
    auto result=a.add(out[0],out[1]); audit(p.certificate().factors.back().runtime[0],result);
    return {std::move(result),p.certificate()};
}
}
