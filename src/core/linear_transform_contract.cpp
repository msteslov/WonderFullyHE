#include "m2424/linear_transform_contract.hpp"
#include <cmath>
#include <cfenv>
namespace m2424 {
namespace {
bool known(const BootstrapBound& x) { return x.kind==BootstrapBoundKind::Deterministic && x.failureEventIds.empty() && std::isfinite(x.upperBound) && x.upperBound>=0 && !x.provenance.empty(); }
double add(double x,double y) { if(!x) return y; if(!y) return x; return std::nextafter(x+y,INFINITY); }
double mul(double x,double y) { if(!x||!y) return 0; return std::nextafter(x*y,INFINITY); }
BootstrapBound b(double x,const char* why) { return {x,BootstrapBoundKind::Deterministic,why,{}}; }
}
LinearTransformPropagation propagateLinearTransform(const BootstrapBound& M,const BootstrapBound& E,const LinearTransformFactorBound& f) {
    LinearTransformPropagation out;
    out.result={BootstrapCertificationStatus::RequiredBoundUnavailable,"LinearTransform","Known deterministic magnitude, error, kappa, delta and local bounds required"};
    if(std::fegetround()!=FE_TONEAREST || !known(M)||!known(E)||!known(f.kappa)||!known(f.delta)||!known(f.localArithmeticError)) return out;
    const double m=mul(f.kappa.upperBound,M.upperBound),p=mul(f.kappa.upperBound,E.upperBound);
    const double d=mul(f.delta.upperBound,add(M.upperBound,E.upperBound));
    const double e=add(add(p,d),f.localArithmeticError.upperBound);
    if(!std::isfinite(m)||!std::isfinite(e)) return out;
    out.idealMagnitude=b(m,"v9: M_next <= kappa*M");
    out.propagatedError=b(p,"v9: kappa*incoming semantic error");
    out.operatorPerturbationError=b(d,"v9: delta*(M+E); encoded operator perturbation only");
    out.semanticError=b(e,"v9: kappa*E + delta*(M+E) + B_local");
    out.result={BootstrapCertificationStatus::Certified,"LinearTransform","Outward-rounded nonnegative v9 recurrence"};
    return out;
}
BootstrapBound linearTransformGain(const std::vector<LinearTransformFactorBound>& factors) {
    if(factors.empty() || std::fegetround()!=FE_TONEAREST) return {};
    double gain=1;
    for(const auto& f:factors) {
        if(!known(f.kappa)||!known(f.delta)) return {};
        gain=mul(gain,add(f.kappa.upperBound,f.delta.upperBound));
    }
    if(!std::isfinite(gain)||gain<=0) return {};
    return b(gain,"v9: product of (kappa+delta) of prepared factors, including pair addition; excludes additive runtime errors");
}
}
