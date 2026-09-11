#pragma once
#include "m2424/experimental/evalmod_analysis/finite_support_arithmetic.hpp"
#include "m2424/linear_transform_contract.hpp"
#include <cmath>
#include <cstring>
namespace m2424::linear_certificate {

using Status=BootstrapCertificationStatus;
struct Failure { Status status; std::string why; };
inline mpq_class q(double x) { return mpq_class(x); }
inline mpq_class absq(mpq_class x) { return x<0?-x:x; }
inline mpz_class z(std::uint64_t x) { return mpz_class(std::to_string(x)); }
inline double up(const mpq_class& x) { double d=x.get_d(); if(!std::isfinite(d)) throw Failure{Status::RequiredBoundUnavailable,"Nonfinite linear-transform bound"}; if(q(d)<x) d=std::nextafter(d,INFINITY); return d; }
inline BootstrapBound bound(const mpq_class& x,std::string p) { return {up(x),BootstrapBoundKind::Deterministic,std::move(p),{}}; }
inline bool known(const BootstrapBound& x) { return x.kind==BootstrapBoundKind::Deterministic&&x.failureEventIds.empty()&&std::isfinite(x.upperBound)&&x.upperBound>=0&&!x.provenance.empty(); }
inline std::uint64_t bits(double x) { std::uint64_t b; std::memcpy(&b,&x,8); return b; }
inline LinearTransformScale scale(const mpq_class& x,double d) { return {x.get_num().get_str(),x.get_den().get_str(),bits(d)}; }
inline std::string headroom(const std::vector<std::uint64_t>& primes,const mpq_class& S,const mpq_class& magnitude,const std::string& label) {
    mpz_class Q=1; for(auto p:primes) Q*=z(p);
    const mpq_class scaled=S*magnitude; mpz_class ceiling; mpz_cdiv_q(ceiling.get_mpz_t(),scaled.get_num_mpz_t(),scaled.get_den_mpz_t());
    const mpz_class margin=Q-2*ceiling;
    if(margin<=0) throw Failure{Status::HeadroomViolation,label+": no centered headroom"};
    return margin.get_str();
}

}
