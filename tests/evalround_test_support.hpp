#pragma once
#include "m2424/experimental/evalmod_analysis/evalround_execution.hpp"
#include <gmpxx.h>
#include <cmath>
#include <stdexcept>
namespace evalround_test {
inline m2424::BootstrapBound inputBound(const m2424::CkksProfile& p) {
    // SEAL asymmetric encryption: nu=e_pk*u+e0+e1*s, ternary s,u,
    // CBD |e_i|<=21. Bound embeddings by N times coefficient support.
    // Encryptor also divides the zero encryption by the preceding prime;
    // ignoring that reduction is conservative, then add its two-component
    // rounding bound N(1+N)/2. Exact scalar encoding below has zero error.
    const mpz_class N(std::to_string(p.polyModulusDegree));
    const mpq_class e=(mpq_class(21)*(2*N*N+N)+mpq_class(N*(1+N))/2)/mpq_class(p.scale);
    double d=e.get_d(); if(mpq_class(d)<e) d=std::nextafter(d,INFINITY);
    return {d,m2424::BootstrapBoundKind::Deterministic,
        "Exact scalar plaintext; SEAL ternary s,u / CBD support 21: [21(2N^2+N)+N(1+N)/2]/scale",{}};
}
inline m2424::Cipher encryptExactScalar(m2424::SealAdapter& a,const m2424::Cipher& context,double z,double scale) {
    const mpq_class scaled=mpq_class(z)*mpq_class(scale);
    if(scaled.get_den()!=1) throw std::invalid_argument("Test scalar is not exactly encodable");
    std::vector<std::uint64_t> residues;
    for(auto prime:a.coeffModulusValues(context)) {
        const mpz_class modulus(std::to_string(prime)); mpz_class r;
        mpz_mod(r.get_mpz_t(),scaled.get_num_mpz_t(),modulus.get_mpz_t()); residues.push_back(std::stoull(r.get_str()));
    }
    return a.encrypt(a.encodeScalarRnsAtScaleFor(residues,scale,context,0));
}
}
