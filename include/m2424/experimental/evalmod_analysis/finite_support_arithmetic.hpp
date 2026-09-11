#pragma once
#include <gmpxx.h>
#include "m2424/bootstrap_contract.hpp"
#include <seal/util/defines.h>
#include <cstdint>
#include <vector>
#include <string>
namespace m2424::experimental {
inline BootstrapBound finiteSupportBackendKeyNoise() {
#ifdef SEAL_USE_GAUSSIAN_NOISE
    return {std::numeric_limits<double>::infinity(),BootstrapBoundKind::Unknown,"Gaussian evaluation-key finite support has not been certified",{}};
#else
    return {21,BootstrapBoundKind::Deterministic,"SEAL sample_poly_cbd: difference of two 21-bit Hamming weights",{}};
#endif
}
inline mpz_class supportInteger(std::uint64_t x) { return mpz_class(std::to_string(x)); }
inline mpq_class finiteSupportDivideRound(std::size_t N,std::size_t secretSupport,std::size_t components,const mpq_class& scale) {
    mpz_class power=1,sum=0;
    for(std::size_t i=0;i<components;++i) { sum+=power; power*=supportInteger(N)*supportInteger(secretSupport); }
    return mpq_class(supportInteger(N)*sum)/2/scale;
}
// Noise alone, before any componentwise special-prime ModDown rounding.
inline mpq_class finiteSupportKeyNoise(std::size_t N,const mpq_class& noiseSupport,const std::vector<std::uint64_t>& primes,std::uint64_t P,const mpq_class& scale) {
    mpz_class sum=0; for(auto prime:primes) sum+=supportInteger(prime-1);
    return mpq_class(supportInteger(N)*supportInteger(N)*sum)*noiseSupport/supportInteger(P)/scale;
}
}
