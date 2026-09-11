#include "m2424/experimental/evalmod_analysis/certified_diagonal.hpp"
#include "m2424/canonical_embedding_reference.hpp"
#include <mpfr.h>
#include <cmath>
#include <stdexcept>
namespace m2424::experimental {
namespace {
constexpr mpfr_prec_t precision=192;
struct Number {
    mpfr_t v;
    Number() { mpfr_init2(v,precision); mpfr_set_zero(v,0); }
    Number(const Number& n):Number() { mpfr_set(v,n.v,MPFR_RNDN); }
    Number& operator=(const Number& n) { mpfr_set(v,n.v,MPFR_RNDN); return *this; }
    ~Number() { mpfr_clear(v); }
};
struct C { Number re,im; };
C add(const C& a,const C& b,bool subtract=false) {
    C c; if(subtract) { mpfr_sub(c.re.v,a.re.v,b.re.v,MPFR_RNDN); mpfr_sub(c.im.v,a.im.v,b.im.v,MPFR_RNDN); }
    else { mpfr_add(c.re.v,a.re.v,b.re.v,MPFR_RNDN); mpfr_add(c.im.v,a.im.v,b.im.v,MPFR_RNDN); } return c;
}
C multiply(const C& a,const C& b) {
    C c; Number x,y;
    mpfr_mul(x.v,a.re.v,b.re.v,MPFR_RNDN); mpfr_mul(y.v,a.im.v,b.im.v,MPFR_RNDN); mpfr_sub(c.re.v,x.v,y.v,MPFR_RNDN);
    mpfr_mul(x.v,a.re.v,b.im.v,MPFR_RNDN); mpfr_mul(y.v,a.im.v,b.re.v,MPFR_RNDN); mpfr_add(c.im.v,x.v,y.v,MPFR_RNDN); return c;
}
mpq_class absq(mpq_class x) { return x<0?-x:x; }
double upper(mpq_class x) { double d=x.get_d(); if(mpq_class(d)<x) d=std::nextafter(d,INFINITY); return d; }
}
struct CertifiedRootDiagonalEncoder::Impl {
    std::size_t N; std::vector<C> roots; mpq_class coefficientError;
    explicit Impl(std::size_t n):N(n),roots(2*n) {
        if(mpfr_get_emin()>-8192||mpfr_get_emax()<8192) throw std::invalid_argument("Certified FFT requires a wide MPFR exponent range");
        if(N<4||(N&(N-1))||N>(1u<<20)) throw std::invalid_argument("Unsupported certified FFT degree");
        Number pi,angle; mpfr_const_pi(pi.v,MPFR_RNDN);
        for(std::size_t k=0;k<2*N;++k) {
            mpfr_mul_ui(angle.v,pi.v,k,MPFR_RNDN); mpfr_div_ui(angle.v,angle.v,N,MPFR_RNDN);
            mpfr_sin_cos(roots[k].im.v,roots[k].re.v,angle.v,MPFR_RNDN);
        }
        // u=2^-p. Correctly rounded pi, angle and trig give root error <=64u
        // in complex modulus. Butterfly complex products/additions contribute
        // <=64u(M+E)(1+rootError), with no MPFR exponent under/overflow here.
        mpq_class u=mpq_class(1)/(mpz_class(1)<<precision),w=64*u,E=w,M=1;
        for(std::size_t size=2;size<=N;size*=2) { E=2*E+w*(M+E)+64*u*(M+E)*(1+w); M*=2; }
        coefficientError=(E+w*(M+E)+64*u*(M+E)*(1+w))/N;
    }
};
CertifiedRootDiagonalEncoder::CertifiedRootDiagonalEncoder(std::size_t n):impl_(std::make_unique<Impl>(n)) {}
CertifiedRootDiagonalEncoder::~CertifiedRootDiagonalEncoder()=default;
CertifiedDiagonal CertifiedRootDiagonalEncoder::encode(SealAdapter& adapter,const std::vector<int>& exponents,double scale,const mpq_class& prefactor) const {
    const auto N=impl_->N;
    if(adapter.slotCount()*2!=N||exponents.size()!=N/2||!std::isfinite(scale)||scale<=0) throw std::invalid_argument("Certified diagonal shape/scale mismatch");
    mpq_class gain=prefactor; gain.canonicalize();
    std::vector<C> values(N); auto r=canonicalEmbeddingRootExponents(N);
    for(std::size_t j=0;j<N/2;++j) {
        if(exponents[j]<-1||exponents[j]>=int(2*N)) throw std::invalid_argument("Invalid exact root exponent");
        if(exponents[j]<0) continue;
        auto v=impl_->roots[exponents[j]]; values[(r[j]-1)/2]=v;
        mpfr_neg(v.im.v,v.im.v,MPFR_RNDN); values[(2*N-r[j]-1)/2]=v;
    }
    for(std::size_t i=1,j=0;i<N;++i) { std::size_t bit=N/2; for(;j&bit;bit/=2) j^=bit; j^=bit; if(i<j) std::swap(values[i],values[j]); }
    for(std::size_t size=2;size<=N;size*=2) for(std::size_t base=0;base<N;base+=size) for(std::size_t j=0;j<size/2;++j) {
        auto v=multiply(values[base+j+size/2],impl_->roots[(2*N-2*N*j/size)%(2*N)]),u=values[base+j];
        values[base+j]=add(u,v); values[base+j+size/2]=add(u,v,true);
    }
    auto primes=adapter.dataModulusValues(); primes.push_back(adapter.specialKeyModulusValue());
    std::vector<std::uint64_t> residues(N*primes.size());
    mpq_class error=N*impl_->coefficientError*absq(gain);
    mpz_class Q=1; for(auto p:primes) Q*=mpz_class(std::to_string(p));
    for(std::size_t k=0;k<N;++k) {
        auto v=multiply(values[k],impl_->roots[(2*N-k)%(2*N)]);
        mpq_class midpoint; mpfr_get_q(midpoint.get_mpq_t(),v.re.v); midpoint=midpoint/N*gain;
        const mpq_class scaled=midpoint*mpq_class(scale);
        mpq_class shifted=scaled+mpq_class(1,2); mpz_class encoded;
        mpz_fdiv_q(encoded.get_mpz_t(),shifted.get_num_mpz_t(),shifted.get_den_mpz_t());
        if(absq(mpq_class(encoded))*2>=mpq_class(Q)) throw std::invalid_argument("Encoded diagonal exceeds centered key modulus");
        error+=absq(mpq_class(encoded)/mpq_class(scale)-midpoint);
        for(std::size_t i=0;i<primes.size();++i) {
            mpz_class v=encoded%mpz_class(std::to_string(primes[i])); if(v<0) v+=mpz_class(std::to_string(primes[i]));
            residues[i*N+k]=std::stoull(v.get_str());
        }
    }
    return {adapter.encodePolynomialRnsAtKeyScale(residues,scale),{upper(error),BootstrapBoundKind::Deterministic,
        "Exact-root inverse FFT at MPFR p=192; analytic 64u root/butterfly bounds plus exact coefficient rounding; slot sup error <= sum coefficient errors",{}}};
}
}
