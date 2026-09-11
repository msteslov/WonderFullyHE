#pragma once
#include "m2424/bootstrap_contract.hpp"
#include <gmpxx.h>
#include <memory>
namespace m2424::experimental {
struct CertifiedDiagonal { Plain plaintext; BootstrapBound perturbation; };
// Reusable exact-root diagonal encoder. Entry e means exp(2*pi*i*e/(2N));
// -1 means zero. Rational prefactor supports future normalized CtS factors.
class CertifiedRootDiagonalEncoder {
public:
    explicit CertifiedRootDiagonalEncoder(std::size_t degree);
    ~CertifiedRootDiagonalEncoder();
    CertifiedDiagonal encode(SealAdapter&,const std::vector<int>& exponents,double scale,const mpq_class& prefactor=1) const;
private:
    struct Impl; std::unique_ptr<Impl> impl_;
};
}
