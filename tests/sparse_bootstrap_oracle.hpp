#pragma once
#include "m2424/seal_adapter.hpp"
namespace m2424::test {
// Test-only raw secret/semantic access, never linked into non-testing builds.
class SparseBootstrapOracle {
public:
    static std::vector<int> secret(const SealAdapter&);
    static std::vector<double> source(SealAdapter&,const SparseCipher&);
    static std::vector<double> raised(SealAdapter&,const SparseRaisedCipher&);
    static std::vector<double> original(SealAdapter&,const Cipher&);
    static void removeKey(SealAdapter&,bool restoration);
};
}
