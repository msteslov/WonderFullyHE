#pragma once
#include "m2424/seal_adapter.hpp"
namespace m2424::test {
// Compiled only with BUILD_TESTING. Deliberately insecure tiny-N context and
// synthetic c=(m,1), not encryption. No secret access or lift measurements.
class BootstrapFixture {
public:
    static SealAdapter create(const CkksProfile&);
    static Cipher input(SealAdapter&,const Plain&);
};
}
