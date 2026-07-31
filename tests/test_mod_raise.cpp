#include "m2424/canonical_embedding_reference.hpp"
#include "m2424/seal_adapter.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <vector>

int main() {
    const m2424::CkksProfile profile{8192, {60, 49, 49, 60}, std::exp2(49.0), 0};
    auto adapter = m2424::SealAdapter::create(profile);
    adapter.generateKeys(false, false);
    const auto input = adapter.encrypt(adapter.encode({0.25, -0.5, 0.75}));
    const auto lowered = adapter.rescaleToNext(adapter.multiplyPlain(
        input, adapter.encodeScalarAtScaleFor(1.0, std::exp2(49.0), input)));
    const auto lowest = adapter.rescaleToNext(adapter.multiplyPlain(
        lowered, adapter.encodeScalarAtScaleFor(1.0, std::exp2(49.0), lowered)));
    const auto raised = adapter.modRaiseToTop(lowered);
    const auto raisedLowest = adapter.modRaiseToTop(lowest);
    const auto sourceCoefficients = adapter.decryptRaisedCoefficientsAtSourceModulus(raised);
    const auto oracleSlots = m2424::coeffToSlotReference(sourceCoefficients);
    const auto inputInfo = adapter.info(lowered);
    const auto lowestInfo = adapter.info(lowest);
    const auto resultInfo = adapter.info(raised);
    const auto lowestResultInfo = adapter.info(raisedLowest);

    bool rejectsTopLevel = false;
    try {
        (void)adapter.modRaiseToTop(input);
    } catch (const std::invalid_argument&) {
        rejectsTopLevel = true;
    }
    const bool ok = resultInfo.chainIndex > inputInfo.chainIndex
        && resultInfo.scale == inputInfo.scale
        && resultInfo.ciphertextSize == inputInfo.ciphertextSize
        && lowestResultInfo.chainIndex > lowestInfo.chainIndex
        && lowestResultInfo.scale == lowestInfo.scale
        && lowestResultInfo.ciphertextSize == lowestInfo.ciphertextSize
        && sourceCoefficients.size() == profile.polyModulusDegree
        && oracleSlots.size() == profile.polyModulusDegree / 2
        && std::all_of(sourceCoefficients.begin(), sourceCoefficients.end(), [](double value) {
            return std::isfinite(value);
        })
        && rejectsTopLevel;
    std::printf("[test_mod_raise] source_levels=%zu/%zu target_level=%zu %s\n",
                inputInfo.chainIndex, lowestInfo.chainIndex, resultInfo.chainIndex, ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
