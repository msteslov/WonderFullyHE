#include "m2424/homomorphic_linear_transform.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace m2424 {

HomomorphicLinearTransform::HomomorphicLinearTransform(
    std::size_t physicalSlotCount,
    std::vector<HomomorphicDiagonalTerm> terms)
    : physicalSlotCount_(physicalSlotCount), terms_(std::move(terms)) {
    if (physicalSlotCount_ == 0 || terms_.empty()) {
        throw std::invalid_argument("homomorphic linear transform requires slots and terms");
    }
    for (const auto& term : terms_) {
        if (term.diagonal.size() != physicalSlotCount_) {
            throw std::invalid_argument("linear transform diagonal size must equal physical slot count");
        }
    }
}

std::vector<int> HomomorphicLinearTransform::rotationSteps() const {
    std::vector<int> steps;
    for (const auto& term : terms_) {
        if (term.rotation != 0) {
            steps.push_back(term.rotation);
        }
    }
    std::sort(steps.begin(), steps.end());
    steps.erase(std::unique(steps.begin(), steps.end()), steps.end());
    return steps;
}

Cipher HomomorphicLinearTransform::apply(SealAdapter& adapter, const Cipher& input) const {
    if (adapter.slotCount() != physicalSlotCount_) {
        throw std::invalid_argument("linear transform slot count does not match ciphertext context");
    }
    const auto inputInfo = adapter.info(input);
    const auto coeffModulusBits = adapter.coeffModulusBits();
    if (inputInfo.chainIndex == 0 || inputInfo.coeffModulusSize == 0
        || inputInfo.coeffModulusSize > coeffModulusBits.size()) {
        throw std::runtime_error("linear transform requires one remaining rescale level");
    }
    const double diagonalScale = std::exp2(static_cast<double>(coeffModulusBits[inputInfo.coeffModulusSize - 1]));

    Cipher result;
    bool hasResult = false;
    for (const auto& term : terms_) {
        const Cipher rotated = term.rotation == 0 ? input : adapter.rotate(input, term.rotation);
        const Plain diagonal = adapter.encodeComplexAtScaleFor(term.diagonal, diagonalScale, rotated);
        const Cipher weighted = adapter.rescaleToNext(adapter.multiplyPlain(rotated, diagonal));
        if (!hasResult) {
            result = weighted;
            hasResult = true;
        } else {
            result = adapter.add(result, weighted);
        }
    }
    return result;
}

} // namespace m2424
