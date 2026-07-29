#include "m2424/coeff_to_slot_bsgs_plan.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <utility>

namespace m2424 {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

bool isPowerOfTwo(std::size_t value) {
    return value >= 2 && (value & (value - 1)) == 0;
}

int floorDivide(int value, int divisor) {
    int quotient = value / divisor;
    const int remainder = value % divisor;
    if (remainder < 0) {
        --quotient;
    }
    return quotient;
}

std::size_t wrappedIndex(int value, std::size_t modulus) {
    const int signedModulus = static_cast<int>(modulus);
    int result = value % signedModulus;
    if (result < 0) {
        result += signedModulus;
    }
    return static_cast<std::size_t>(result);
}

ComplexVector rotatePlain(const ComplexVector& values, int steps) {
    ComplexVector result(values.size());
    for (std::size_t index = 0; index < values.size(); ++index) {
        result[index] = values[wrappedIndex(static_cast<int>(index) + steps, values.size())];
    }
    return result;
}

} // namespace

CoeffToSlotBsgsPlan::CoeffToSlotBsgsPlan(std::size_t transformSlots,
                                         std::size_t physicalSlotCount,
                                         std::size_t babyStep)
    : transformSlots_(transformSlots), physicalSlotCount_(physicalSlotCount), babyStep_(babyStep) {
    if (!isPowerOfTwo(transformSlots_) || transformSlots_ > physicalSlotCount_ || babyStep_ == 0) {
        throw std::invalid_argument("invalid BSGS transform dimensions");
    }

    std::map<int, std::map<int, ComplexVector>> groupedDiagonals;
    for (std::size_t row = 0; row < transformSlots_; ++row) {
        for (std::size_t column = 0; column < transformSlots_; ++column) {
            const int rotation = static_cast<int>(column) - static_cast<int>(row);
            const int giantIndex = floorDivide(rotation, static_cast<int>(babyStep_));
            const int babyRotation = rotation - giantIndex * static_cast<int>(babyStep_);
            const int giantRotation = giantIndex * static_cast<int>(babyStep_);
            auto& diagonal = groupedDiagonals[giantRotation][babyRotation];
            if (diagonal.empty()) {
                diagonal.resize(physicalSlotCount_);
            }
            const double angle = -2.0 * kPi * static_cast<double>(row * column)
                / static_cast<double>(transformSlots_);
            const Complex dftValue{std::cos(angle), std::sin(angle)};
            // После giant rotation значение из позиции row + giantRotation попадает в row.
            // Поэтому коэффициент строки row храним в этой исходной позиции.
            const std::size_t shiftedRow = wrappedIndex(static_cast<int>(row) + giantRotation, physicalSlotCount_);
            diagonal[shiftedRow] += dftValue;
        }
    }

    groups_.reserve(groupedDiagonals.size());
    for (auto& [giantRotation, babyTerms] : groupedDiagonals) {
        GiantGroup group;
        group.giantRotation = giantRotation;
        group.terms.reserve(babyTerms.size());
        for (auto& [babyRotation, diagonal] : babyTerms) {
            group.terms.push_back({babyRotation, std::move(diagonal)});
        }
        groups_.push_back(std::move(group));
    }
}

std::size_t CoeffToSlotBsgsPlan::requiredLevels() const noexcept {
    return 1;
}

std::vector<int> CoeffToSlotBsgsPlan::rotationSteps() const {
    std::vector<int> steps;
    for (std::size_t groupIndex = 0; groupIndex < groups_.size(); ++groupIndex) {
        const auto& group = groups_[groupIndex];
        if (group.giantRotation != 0) {
            steps.push_back(group.giantRotation);
        }
        for (const auto& term : group.terms) {
            if (term.babyRotation != 0) {
                steps.push_back(term.babyRotation);
            }
        }
    }
    std::sort(steps.begin(), steps.end());
    steps.erase(std::unique(steps.begin(), steps.end()), steps.end());
    return steps;
}

ComplexVector CoeffToSlotBsgsPlan::applyPlain(const ComplexVector& input) const {
    if (input.size() != transformSlots_) {
        throw std::invalid_argument("BSGS plaintext input size must match transform slots");
    }
    ComplexVector padded(physicalSlotCount_);
    std::copy(input.begin(), input.end(), padded.begin());
    std::map<int, ComplexVector> babyValues;
    babyValues.emplace(0, padded);
    for (const auto& group : groups_) {
        for (const auto& term : group.terms) {
            if (babyValues.find(term.babyRotation) == babyValues.end()) {
                babyValues.emplace(term.babyRotation, rotatePlain(padded, term.babyRotation));
            }
        }
    }

    ComplexVector result(physicalSlotCount_);
    for (const auto& group : groups_) {
        ComplexVector inner(physicalSlotCount_);
        for (const auto& term : group.terms) {
            const auto& baby = babyValues.at(term.babyRotation);
            for (std::size_t index = 0; index < physicalSlotCount_; ++index) {
                inner[index] += term.diagonal[index] * baby[index];
            }
        }
        const auto shifted = rotatePlain(inner, group.giantRotation);
        for (std::size_t index = 0; index < physicalSlotCount_; ++index) {
            result[index] += shifted[index];
        }
    }
    return {result.begin(), result.begin() + static_cast<std::ptrdiff_t>(transformSlots_)};
}

void CoeffToSlotBsgsPlan::prepare(SealAdapter& adapter, const Cipher& input) const {
    if (adapter.slotCount() != physicalSlotCount_) {
        throw std::invalid_argument("BSGS slot count does not match ciphertext context");
    }
    const auto inputInfo = adapter.info(input);
    if (preparedAdapter_ == &adapter && preparedChainIndex_ == inputInfo.chainIndex
        && preparedInputScale_ == inputInfo.scale && preparedDiagonals_.size() == groups_.size()) {
        return;
    }
    const auto coeffModulusBits = adapter.coeffModulusBits();
    if (inputInfo.chainIndex == 0 || inputInfo.coeffModulusSize == 0
        || inputInfo.coeffModulusSize > coeffModulusBits.size()) {
        throw std::runtime_error("BSGS transform requires one remaining rescale level");
    }
    const double diagonalScale = std::exp2(static_cast<double>(coeffModulusBits[inputInfo.coeffModulusSize - 1]));

    std::vector<std::vector<Plain>> diagonals;
    diagonals.reserve(groups_.size());
    for (const auto& group : groups_) {
        std::vector<Plain> groupDiagonals;
        groupDiagonals.reserve(group.terms.size());
        for (const auto& term : group.terms) {
            groupDiagonals.push_back(adapter.encodeComplexAtScaleFor(term.diagonal, diagonalScale, input));
        }
        diagonals.push_back(std::move(groupDiagonals));
    }
    preparedAdapter_ = &adapter;
    preparedChainIndex_ = inputInfo.chainIndex;
    preparedInputScale_ = inputInfo.scale;
    preparedDiagonals_ = std::move(diagonals);
}

std::size_t CoeffToSlotBsgsPlan::preparedPlaintextBytes(SealAdapter& adapter) const {
    std::size_t result = 0;
    for (const auto& group : preparedDiagonals_) {
        for (const auto& diagonal : group) {
            result += adapter.serializedSize(diagonal);
        }
    }
    return result;
}

Cipher CoeffToSlotBsgsPlan::apply(SealAdapter& adapter, const Cipher& input) const {
    prepare(adapter, input);

    std::map<int, Cipher> babyValues;
    babyValues.emplace(0, input);
    for (std::size_t groupIndex = 0; groupIndex < groups_.size(); ++groupIndex) {
        const auto& group = groups_[groupIndex];
        for (const auto& term : group.terms) {
            if (babyValues.find(term.babyRotation) == babyValues.end()) {
                babyValues.emplace(term.babyRotation, adapter.rotate(input, term.babyRotation));
            }
        }
    }

    Cipher result;
    bool hasResult = false;
    for (std::size_t groupIndex = 0; groupIndex < groups_.size(); ++groupIndex) {
        const auto& group = groups_[groupIndex];
        Cipher inner;
        bool hasInner = false;
        for (std::size_t termIndex = 0; termIndex < group.terms.size(); ++termIndex) {
            const auto& term = group.terms[termIndex];
            const auto& baby = babyValues.at(term.babyRotation);
            const Cipher weighted = adapter.rescaleToNext(
                adapter.multiplyPlain(baby, preparedDiagonals_[groupIndex][termIndex]));
            if (!hasInner) {
                inner = weighted;
                hasInner = true;
            } else {
                inner = adapter.add(inner, weighted);
            }
        }
        const Cipher shifted = group.giantRotation == 0 ? inner : adapter.rotate(inner, group.giantRotation);
        if (!hasResult) {
            result = shifted;
            hasResult = true;
        } else {
            result = adapter.add(result, shifted);
        }
    }
    return result;
}

CoeffToSlotPlanRequirements CoeffToSlotBsgsPlan::requirements() const {
    return {requiredLevels(), rotationSteps()};
}

} // namespace m2424
