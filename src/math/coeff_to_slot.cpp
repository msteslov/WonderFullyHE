#include "m2424/coeff_to_slot.hpp"
#include "m2424/evalround_plus_coeff_to_slot.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <functional>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <utility>

namespace m2424 {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
using DiagonalMap = std::map<std::size_t, ComplexVector>;

std::size_t validatedDegree(std::size_t degree) {
    if (degree < 4 || (degree & (degree - 1)) != 0) {
        throw std::invalid_argument("CoeffToSlot degree must be a power of two >= 4");
    }
    return degree;
}

std::size_t log2Exact(std::size_t value) {
    std::size_t result = 0;
    while ((std::size_t{1} << result) < value) {
        ++result;
    }
    return result;
}

Complex root(std::size_t order, std::size_t exponent, bool inverse = true) {
    const double sign = inverse ? -1.0 : 1.0;
    const double angle = sign * 2.0 * kPi * static_cast<double>(exponent)
        / static_cast<double>(order);
    return {std::cos(angle), std::sin(angle)};
}

void addEntry(DiagonalMap& matrix, std::size_t slots, std::size_t row,
              std::size_t column, Complex value) {
    const std::size_t diagonal = (column + slots - row) % slots;
    auto [it, inserted] = matrix.emplace(diagonal, ComplexVector{});
    if (inserted) {
        it->second.assign(slots, Complex{});
    }
    it->second[row] += value;
}

DiagonalMap butterflyStage(std::size_t slots, std::size_t blockSize) {
    DiagonalMap result;
    const std::size_t half = blockSize / 2;
    for (std::size_t base = 0; base < slots; base += blockSize) {
        if (blockSize == 2) {
            addEntry(result, slots, base, base, 1.0);
            addEntry(result, slots, base, base + 1, 1.0);
            addEntry(result, slots, base + 1, base, root(8, 1));
            addEntry(result, slots, base + 1, base + 1, root(8, 3));
            continue;
        }
        std::size_t power = 1;
        for (std::size_t index = 0; index < half; ++index) {
            const Complex twiddle = root(4 * blockSize, power);
            addEntry(result, slots, base + index, base + index, 1.0);
            addEntry(result, slots, base + index, base + half + index, 1.0);
            addEntry(result, slots, base + half + index, base + index, twiddle);
            addEntry(result, slots, base + half + index, base + half + index, -twiddle);
            power = (power * 3) % (4 * blockSize);
        }
    }
    return result;
}

std::size_t swapBits(std::size_t value, std::size_t first, std::size_t second) {
    const bool a = ((value >> first) & 1U) != 0;
    const bool b = ((value >> second) & 1U) != 0;
    if (a != b) {
        value ^= (std::size_t{1} << first) | (std::size_t{1} << second);
    }
    return value;
}

DiagonalMap bitSwapStage(std::size_t slots, std::size_t first, std::size_t second) {
    DiagonalMap result;
    for (std::size_t row = 0; row < slots; ++row) {
        addEntry(result, slots, row, swapBits(row, first, second), 1.0);
    }
    return result;
}

DiagonalMap compose(const DiagonalMap& outer, const DiagonalMap& inner,
                    std::size_t slots) {
    DiagonalMap result;
    for (const auto& [outerRotation, outerDiagonal] : outer) {
        for (const auto& [innerRotation, innerDiagonal] : inner) {
            const std::size_t rotation = (outerRotation + innerRotation) % slots;
            auto [it, inserted] = result.emplace(rotation, ComplexVector{});
            if (inserted) {
                it->second.assign(slots, Complex{});
            }
            for (std::size_t row = 0; row < slots; ++row) {
                it->second[row] += outerDiagonal[row]
                    * innerDiagonal[(row + outerRotation) % slots];
            }
        }
    }
    for (auto it = result.begin(); it != result.end();) {
        const bool zero = std::all_of(it->second.begin(), it->second.end(),
            [](Complex value) { return std::abs(value) < 1e-15; });
        if (zero) {
            it = result.erase(it);
        } else {
            ++it;
        }
    }
    return result;
}

ComplexVector applyMatrix(const DiagonalMap& matrix, const ComplexVector& input) {
    ComplexVector output(input.size());
    for (const auto& [rotation, diagonal] : matrix) {
        for (std::size_t row = 0; row < input.size(); ++row) {
            output[row] += diagonal[row] * input[(row + rotation) % input.size()];
        }
    }
    return output;
}

std::vector<std::size_t> balancedRadices(std::size_t factors, std::size_t depth) {
    depth = std::max<std::size_t>(1, std::min(depth, factors));
    std::vector<std::size_t> result(depth, factors / depth);
    for (std::size_t index = 0; index < factors % depth; ++index) {
        ++result[index];
    }
    return result;
}

std::vector<std::set<std::size_t>> rawFactorSupports(std::size_t slots) {
    std::vector<std::set<std::size_t>> result;
    for (std::size_t block = slots; block >= 2; block /= 2) {
        const std::size_t offset = block / 2;
        result.push_back({0, offset, (slots - offset) % slots});
    }
    const std::size_t bits = log2Exact(slots);
    for (std::size_t low = 0; low < bits / 2; ++low) {
        const std::size_t high = bits - 1 - low;
        const std::size_t offset = (std::size_t{1} << high) - (std::size_t{1} << low);
        result.push_back({0, offset, slots - offset});
    }
    return result;
}

std::set<std::size_t> composeSupports(const std::set<std::size_t>& outer,
                                      const std::set<std::size_t>& inner,
                                      std::size_t slots) {
    std::set<std::size_t> result;
    for (const std::size_t lhs : outer) {
        for (const std::size_t rhs : inner) result.insert((lhs + rhs) % slots);
    }
    return result;
}

std::set<int> supportBsgsRotationSteps(const std::set<std::size_t>& support,
                                       std::size_t babyStep) {
    std::set<int> result;
    for (const std::size_t rotation : support) {
        const std::size_t baby = rotation % babyStep;
        const std::size_t giant = rotation - baby;
        if (baby != 0) result.insert(static_cast<int>(baby));
        if (giant != 0) result.insert(static_cast<int>(giant));
    }
    return result;
}

std::size_t supportBsgsRotationOperationCount(
    const std::set<std::size_t>& support,
    std::size_t babyStep) {
    std::set<std::size_t> babies;
    std::set<std::size_t> giants;
    for (const std::size_t rotation : support) {
        const std::size_t baby = rotation % babyStep;
        const std::size_t giant = rotation - baby;
        if (baby != 0) babies.insert(baby);
        if (giant != 0) giants.insert(giant);
    }
    return babies.size() + giants.size();
}

std::size_t selectSupportBabyStep(const std::set<std::size_t>& support,
                                  std::size_t slots) {
    std::size_t bestStep = 1;
    std::size_t bestCost = support.size();
    for (std::size_t candidate = 1; candidate <= slots; candidate *= 2) {
        const std::size_t cost = supportBsgsRotationOperationCount(support, candidate);
        if (cost < bestCost) {
            bestCost = cost;
            bestStep = candidate;
        }
        if (candidate > slots / 2) break;
    }
    return bestStep;
}

std::vector<DiagonalMap> mergeFactors(const std::vector<DiagonalMap>& raw,
                                      const std::vector<std::size_t>& radices,
                                      std::size_t slots) {
    std::vector<DiagonalMap> result;
    std::size_t cursor = 0;
    for (const std::size_t radix : radices) {
        DiagonalMap grouped = raw[cursor++];
        for (std::size_t index = 1; index < radix; ++index) {
            grouped = compose(raw[cursor++], grouped, slots);
        }
        result.push_back(std::move(grouped));
    }
    return result;
}

DiagonalMap withInputMultiplier(DiagonalMap matrix, const ComplexVector& multiplier) {
    for (auto& [rotation, diagonal] : matrix) {
        for (std::size_t row = 0; row < diagonal.size(); ++row) {
            diagonal[row] *= multiplier[(row + rotation) % diagonal.size()];
        }
    }
    return matrix;
}

DiagonalMap prefactored(DiagonalMap factor, const CoeffToSlotPrefactor& scalar) {
    if (scalar.isIdentity()) return factor;
    for (auto& entry : factor) for (auto& value : entry.second) {
        const auto real = scalar.multiplyRounded(value.real());
        const auto imag = scalar.multiplyRounded(value.imag());
        if ((value.real() != 0 && real == 0) || (value.imag() != 0 && imag == 0))
            throw std::underflow_error("folded diagonal is not representable in binary64");
        value = {real, imag};
    }
    return factor;
}

struct PreparedBsgsTerm {
    std::size_t babyRotation{};
    Plain extendedDiagonal;
};

struct PreparedBsgsGroup {
    std::size_t giantRotation{};
    std::vector<PreparedBsgsTerm> terms;
};

using PreparedFactor = std::vector<PreparedBsgsGroup>;

std::vector<int> bsgsRotationSteps(const DiagonalMap& factor,
                                   std::size_t babyStep) {
    std::set<int> result;
    for (const auto& [rotation, diagonal] : factor) {
        (void)diagonal;
        const std::size_t baby = rotation % babyStep;
        const std::size_t giant = rotation - baby;
        if (baby != 0) result.insert(static_cast<int>(baby));
        if (giant != 0) result.insert(static_cast<int>(giant));
    }
    return {result.begin(), result.end()};
}

std::size_t bsgsRotationOperationCount(const DiagonalMap& factor,
                                       std::size_t babyStep) {
    std::set<std::size_t> babies;
    std::set<std::size_t> giants;
    for (const auto& [rotation, diagonal] : factor) {
        (void)diagonal;
        const std::size_t baby = rotation % babyStep;
        const std::size_t giant = rotation - baby;
        if (baby != 0) babies.insert(baby);
        if (giant != 0) giants.insert(giant);
    }
    return babies.size() + giants.size();
}

std::size_t selectBabyStep(const DiagonalMap& factor, std::size_t slots) {
    std::size_t bestStep = 1;
    std::size_t bestCost = factor.size();
    for (std::size_t candidate = 1; candidate <= slots; ++candidate) {
        const std::size_t cost = bsgsRotationOperationCount(factor, candidate);
        if (cost < bestCost) {
            bestCost = cost;
            bestStep = candidate;
        }
    }
    return bestStep;
}

Cipher applyCipherMatrix(SealAdapter& adapter, const PreparedFactor& factor,
                         const Cipher& input) {
    std::vector<HoistedBsgsGroup> groups;
    groups.reserve(factor.size());
    for (const auto& group : factor) {
        HoistedBsgsGroup backendGroup;
        backendGroup.giantRotation = static_cast<int>(group.giantRotation);
        backendGroup.terms.reserve(group.terms.size());
        for (const auto& term : group.terms) {
            backendGroup.terms.push_back({
                static_cast<int>(term.babyRotation), &term.extendedDiagonal});
        }
        groups.push_back(std::move(backendGroup));
    }
    return adapter.rescaleToNext(
        adapter.applyBsgsDoubleHoisted(input, groups));
}

} // namespace

struct CoeffToSlotPlan::Impl {
    explicit Impl(std::size_t requestedDegree,
                  std::size_t targetDepth,
                  CoeffToSlotFactorization requestedFactorization = {})
        : degree(validatedDegree(requestedDegree)), slots(degree / 2) {
        const std::size_t logSlots = log2Exact(slots);
        std::vector<DiagonalMap> raw;
        for (std::size_t block = slots; block >= 2; block /= 2) {
            raw.push_back(butterflyStage(slots, block));
        }
        for (std::size_t low = 0; low < logSlots / 2; ++low) {
            raw.push_back(bitSwapStage(slots, low, logSlots - 1 - low));
        }
        rawStageCount = raw.size();
        if (requestedFactorization.radices.empty()) {
            const auto tuned = CoeffToSlotPlan::knownTunedFactorization(degree, targetDepth);
            factorization = tuned.value_or(
                CoeffToSlotFactorization{balancedRadices(raw.size(), targetDepth)});
        } else {
            const bool valid = std::all_of(requestedFactorization.radices.begin(),
                                           requestedFactorization.radices.end(),
                                           [](std::size_t radix) { return radix != 0; })
                && std::accumulate(requestedFactorization.radices.begin(),
                                   requestedFactorization.radices.end(), std::size_t{0})
                    == raw.size();
            if (!valid) {
                throw std::invalid_argument("invalid CoeffToSlot factorization");
            }
            factorization = std::move(requestedFactorization);
        }
        factors = mergeFactors(raw, factorization.radices, slots);
        babySteps.reserve(factors.size());
        for (const auto& factor : factors) {
            babySteps.push_back(selectBabyStep(factor, slots));
        }

        ComplexVector secondMultiplier(slots);
        std::size_t exponent = 1;
        for (std::size_t index = 0; index < slots; ++index) {
            secondMultiplier[index] = root(4, exponent % 4);
            exponent = (exponent * 3) % (4 * slots);
        }
        secondFirstFactor = withInputMultiplier(factors.front(), secondMultiplier);

        const double normalization = 1.0 / static_cast<double>(degree);
        for (auto& [rotation, diagonal] : factors.front()) {
            (void)rotation;
            for (Complex& value : diagonal) {
                value *= normalization;
            }
        }
        for (auto& [rotation, diagonal] : secondFirstFactor) {
            (void)rotation;
            for (Complex& value : diagonal) {
                value *= normalization;
            }
        }
    }

    std::size_t degree{};
    std::size_t slots{};
    std::size_t rawStageCount{};
    CoeffToSlotFactorization factorization;
    std::vector<DiagonalMap> factors;
    std::vector<std::size_t> babySteps;
    DiagonalMap secondFirstFactor;
};

struct PreparedCoeffToSlotPlan::Impl {
    CoeffToSlotPrefactor prefactor;
    std::size_t planDegree{};
    CoeffToSlotFactorization factorization;
    std::array<std::uint64_t, 4> contextFingerprint{};
    std::array<std::uint64_t, 4> parmsFingerprint{};
    std::size_t startChainIndex{};
    double inputScale{};
    double contractInputScaleLog2{};
    double contractOutputScaleLog2{};
    std::vector<PreparedFactor> factors;
    PreparedFactor secondFirstFactor;
};

PreparedCoeffToSlotPlan::PreparedCoeffToSlotPlan(std::unique_ptr<Impl> implementation)
    : pimpl_(std::move(implementation)) {}
PreparedCoeffToSlotPlan::~PreparedCoeffToSlotPlan() = default;
PreparedCoeffToSlotPlan::PreparedCoeffToSlotPlan(PreparedCoeffToSlotPlan&&) noexcept = default;
PreparedCoeffToSlotPlan& PreparedCoeffToSlotPlan::operator=(PreparedCoeffToSlotPlan&&) noexcept = default;

std::size_t PreparedCoeffToSlotPlan::plaintextCount() const {
    std::size_t result = 0;
    const auto countFactor = [](const PreparedFactor& factor) {
        std::size_t count = 0;
        for (const auto& group : factor) count += group.terms.size();
        return count;
    };
    result += countFactor(pimpl_->secondFirstFactor);
    for (const auto& factor : pimpl_->factors) result += countFactor(factor);
    return result;
}

std::size_t PreparedCoeffToSlotPlan::serializedPlaintextBytes(const SealAdapter& adapter) const {
    std::size_t result = 0;
    for (const auto& factor : pimpl_->factors) {
        for (const auto& group : factor) {
            for (const auto& term : group.terms) {
                result += adapter.serializedSize(term.extendedDiagonal);
            }
        }
    }
    for (const auto& group : pimpl_->secondFirstFactor) {
        for (const auto& term : group.terms) {
            result += adapter.serializedSize(term.extendedDiagonal);
        }
    }
    return result;
}

CoeffToSlotPlan::CoeffToSlotPlan(std::size_t degree, std::size_t targetDepth)
    : pimpl_(std::make_unique<Impl>(degree, targetDepth)) {}
CoeffToSlotPlan::CoeffToSlotPlan(std::size_t degree, CoeffToSlotFactorization factorization)
    : pimpl_(std::make_unique<Impl>(degree, factorization.radices.size(),
                                    std::move(factorization))) {}
CoeffToSlotPlan::~CoeffToSlotPlan() = default;
CoeffToSlotPlan::CoeffToSlotPlan(CoeffToSlotPlan&&) noexcept = default;
CoeffToSlotPlan& CoeffToSlotPlan::operator=(CoeffToSlotPlan&&) noexcept = default;

std::size_t CoeffToSlotPlan::polyModulusDegree() const { return pimpl_->degree; }
std::size_t CoeffToSlotPlan::butterflyStageCount() const {
    return log2Exact(pimpl_->slots);
}
std::size_t CoeffToSlotPlan::rawStageCount() const { return pimpl_->rawStageCount; }
std::size_t CoeffToSlotPlan::depth() const { return pimpl_->factors.size(); }
const CoeffToSlotFactorization& CoeffToSlotPlan::factorization() const {
    return pimpl_->factorization;
}

std::vector<std::vector<int>> CoeffToSlotPlan::stageRotationSteps() const {
    std::vector<std::vector<int>> result;
    result.reserve(pimpl_->factors.size());
    for (std::size_t index = 0; index < pimpl_->factors.size(); ++index) {
        result.push_back(bsgsRotationSteps(
            pimpl_->factors[index], pimpl_->babySteps[index]));
    }
    return result;
}

CoeffToSlotPlanRequirements CoeffToSlotPlan::requirements() const {
    std::set<int> rotations;
    for (std::size_t index = 0; index < pimpl_->factors.size(); ++index) {
        const auto steps = bsgsRotationSteps(
            pimpl_->factors[index], pimpl_->babySteps[index]);
        rotations.insert(steps.begin(), steps.end());
    }
    return {depth(), std::vector<int>(rotations.begin(), rotations.end()), true};
}

CoeffToSlotPlanMetrics CoeffToSlotPlan::metrics() const {
    CoeffToSlotPlanMetrics result;
    result.butterflyStages = butterflyStageCount();
    result.permutationStages = rawStageCount() - result.butterflyStages;
    result.depth = depth();
    for (std::size_t index = 0; index < pimpl_->factors.size(); ++index) {
        const auto& factor = pimpl_->factors[index];
        const std::size_t diagonals = factor.size();
        result.diagonalsPerStage.push_back(diagonals);
        result.plaintextMultiplicationsPerApply += 2 * diagonals;
        result.rotationsPerApply += 2 * bsgsRotationOperationCount(
            factor, pimpl_->babySteps[index]);
        std::set<std::size_t> babies;
        std::set<std::size_t> giants;
        for (const auto& [rotation, diagonal] : factor) {
            (void)diagonal;
            const auto baby = rotation % pimpl_->babySteps[index];
            const auto giant = rotation - baby;
            if (baby != 0) babies.insert(baby);
            if (giant != 0) giants.insert(giant);
        }
        result.hoistedDecompositionsPerApply += 2 * (1 + giants.size());
        result.hoistedAutomorphismsPerApply += 2 * (babies.size() + giants.size());
        result.innerModDownsPerApply += 2 * (giants.size() + 1);
        result.finalModDownsPerApply += 2;
        result.additionsPerApply += 2 * (diagonals - 1);
        result.storedComplexValues += diagonals * pimpl_->slots;
    }
    result.additionsPerApply += 2;
    result.rescalesPerApply = 2 * depth();
    result.uniqueEvaluationKeys = requirements().rotationSteps.size() + 1;
    result.storedComplexValues += pimpl_->secondFirstFactor.size() * pimpl_->slots;
    return result;
}

std::optional<CoeffToSlotFactorization> CoeffToSlotPlan::knownTunedFactorization(
    std::size_t polyModulusDegree, std::size_t depth) {
    if (polyModulusDegree == 16384 && depth == 4) {
        return CoeffToSlotFactorization{{6, 5, 5, 3}};
    }
    return std::nullopt;
}

std::vector<CoeffToSlotFactorizationEstimate> CoeffToSlotPlan::estimateFactorizations(
    std::size_t polyModulusDegree,
    std::size_t requestedDepth,
    std::size_t limit) {
    const std::size_t degree = validatedDegree(polyModulusDegree);
    const std::size_t slots = degree / 2;
    const auto supports = rawFactorSupports(slots);
    if (requestedDepth == 0 || requestedDepth > supports.size() || limit == 0) {
        throw std::invalid_argument("invalid CoeffToSlot factorization search parameters");
    }

    std::vector<CoeffToSlotFactorizationEstimate> estimates;
    std::vector<std::size_t> current;
    const auto enumerate = [&](const auto& self, std::size_t remaining,
                               std::size_t parts) -> void {
        if (parts == 1) {
            current.push_back(remaining);
            std::size_t cursor = 0;
            std::size_t diagonalCount = 0;
            std::size_t rotationCount = 0;
            std::set<int> rotations;
            for (const std::size_t radix : current) {
                auto group = supports[cursor++];
                for (std::size_t index = 1; index < radix; ++index) {
                    group = composeSupports(supports[cursor++], group, slots);
                }
                diagonalCount += group.size();
                const auto babyStep = selectSupportBabyStep(group, slots);
                const auto groupRotations = supportBsgsRotationSteps(group, babyStep);
                rotationCount += supportBsgsRotationOperationCount(group, babyStep);
                rotations.insert(groupRotations.begin(), groupRotations.end());
            }
            estimates.push_back({CoeffToSlotFactorization{current},
                                 diagonalCount, rotationCount, rotations.size()});
            current.pop_back();
            return;
        }
        for (std::size_t value = 1; value <= remaining - (parts - 1); ++value) {
            current.push_back(value);
            self(self, remaining - value, parts - 1);
            current.pop_back();
        }
    };
    enumerate(enumerate, supports.size(), requestedDepth);
    std::sort(estimates.begin(), estimates.end(),
        [](const auto& lhs, const auto& rhs) {
            constexpr std::size_t rotationWeight = 3;
            const std::size_t lhsCost = lhs.estimatedDiagonals
                + rotationWeight * lhs.estimatedRotationsPerOutput;
            const std::size_t rhsCost = rhs.estimatedDiagonals
                + rotationWeight * rhs.estimatedRotationsPerOutput;
            if (lhsCost != rhsCost) {
                return lhsCost < rhsCost;
            }
            if (lhs.estimatedRotationsPerOutput != rhs.estimatedRotationsPerOutput) {
                return lhs.estimatedRotationsPerOutput < rhs.estimatedRotationsPerOutput;
            }
            if (lhs.estimatedUniqueRotations != rhs.estimatedUniqueRotations) {
                return lhs.estimatedUniqueRotations < rhs.estimatedUniqueRotations;
            }
            return lhs.factorization.radices < rhs.factorization.radices;
        });
    if (estimates.size() > limit) estimates.resize(limit);
    return estimates;
}

std::pair<ComplexVector, ComplexVector>
CoeffToSlotPlan::applyPlain(const ComplexVector& slots) const {
    return applyPlain(slots, CoeffToSlotPrefactor{});
}

std::pair<ComplexVector, ComplexVector>
CoeffToSlotPlan::applyPlain(const ComplexVector& slots, const CoeffToSlotPrefactor& scalar) const {
    if (slots.size() != pimpl_->slots) {
        throw std::invalid_argument("CoeffToSlotPlan plaintext slot count mismatch");
    }
    ComplexVector first = applyMatrix(prefactored(pimpl_->factors.front(), scalar), slots);
    ComplexVector second = applyMatrix(prefactored(pimpl_->secondFirstFactor, scalar), slots);
    for (std::size_t index = 1; index < pimpl_->factors.size(); ++index) {
        first = applyMatrix(pimpl_->factors[index], first);
        second = applyMatrix(pimpl_->factors[index], second);
    }
    for (std::size_t index = 0; index < pimpl_->slots; ++index) {
        first[index] += std::conj(first[index]);
        second[index] += std::conj(second[index]);
    }
    return {std::move(first), std::move(second)};
}

PreparedCoeffToSlotPlan CoeffToSlotPlan::prepare(
    SealAdapter& adapter,
    const RaisedCipher& input,
    const CoeffToSlotContract& contract) const {
    return prepare(adapter, input, contract, CoeffToSlotPrefactor{});
}

PreparedCoeffToSlotPlan CoeffToSlotPlan::prepare(
    SealAdapter& adapter, const RaisedCipher& input, const CoeffToSlotContract& contract,
    const CoeffToSlotPrefactor& scalar) const {
    if (pimpl_->degree != contract.polyModulusDegree) {
        throw std::invalid_argument("CoeffToSlot plan degree does not match contract");
    }
    const auto requirements = this->requirements();
    const auto preflight = preflightCoeffToSlot(adapter, input, contract, requirements);
    if (!preflight.ready) {
        throw std::invalid_argument("CoeffToSlot preflight failed: " + preflight.blocker);
    }

    const auto info = adapter.info(input);
    std::vector<PreparedFactor> preparedFactors;
    PreparedFactor preparedSecond;
    preparedFactors.reserve(pimpl_->factors.size());
    const auto encodeFactor = [&](const DiagonalMap& factor,
                                  std::size_t babyStep,
                                  std::size_t chainIndex) {
        std::map<std::size_t, PreparedBsgsGroup> grouped;
        const double plaintextScale =
            adapter.rescalePlaintextScaleAtChainIndex(chainIndex);
        for (const auto& [rotation, diagonal] : factor) {
            const std::size_t baby = rotation % babyStep;
            const std::size_t giant = rotation - baby;
            ComplexVector adjusted(diagonal.size());
            for (std::size_t row = 0; row < diagonal.size(); ++row) {
                adjusted[(row + giant) % diagonal.size()] = diagonal[row];
            }
            auto& group = grouped[giant];
            group.giantRotation = giant;
            group.terms.push_back({
                baby,
                adapter.encodeComplexAtKeyScale(adjusted, plaintextScale)});
        }
        PreparedFactor result;
        result.reserve(grouped.size());
        for (auto& entry : grouped) result.push_back(std::move(entry.second));
        return result;
    };
    for (std::size_t index = 0; index < pimpl_->factors.size(); ++index) {
        preparedFactors.push_back(encodeFactor(
            index == 0 ? prefactored(pimpl_->factors[index], scalar) : pimpl_->factors[index],
            pimpl_->babySteps[index], info.chainIndex - index));
    }
    preparedSecond = encodeFactor(
        prefactored(pimpl_->secondFirstFactor, scalar), pimpl_->babySteps.front(), info.chainIndex);

    auto result = std::make_unique<PreparedCoeffToSlotPlan::Impl>();
    result->prefactor = scalar;
    result->planDegree = pimpl_->degree;
    result->factorization = pimpl_->factorization;
    result->factors = std::move(preparedFactors);
    result->secondFirstFactor = std::move(preparedSecond);
    result->contextFingerprint = adapter.contextFingerprint();
    result->parmsFingerprint = adapter.parmsFingerprint(input);
    result->startChainIndex = info.chainIndex;
    result->inputScale = info.scale;
    result->contractInputScaleLog2 = contract.inputScaleLog2;
    result->contractOutputScaleLog2 = contract.outputScaleLog2;
    return PreparedCoeffToSlotPlan(std::move(result));
}

bool CoeffToSlotPlan::isPreparedFor(const PreparedCoeffToSlotPlan& prepared,
                                    const SealAdapter& adapter,
                                    const RaisedCipher& input,
                                    const CoeffToSlotContract& contract) const {
    return isPreparedFor(prepared, adapter, input, contract, CoeffToSlotPrefactor{});
}

bool CoeffToSlotPlan::isPreparedFor(const PreparedCoeffToSlotPlan& prepared,
    const SealAdapter& adapter, const RaisedCipher& input, const CoeffToSlotContract& contract,
    const CoeffToSlotPrefactor& scalar) const {
    if (!prepared.pimpl_ || pimpl_->degree != contract.polyModulusDegree) {
        return false;
    }
    const auto& state = *prepared.pimpl_;
    if (!(state.prefactor == scalar) || state.planDegree != pimpl_->degree
        || state.factorization.radices != pimpl_->factorization.radices
        || state.contextFingerprint != adapter.contextFingerprint()) {
        return false;
    }
    const auto info = adapter.info(input);
    return state.parmsFingerprint == adapter.parmsFingerprint(input)
        && state.startChainIndex == info.chainIndex
        && std::abs(std::log2(state.inputScale) - std::log2(info.scale)) < 1e-9
        && state.contractInputScaleLog2 == contract.inputScaleLog2
        && state.contractOutputScaleLog2 == contract.outputScaleLog2;
}

CoeffToSlot::CoeffToSlot(std::size_t degree, std::size_t targetDepth)
    : plan_(degree, targetDepth) {}
CoeffToSlot::CoeffToSlot(std::size_t degree, CoeffToSlotFactorization factorization)
    : plan_(degree, std::move(factorization)) {}
CoeffToSlot::~CoeffToSlot() = default;
CoeffToSlot::CoeffToSlot(CoeffToSlot&&) noexcept = default;
CoeffToSlot& CoeffToSlot::operator=(CoeffToSlot&&) noexcept = default;

const CoeffToSlotPlan& CoeffToSlot::plan() const { return plan_; }
CoeffToSlotPlanRequirements CoeffToSlot::requirements() const {
    return plan_.requirements();
}

PreparedCoeffToSlotPlan CoeffToSlot::prepare(
    SealAdapter& adapter,
    const RaisedCipher& input,
    const CoeffToSlotContract& contract) const {
    return plan_.prepare(adapter, input, contract);
}

CoeffToSlotResult CoeffToSlot::apply(SealAdapter& adapter,
                                     RaisedCipher&& input,
                                     const CoeffToSlotContract& contract,
                                     const PreparedCoeffToSlotPlan& prepared) const {
    const auto& plan = *plan_.pimpl_;
    if (plan.degree != contract.polyModulusDegree) {
        throw std::invalid_argument("CoeffToSlot plan degree does not match contract");
    }
    const auto preflight = preflightCoeffToSlot(adapter, input, contract, requirements());
    if (!preflight.ready) {
        throw std::invalid_argument("CoeffToSlot preflight failed: " + preflight.blocker);
    }
    if (!plan_.isPreparedFor(prepared, adapter, input, contract)) {
        throw std::logic_error("CoeffToSlot is not prepared for this context, level and scale");
    }
    Cipher source = std::move(input.cipher_);
    Cipher first = source;
    Cipher second = source;
    for (std::size_t index = 0; index < plan.factors.size(); ++index) {
        first = applyCipherMatrix(adapter, prepared.pimpl_->factors[index], first);
        second = applyCipherMatrix(adapter,
            index == 0 ? prepared.pimpl_->secondFirstFactor : prepared.pimpl_->factors[index],
            second);
    }
    first = adapter.add(first, adapter.conjugate(first));
    second = adapter.add(second, adapter.conjugate(second));
    const auto firstInfo = adapter.info(first);
    const auto secondInfo = adapter.info(second);
    if (std::abs(std::log2(firstInfo.scale) - contract.outputScaleLog2)
            > contract.inputScaleToleranceLog2
        || std::abs(std::log2(secondInfo.scale) - contract.outputScaleLog2)
            > contract.inputScaleToleranceLog2) {
        throw std::runtime_error("CoeffToSlot output scale violates contract");
    }
    return {std::move(first), std::move(second)};
}

namespace {
BootstrapBound unavailable(const std::string& reason) {
    BootstrapBound result;
    result.provenance = reason;
    return result;
}
BootstrapBound deterministicBound(double value, const std::string& reason) {
    return {value, BootstrapBoundKind::Deterministic, reason, {}};
}
bool finiteBound(const BootstrapBound& bound) {
    return bound.kind == BootstrapBoundKind::Deterministic && std::isfinite(bound.upperBound)
        && bound.upperBound >= 0 && !bound.provenance.empty() && bound.failureEventIds.empty();
}
double roundUp(double value) {
    return std::nextafter(value, std::numeric_limits<double>::infinity());
}
BootstrapBound magnitudeStep(const BootstrapBound& kappa, const BootstrapBound& magnitude) {
    if (!finiteBound(kappa) || !finiteBound(magnitude))
        return unavailable("v9 (6): missing deterministic input magnitude or ideal operator norm");
    const double value = roundUp(kappa.upperBound * magnitude.upperBound);
    if (!std::isfinite(value)) return unavailable("v9 (6): magnitude upper bound overflow");
    return deterministicBound(value, "v9 (6): M_l <= kappa_l M_(l-1), rounded upward");
}
BootstrapBound errorStep(const BootstrapBound& kappa, const BootstrapBound& delta,
    const BootstrapBound& magnitude, const BootstrapBound& error, const BootstrapBound& local) {
    if (!finiteBound(kappa) || !finiteBound(delta) || !finiteBound(magnitude)
        || !finiteBound(error) || !finiteBound(local))
        return unavailable("v9 (6): kappa*E + delta*(M+E) + B_loc; "
            "rigorous diagonal-encoding and double-hoisted BSGS local bounds unavailable");
    const double encoded = roundUp(delta.upperBound * roundUp(magnitude.upperBound + error.upperBound));
    const double propagated = roundUp(kappa.upperBound * error.upperBound);
    const double result = roundUp(roundUp(propagated + encoded) + local.upperBound);
    if (!std::isfinite(result)) return unavailable("v9 (6): propagated error overflow");
    return deterministicBound(result, "v9 (6): kappa*E + delta*(M+E) + B_loc, rounded upward");
}
std::uint64_t binaryScale(double scale) {
    std::uint64_t bits;
    std::memcpy(&bits, &scale, sizeof(bits));
    return bits;
}
} // namespace

PreparedEvalRoundPlusCoeffToSlot::PreparedEvalRoundPlusCoeffToSlot(
    PreparedCoeffToSlotPlan hp, PreparedCoeffToSlotPlan lp, BootstrapInputContext source,
    CoeffToSlotContract hpContract, CoeffToSlotContract lpContract)
    : hp_(std::move(hp)), lp_(std::move(lp)), source_(std::move(source)),
      hpContract_(std::move(hpContract)), lpContract_(std::move(lpContract)) {}

EvalRoundPlusCoeffToSlot::EvalRoundPlusCoeffToSlot(std::size_t degree, std::size_t depth)
    : plan_(degree, depth) {}

BootstrapContractResult EvalRoundPlusCoeffToSlot::preflight(const SealAdapter& adapter,
    const RaisedCipher& input, const BootstrapInputContext& source,
    const CoeffToSlotContract& hp, const CoeffToSlotContract& lp) const {
    using Status = BootstrapCertificationStatus;
    try {
        if (source.contextFingerprint != adapter.contextFingerprint())
            return {Status::MissingExactModulusContext, "coeff_to_slot.context", "Source context fingerprint mismatch"};
        const auto topPrimes = adapter.dataModulusValues();
        if (source.sourcePrimes.empty() || source.sourcePrimes.size() != input.sourceCoeffModulusSize_
            || source.sourcePrimes.size() > topPrimes.size() || source.raisedPrimes != topPrimes
            || !std::equal(source.sourcePrimes.begin(), source.sourcePrimes.end(), topPrimes.begin())
            || source.chainIndex != source.sourcePrimes.size() - 1
            || source.specialPrime != adapter.specialKeyModulusValue())
            return {Status::MissingExactModulusContext, "coeff_to_slot.source", "Source modulus does not match the ModRaise source base"};
        const auto state = adapter.info(input);
        if (source.scaleBinary64Bits != binaryScale(state.scale))
            return {Status::InputScaleMismatch, "coeff_to_slot.scale", "Exact Delta0 differs from the raised input scale"};
        if (state.coeffModulusSize != topPrimes.size())
            return {Status::MissingExactModulusContext, "coeff_to_slot.raised", "Input is not at the raised top level"};
        for (const auto* contract : {&hp, &lp}) {
            if (contract->polyModulusDegree != plan_.polyModulusDegree())
                return {Status::InvalidInput, "coeff_to_slot.degree", "Plan and branch contract degrees differ"};
            const auto result = preflightCoeffToSlot(adapter, input, *contract, requirements());
            if (!result.ready) {
                const auto status = !result.inputScaleMatches ? Status::InputScaleMismatch
                    : !result.levelsAvailable ? Status::InsufficientLevels
                    : !result.physicalSlotsAvailable ? Status::InvalidInput : Status::MissingEvaluationKeys;
                return {status, "coeff_to_slot.preflight", result.blocker};
            }
        }
        return {Status::Certified, {}, {}}; // Input compatibility only, not a branch certificate.
    } catch (const std::exception& error) {
        return {Status::InvalidInput, "coeff_to_slot.preflight", error.what()};
    }
}

PreparedEvalRoundPlusCoeffToSlot EvalRoundPlusCoeffToSlot::prepare(SealAdapter& adapter,
    const RaisedCipher& input, const BootstrapInputContext& source,
    const CoeffToSlotContract& hp, const CoeffToSlotContract& lp) const {
    const auto result = preflight(adapter, input, source, hp, lp);
    if (result.status != BootstrapCertificationStatus::Certified)
        throw std::invalid_argument(result.gate + ": " + result.provenance);
    auto hpPlan = plan_.prepare(adapter, input, hp);
    auto lpPlan = plan_.prepare(adapter, input, lp, CoeffToSlotPrefactor::sourceNormalization(source));
    return {std::move(hpPlan), std::move(lpPlan), source, hp, lp};
}

EvalRoundPlusCoeffToSlotResult EvalRoundPlusCoeffToSlot::apply(SealAdapter& adapter,
    const RaisedCipher& input, const PreparedEvalRoundPlusCoeffToSlot& prepared,
    const BootstrapBound& inputMagnitude) const {
    const auto ready = preflight(adapter, input, prepared.source_, prepared.hpContract_, prepared.lpContract_);
    if (ready.status != BootstrapCertificationStatus::Certified)
        throw std::invalid_argument(ready.gate + ": " + ready.provenance);
    const auto alpha = CoeffToSlotPrefactor::sourceNormalization(prepared.source_);
    if (!plan_.isPreparedFor(prepared.hp_, adapter, input, prepared.hpContract_)
        || !plan_.isPreparedFor(prepared.lp_, adapter, input, prepared.lpContract_, alpha))
        throw std::logic_error("EvalRound+ plans do not match context, input, factorization or prefactor");
    const auto start = adapter.info(input);
    const auto branch = [&](const PreparedCoeffToSlotPlan& constants, const CoeffToSlotContract& contract,
                            BootstrapGate gate, const CoeffToSlotPrefactor& scalar) {
        CoeffToSlotBranchTrace trace;
        trace.gate = gate;
        trace.prefactor = scalar;
        trace.inputMagnitude = inputMagnitude;
        std::array<Cipher, 2> outputs{input.cipher_, input.cipher_};
        for (std::size_t half = 0; half < 2; ++half) {
            auto magnitude = inputMagnitude;
            // Exact input semantics y=sigma(u/Delta0) includes nu_b (v9 3.5).
            auto error = deterministicBound(0, "v9 3.5: u includes input noise; no prior transform error");
            std::size_t rawOffset = 0;
            for (std::size_t index = 0; index < plan_.depth(); ++index) {
                const auto rawCount = plan_.factorization().radices[index];
                const auto butterflies = std::min(rawCount,
                    rawOffset < plan_.butterflyStageCount() ? plan_.butterflyStageCount() - rawOffset : 0);
                rawOffset += rawCount;
                double norm = std::ldexp(1.0, static_cast<int>(butterflies));
                if (index == 0) norm /= static_cast<double>(plan_.polyModulusDegree());
                if (index == 0 && !scalar.isIdentity()) norm = roundUp(scalar.multiplyRounded(norm));
                const auto kappa = deterministicBound(norm,
                    "v9 3.2/3.5: product of radix-2 row sums, unit phase and first-factor normalization");
                const auto delta = unavailable("Rigorous ideal-twiddle/composition and SEAL diagonal encoding bound unavailable");
                const auto local = unavailable("Double-hoisted BSGS decomposition/key noise, inner/final ModDown, "
                    "rescale and scale-representation bounds unavailable");
                outputs[half] = applyCipherMatrix(adapter,
                    half == 1 && index == 0 ? constants.pimpl_->secondFirstFactor : constants.pimpl_->factors[index], outputs[half]);
                const auto nextMagnitude = magnitudeStep(kappa, magnitude);
                const auto nextError = errorStep(kappa, delta, magnitude, error, local);
                trace.halves[half].push_back({"bsgs_factor_rescale", adapter.info(outputs[half]),
                    adapter.coeffModulusValues(outputs[half]), kappa, delta, nextMagnitude, nextError, local});
                magnitude = nextMagnitude;
                error = nextError;
                ++trace.rescaleOperations;
            }
            const auto kappa = deterministicBound(2, "v9 (4): norm of x + conjugate(x) is at most two");
            const auto delta = deterministicBound(0, "Exact real projection has no encoded diagonal");
            const auto local = unavailable("Rigorous conjugation key-switch and addition local bounds unavailable");
            outputs[half] = adapter.add(outputs[half], adapter.conjugate(outputs[half]));
            trace.halves[half].push_back({"conjugate_add", adapter.info(outputs[half]),
                adapter.coeffModulusValues(outputs[half]), kappa, delta, magnitudeStep(kappa, magnitude),
                errorStep(kappa, delta, magnitude, error, local), local});
            if (std::abs(std::log2(adapter.scale(outputs[half])) - contract.outputScaleLog2)
                > contract.inputScaleToleranceLog2)
                throw std::runtime_error("EvalRound+ branch output scale violates contract");
        }
        trace.levelsConsumed = start.chainIndex - adapter.chainIndex(outputs[0]);
        trace.outputError = unavailable("v9 (6): branch bound requires rigorous diagonal encoding, "
            "double-hoisted BSGS and final conjugation bounds; CKKS measurements are diagnostics only");
        trace.evidence = {false, trace.outputError.provenance, {trace.outputError}};
        trace.certificate = {BootstrapCertificationStatus::RequiredBoundUnavailable,
            gate == BootstrapGate::CoeffToSlotHP ? "coeff_to_slot_hp" : "coeff_to_slot_lp", trace.outputError.provenance};
        return std::make_pair(std::move(outputs), std::move(trace));
    };
    auto hp = branch(prepared.hp_, prepared.hpContract_, BootstrapGate::CoeffToSlotHP, {});
    auto lp = branch(prepared.lp_, prepared.lpContract_, BootstrapGate::CoeffToSlotLP, alpha);
    return {std::move(hp.first[0]), std::move(hp.first[1]), std::move(lp.first[0]), std::move(lp.first[1]),
            std::move(hp.second), std::move(lp.second)};
}

std::vector<std::size_t> EvalRoundPlusCoeffToSlot::certificationBabySteps() const {
    return plan_.pimpl_->babySteps;
}

} // namespace m2424
