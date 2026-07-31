#include "m2424/coeff_to_slot.hpp"

#include <algorithm>
#include <cmath>
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
        adapter.applyBsgsInnerDoubleHoisted(input, groups));
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
            factorization.radices = balancedRadices(raw.size(), targetDepth);
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
        result.hoistedDecompositionsPerApply += babies.empty() ? 0 : 2;
        result.hoistedAutomorphismsPerApply += 2 * babies.size();
        result.ordinaryGiantRotationsPerApply += 2 * giants.size();
        result.additionsPerApply += 2 * (diagonals - 1);
        result.storedComplexValues += diagonals * pimpl_->slots;
    }
    result.additionsPerApply += 2;
    result.rescalesPerApply = 2 * depth();
    result.uniqueEvaluationKeys = requirements().rotationSteps.size() + 1;
    result.storedComplexValues += pimpl_->secondFirstFactor.size() * pimpl_->slots;
    return result;
}

std::vector<CoeffToSlotFactorizationEstimate> CoeffToSlotPlan::rankFactorizations(
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
    if (slots.size() != pimpl_->slots) {
        throw std::invalid_argument("CoeffToSlotPlan plaintext slot count mismatch");
    }
    ComplexVector first = applyMatrix(pimpl_->factors.front(), slots);
    ComplexVector second = applyMatrix(pimpl_->secondFirstFactor, slots);
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
            pimpl_->factors[index], pimpl_->babySteps[index], info.chainIndex - index));
    }
    preparedSecond = encodeFactor(
        pimpl_->secondFirstFactor, pimpl_->babySteps.front(), info.chainIndex);

    auto result = std::make_unique<PreparedCoeffToSlotPlan::Impl>();
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
    if (!prepared.pimpl_) {
        return false;
    }
    const auto& state = *prepared.pimpl_;
    if (state.contextFingerprint != adapter.contextFingerprint()) {
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

} // namespace m2424
