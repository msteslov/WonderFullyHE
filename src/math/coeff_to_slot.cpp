#include "m2424/coeff_to_slot.hpp"

#include "m2424/canonical_embedding_reference.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <utility>

namespace m2424 {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

std::size_t validatedDegree(std::size_t degree) {
    if (degree < 2 || (degree & (degree - 1)) != 0) {
        throw std::invalid_argument("CoeffToSlot polyModulusDegree must be a power of two");
    }
    return degree;
}

std::size_t wrappedIndex(std::size_t value, std::size_t modulus) {
    return value % modulus;
}

struct BsgsTerm {
    int babyRotation{};
};

struct BsgsGroup {
    int giantRotation{};
    std::vector<BsgsTerm> terms;
};

class CanonicalBsgsTransform {
public:
    CanonicalBsgsTransform(std::size_t degree, std::size_t coefficientOffset)
        : degree_(degree),
          coefficientOffset_(coefficientOffset),
          slots_(degree / 2),
          rootOrder_(2 * degree),
          exponents_(canonicalEmbeddingRootExponents(degree)),
          babyStep_(static_cast<std::size_t>(std::ceil(std::sqrt(static_cast<double>(slots_))))) {
        std::map<int, BsgsGroup> grouped;
        for (std::size_t rotation = 0; rotation < slots_; ++rotation) {
            const int giant = static_cast<int>((rotation / babyStep_) * babyStep_);
            const int baby = static_cast<int>(rotation) - giant;
            auto& group = grouped[giant];
            group.giantRotation = giant;
            BsgsTerm term;
            term.babyRotation = baby;
            group.terms.push_back(std::move(term));
        }
        groups_.reserve(grouped.size());
        for (auto& entry : grouped) {
            groups_.push_back(std::move(entry.second));
        }
    }

    Cipher apply(SealAdapter& adapter, const Cipher& input) const {
        const auto info = adapter.info(input);
        const auto bits = adapter.coeffModulusBits();
        if (info.coeffModulusSize == 0 || info.coeffModulusSize > bits.size()) {
            throw std::runtime_error("CoeffToSlot has no rescale modulus");
        }
        const double plaintextScale = std::exp2(
            static_cast<double>(bits[info.coeffModulusSize - 1]));
        std::map<int, Cipher> babies;
        babies.emplace(0, input);
        for (std::size_t baby = 1; baby < babyStep_; ++baby) {
            babies.emplace(static_cast<int>(baby), adapter.rotate(input, static_cast<int>(baby)));
        }

        Cipher result;
        bool hasResult = false;
        for (std::size_t groupIndex = 0; groupIndex < groups_.size(); ++groupIndex) {
            const auto& group = groups_[groupIndex];
            Cipher inner;
            bool hasInner = false;
            for (std::size_t termIndex = 0; termIndex < group.terms.size(); ++termIndex) {
                const auto& term = group.terms[termIndex];
                const std::size_t rotation = static_cast<std::size_t>(
                    group.giantRotation + term.babyRotation);
                ComplexVector diagonal(slots_);
                for (std::size_t row = 0; row < slots_; ++row) {
                    const std::size_t column = (row + rotation) % slots_;
                    const std::size_t coefficient = row + coefficientOffset_;
                    const std::size_t exponent =
                        (rootOrder_ - ((exponents_[column] * coefficient) % rootOrder_)) % rootOrder_;
                    const double angle = 2.0 * kPi * static_cast<double>(exponent)
                        / static_cast<double>(rootOrder_);
                    diagonal[wrappedIndex(
                        row + static_cast<std::size_t>(group.giantRotation), slots_)] =
                        Complex{std::cos(angle), std::sin(angle)} / static_cast<double>(degree_);
                }
                const Plain plaintext =
                    adapter.encodeComplexAtScaleFor(diagonal, plaintextScale, input);
                Cipher weighted = adapter.rescaleToNext(adapter.multiplyPlain(
                    babies.at(term.babyRotation), plaintext));
                if (!hasInner) {
                    inner = std::move(weighted);
                    hasInner = true;
                } else {
                    inner = adapter.add(inner, weighted);
                }
            }
            Cipher shifted = group.giantRotation == 0
                ? std::move(inner)
                : adapter.rotate(inner, group.giantRotation);
            if (!hasResult) {
                result = std::move(shifted);
                hasResult = true;
            } else {
                result = adapter.add(result, shifted);
            }
        }
        return result;
    }

private:
    std::size_t degree_{};
    std::size_t coefficientOffset_{};
    std::size_t slots_{};
    std::size_t rootOrder_{};
    std::vector<std::size_t> exponents_;
    std::size_t babyStep_{};
    std::vector<BsgsGroup> groups_;
};

} // namespace

struct CoeffToSlot::Impl {
    explicit Impl(std::size_t degree)
        : degree(validatedDegree(degree)),
          slots(degree / 2),
          babyStep(static_cast<std::size_t>(std::ceil(std::sqrt(static_cast<double>(slots))))),
          first(degree, 0),
          second(degree, slots) {}

    std::size_t degree{};
    std::size_t slots{};
    std::size_t babyStep{};
    CanonicalBsgsTransform first;
    CanonicalBsgsTransform second;
};

CoeffToSlot::CoeffToSlot(std::size_t polyModulusDegree)
    : pimpl_(std::make_unique<Impl>(polyModulusDegree)) {}
CoeffToSlot::~CoeffToSlot() = default;
CoeffToSlot::CoeffToSlot(CoeffToSlot&&) noexcept = default;
CoeffToSlot& CoeffToSlot::operator=(CoeffToSlot&&) noexcept = default;

CoeffToSlotPlanRequirements CoeffToSlot::requirements() const {
    std::vector<int> rotations;
    for (std::size_t step = 1; step < pimpl_->babyStep; ++step) {
        rotations.push_back(static_cast<int>(step));
    }
    for (std::size_t step = pimpl_->babyStep; step < pimpl_->slots; step += pimpl_->babyStep) {
        rotations.push_back(static_cast<int>(step));
    }
    return {1, std::move(rotations), true};
}

CoeffToSlotResult CoeffToSlot::apply(SealAdapter& adapter, RaisedCipher&& input) const {
    if (adapter.slotCount() != pimpl_->slots) {
        throw std::invalid_argument("CoeffToSlot degree does not match RaisedCipher context");
    }
    Cipher coefficientCipher = std::move(input.cipher_);
    Cipher first = pimpl_->first.apply(adapter, coefficientCipher);
    Cipher second = pimpl_->second.apply(adapter, coefficientCipher);
    first = adapter.add(first, adapter.conjugate(first));
    second = adapter.add(second, adapter.conjugate(second));
    return {std::move(first), std::move(second)};
}

} // namespace m2424
