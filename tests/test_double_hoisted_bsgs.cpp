#include "m2424/seal_adapter.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <functional>
#include <random>
#include <utility>
#include <vector>

namespace {

struct OwnedTerm { int baby{}; m2424::Plain ordinary; m2424::Plain extended; };
struct OwnedGroup { int giant{}; std::vector<OwnedTerm> terms; };

double maxError(const std::vector<std::complex<double>>& expected,
                const std::vector<std::complex<double>>& actual) {
    double result = 0.0;
    for (std::size_t index = 0; index < expected.size(); ++index) {
        result = std::max(result, std::abs(expected[index] - actual[index]));
    }
    return result;
}

bool runProfile(const m2424::CkksProfile& profile, unsigned seed,
                std::size_t trials, double errorLimit) {
    auto adapter = m2424::SealAdapter::create(profile);
    adapter.generateKeys(std::vector<int>{-1, 1, 2}, false);
    std::mt19937_64 random(seed);
    std::uniform_real_distribution<double> distribution(-0.25, 0.25);
    double worstError = 0.0;

    for (std::size_t levelDrop = 0; levelDrop < 2; ++levelDrop) {
        for (std::size_t trial = 0; trial < trials; ++trial) {
            std::vector<std::complex<double>> values(profile.slots);
            for (auto& value : values) value = {distribution(random), distribution(random)};
            auto input = adapter.encrypt(adapter.encodeComplex(values));
            for (std::size_t level = 0; level < levelDrop; ++level) {
                input = adapter.rescaleToNext(adapter.multiplyPlain(
                    input, adapter.encodeScalarAtScaleFor(
                        1.0, std::exp2(static_cast<double>(
                            profile.coeffModulusBits[profile.coeffModulusBits.size() - 2 - level])), input)));
            }

            std::vector<OwnedGroup> owned;
            const double diagonalScale = profile.scale;
            for (const auto& specification :
                 std::vector<std::pair<int, std::vector<int>>>{{0, {0, 1}}, {2, {0, -1}}}) {
                OwnedGroup group;
                group.giant = specification.first;
                for (int baby : specification.second) {
                    std::vector<std::complex<double>> diagonal(profile.slots);
                    for (auto& value : diagonal) value = {distribution(random), distribution(random)};
                    group.terms.push_back({
                        baby,
                        adapter.encodeComplexAtScaleFor(diagonal, diagonalScale, input),
                        adapter.encodeComplexAtKeyScale(diagonal, diagonalScale)});
                }
                owned.push_back(std::move(group));
            }

            std::vector<m2424::HoistedBsgsGroup> backendGroups;
            m2424::Cipher ordinary;
            bool ordinaryInitialized = false;
            for (auto& group : owned) {
                m2424::HoistedBsgsGroup backendGroup;
                backendGroup.giantRotation = group.giant;
                m2424::Cipher inner;
                bool innerInitialized = false;
                for (auto& term : group.terms) {
                    backendGroup.terms.push_back({term.baby, &term.extended});
                    const auto rotated = term.baby == 0 ? input : adapter.rotate(input, term.baby);
                    auto weighted = adapter.multiplyPlain(rotated, term.ordinary);
                    inner = innerInitialized ? adapter.add(inner, weighted) : std::move(weighted);
                    innerInitialized = true;
                }
                auto shifted = group.giant == 0 ? inner : adapter.rotate(inner, group.giant);
                ordinary = ordinaryInitialized ? adapter.add(ordinary, shifted) : std::move(shifted);
                ordinaryInitialized = true;
                backendGroups.push_back(std::move(backendGroup));
            }

            const auto hoisted = adapter.applyBsgsDoubleHoisted(input, backendGroups);
            const auto expected = adapter.decodeComplex(adapter.decrypt(ordinary));
            const auto actual = adapter.decodeComplex(adapter.decrypt(hoisted));
            worstError = std::max(worstError, maxError(expected, actual));
            if (adapter.info(ordinary).chainIndex != adapter.info(hoisted).chainIndex ||
                std::abs(std::log2(adapter.scale(ordinary)) - std::log2(adapter.scale(hoisted))) > 1e-9) {
                return false;
            }
        }
    }
    std::printf("[test_double_hoisted_bsgs] N=%zu error=%.3e\n",
                profile.polyModulusDegree, worstError);
    return worstError <= errorLimit;
}

bool rejectsMalformedGroups() {
    auto adapter = m2424::SealAdapter::create(
        {4096, {35, 30, 35}, std::exp2(20.0), 2048});
    adapter.generateKeys(std::vector<int>{1}, false);
    std::vector<std::complex<double>> values(2048, {0.125, -0.0625});
    const auto input = adapter.encrypt(adapter.encodeComplex(values));
    auto first = adapter.encodeComplexAtKeyScale(values, std::exp2(20.0));
    auto second = adapter.encodeComplexAtKeyScale(values, std::exp2(19.0));
    auto ordinary = adapter.encodeComplexAtScaleFor(values, std::exp2(20.0), input);
    const auto rejected = [&](std::vector<m2424::HoistedBsgsGroup> groups) {
        try {
            (void)adapter.applyBsgsDoubleHoisted(input, groups);
            return false;
        } catch (const std::invalid_argument&) {
            return true;
        }
    };

    std::vector<m2424::HoistedBsgsTerm> tooMany(257, {0, &first});
    return rejected({})
        && rejected({{0, {}}})
        && rejected({{0, {{0, &first}, {0, &first}}}})
        && rejected({{0, {{0, &first}, {1, &second}}}})
        && rejected({{0, {{0, &first}}}, {0, {{1, &first}}}})
        && rejected({{0, std::move(tooMany)}})
        && rejected({{0, {{0, &ordinary}}}});
}

} // namespace

int main() {
    const bool medium = runProfile(
        {8192, {50, 35, 35, 35, 50}, std::exp2(30.0), 4096}, 11, 5, 1e-3);
    const bool large = runProfile(
        {16384, {60, 60, 60, 60, 60, 60, 60}, std::exp2(59.5), 8192}, 17, 3, 2e-8);
    return medium && large && rejectsMalformedGroups() ? 0 : 1;
}
