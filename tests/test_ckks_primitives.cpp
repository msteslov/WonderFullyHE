#include "m2424/accuracy.hpp"
#include "m2424/profiles.hpp"
#include "m2424/seal_adapter.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

namespace {

std::vector<double> head(const std::vector<double>& values, std::size_t size) {
    return {values.begin(), values.begin() + static_cast<std::ptrdiff_t>(size)};
}

} // namespace

int main() {
    auto adapter = m2424::SealAdapter::create(m2424::profiles::high_precision_ckks());
    adapter.generateKeys(std::vector<int>{1, 2, 3}, true);

    const std::vector<double> input{0.125, -0.25, 0.375, -0.5};
    const auto encrypted = adapter.encrypt(adapter.encode(input));

    std::vector<double> doubled;
    std::vector<double> squared;
    std::vector<double> rotated;

    const auto added = adapter.add(encrypted, encrypted);
    doubled = head(adapter.decode(adapter.decrypt(added)), input.size());

    const auto multiplied = adapter.rescaleToNext(adapter.relinearize(adapter.multiply(encrypted, encrypted)));
    squared = head(adapter.decode(adapter.decrypt(multiplied)), input.size());
    const auto dataPrimes = adapter.dataModulusValues();
    const auto inputPrimes = adapter.coeffModulusValues(encrypted);
    const auto multipliedPrimes = adapter.coeffModulusValues(multiplied);
    const auto specialPrime = adapter.specialKeyModulusValue();

    const auto rotation = adapter.rotate(encrypted, 1);
    rotated = head(adapter.decode(adapter.decrypt(rotation)), input.size());
    const std::vector<int> hoistedSteps{1, 2, 3};
    const auto hoisted = adapter.rotateManyHoisted(encrypted, hoistedSteps);

    const std::vector<double> doubledExpected{0.25, -0.5, 0.75, -1.0};
    const std::vector<double> squaredExpected{0.015625, 0.0625, 0.140625, 0.25};
    const std::vector<double> rotatedExpected{-0.25, 0.375, -0.5, 0.0};

    const auto addAccuracy = m2424::compare(doubledExpected, doubled, m2424::kTargetAbsoluteError);
    const auto multiplyAccuracy = m2424::compare(squaredExpected, squared, m2424::kTargetAbsoluteError);
    const auto rotateAccuracy = m2424::compare(rotatedExpected, rotated, m2424::kTargetAbsoluteError);
    double hoistedError = 0.0;
    bool hoistedOk = hoisted.size() == hoistedSteps.size();
    for (std::size_t index = 0; index < hoisted.size(); ++index) {
        const auto expected = adapter.decode(adapter.decrypt(adapter.rotate(encrypted, hoistedSteps[index])));
        const auto actual = adapter.decode(adapter.decrypt(hoisted[index]));
        const auto accuracy = m2424::compare(expected, actual, m2424::kTargetAbsoluteError);
        hoistedError = std::max(hoistedError, accuracy.max_abs_error);
        hoistedOk = hoistedOk && accuracy.ok &&
            adapter.info(hoisted[index]).chainIndex == adapter.info(encrypted).chainIndex &&
            adapter.scale(hoisted[index]) == adapter.scale(encrypted);
    }
    const bool exactModuliOk = inputPrimes == dataPrimes && dataPrimes.size() >= 2
        && multipliedPrimes.size() + 1 == inputPrimes.size()
        && std::equal(multipliedPrimes.begin(), multipliedPrimes.end(), inputPrimes.begin())
        && specialPrime != 0
        && std::find(dataPrimes.begin(), dataPrimes.end(), specialPrime) == dataPrimes.end();
    const bool ok = addAccuracy.ok && multiplyAccuracy.ok && rotateAccuracy.ok && hoistedOk
        && exactModuliOk;

    std::printf("[test_ckks_primitives] target=%.1e add=%.3e multiply=%.3e rotate=%.3e hoisted=%.3e %s\n",
                m2424::kTargetAbsoluteError,
                addAccuracy.max_abs_error,
                multiplyAccuracy.max_abs_error,
                rotateAccuracy.max_abs_error,
                hoistedError,
                ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
