#include "m2424/accuracy.hpp"
#include "m2424/coeff_to_slot_fft_plan.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

m2424::ComplexVector directDft(const m2424::ComplexVector& input) {
    m2424::ComplexVector output;
    output.reserve(input.size());
    for (std::size_t row = 0; row < input.size(); ++row) {
        m2424::Complex value{};
        for (std::size_t column = 0; column < input.size(); ++column) {
            const double angle = -2.0 * kPi * static_cast<double>(row * column)
                / static_cast<double>(input.size());
            value += input[column] * m2424::Complex{std::cos(angle), std::sin(angle)};
        }
        output.push_back(value);
    }
    return output;
}

double maxError(const m2424::ComplexVector& expected, const m2424::ComplexVector& actual) {
    double result = 0.0;
    for (std::size_t index = 0; index < expected.size(); ++index) {
        result = std::max(result, std::abs(expected[index] - actual[index]));
    }
    return result;
}

m2424::CkksProfile fftCalibrationProfile() {
    return {16384, {60, 59, 59, 59, 59, 59, 60}, std::exp2(59.0), 8192};
}

} // namespace

int main() {
    const m2424::ComplexVector input{
        {0.25, -0.5}, {-0.75, 0.125}, {0.5, 0.25}, {-0.125, -0.375},
        {0.625, 0.75}, {-0.25, 0.5}, {0.875, -0.125}, {-0.5, -0.75},
        {0.125, 0.375}, {-0.625, 0.25}, {0.75, -0.5}, {-0.375, 0.625},
        {0.5, 0.0}, {-0.25, -0.25}, {0.125, -0.625}, {-0.75, 0.75}
    };
    bool plainOk = true;
    for (const std::size_t size : {std::size_t{4}, std::size_t{8}, std::size_t{16}}) {
        const m2424::CoeffToSlotFftPlan plan(size, 8192);
        const m2424::ComplexVector sample(input.begin(), input.begin() + static_cast<std::ptrdiff_t>(size));
        plainOk = plainOk && maxError(directDft(sample), plan.applyPlain(sample)) <= 1e-12;
    }

    auto adapter = m2424::SealAdapter::create(fftCalibrationProfile());
    const m2424::CoeffToSlotFftPlan encryptedPlan(16, adapter.slotCount());
    adapter.generateKeys(encryptedPlan.rotationSteps(), false);
    const auto encrypted = adapter.encrypt(adapter.encodeComplex(input));
    const auto preflight = m2424::preflightCoeffToSlot(
        adapter,
        encrypted,
        m2424::makeCoeffToSlotContract(m2424::bootstrapCandidateById("precision_8192_s59")),
        encryptedPlan.requirements());
    encryptedPlan.prepare(adapter, encrypted);
    const std::size_t preparedBytes = encryptedPlan.preparedPlaintextBytes(adapter);
    const auto result = encryptedPlan.apply(adapter, encrypted);
    const auto decoded = adapter.decodeComplex(adapter.decrypt(result));
    const auto expected = directDft(input);
    const m2424::ComplexVector actual(decoded.begin(), decoded.begin() + static_cast<std::ptrdiff_t>(input.size()));
    const double encryptedError = maxError(expected, actual);
    const auto inputInfo = adapter.info(encrypted);
    const auto resultInfo = adapter.info(result);

    const bool encryptedOk = preflight.ready
        && preparedBytes > 0
        && encryptedError <= m2424::bootstrapCandidateById("precision_8192_s59").errorBudget.coeffToSlot
        && resultInfo.chainIndex + encryptedPlan.requiredLevels() == inputInfo.chainIndex
        && std::abs(std::log2(resultInfo.scale) - std::log2(inputInfo.scale)) <= 0.25;
    const bool ok = plainOk && encryptedOk;
    std::printf("[test_coeff_to_slot_fft_plan] plain=%s encrypted_error=%.3e levels=%zu prepared_bytes=%zu %s\n",
                plainOk ? "PASS" : "FAIL", encryptedError, encryptedPlan.requiredLevels(), preparedBytes,
                ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
