#include "m2424/accuracy.hpp"
#include "m2424/bootstrap_candidates.hpp"
#include "m2424/coeff_to_slot_bsgs_plan.hpp"
#include "m2424/coeff_to_slot_contract.hpp"
#include "m2424/coeff_to_slot_fft_plan.hpp"
#include "m2424/seal_adapter.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;
constexpr std::size_t kTrials = 7;
constexpr double kPi = 3.141592653589793238462643383279502884;

struct Measurement {
    double keySetupMs{};
    double prepareMs{};
    double minApplyMs{};
    double medianApplyMs{};
    double error{};
    std::size_t levels{};
    std::size_t rotationKeyCount{};
    std::size_t galoisKeyBytes{};
    std::size_t preparedPlaintextBytes{};
    bool passes{};
};

template <class Function>
double elapsedMs(Function&& function) {
    const auto start = Clock::now();
    function();
    const auto finish = Clock::now();
    return std::chrono::duration<double, std::milli>(finish - start).count();
}

m2424::CkksProfile comparisonProfile() {
    return {16384, {60, 59, 59, 59, 59, 59, 60}, std::exp2(59.0), 8192};
}

m2424::ComplexVector makeInput(std::size_t size) {
    m2424::ComplexVector result;
    result.reserve(size);
    for (std::size_t index = 0; index < size; ++index) {
        const double phase = static_cast<double>(index) * 0.6180339887498948;
        result.push_back({0.75 * std::sin(phase), 0.75 * std::cos(phase * 1.4142135623730950)});
    }
    return result;
}

m2424::ComplexVector directDft(const m2424::ComplexVector& input) {
    m2424::ComplexVector result(input.size());
    for (std::size_t row = 0; row < input.size(); ++row) {
        for (std::size_t column = 0; column < input.size(); ++column) {
            const double angle = -2.0 * kPi * static_cast<double>(row * column)
                / static_cast<double>(input.size());
            result[row] += input[column] * m2424::Complex{std::cos(angle), std::sin(angle)};
        }
    }
    return result;
}

double maxError(const m2424::ComplexVector& expected, const m2424::ComplexVector& actual) {
    double result = 0.0;
    for (std::size_t index = 0; index < expected.size(); ++index) {
        result = std::max(result, std::abs(expected[index] - actual[index]));
    }
    return result;
}

double median(std::vector<double> values) {
    std::sort(values.begin(), values.end());
    return values[values.size() / 2];
}

template <class Plan>
Measurement measurePlan(Plan& plan, const m2424::ComplexVector& input) {
    const auto& candidate = m2424::bootstrapCandidateById("precision_8192_s59");
    auto adapter = m2424::SealAdapter::create(comparisonProfile());
    Measurement measurement;
    measurement.rotationKeyCount = plan.rotationSteps().size();
    measurement.levels = plan.requiredLevels();
    measurement.keySetupMs = elapsedMs([&] { adapter.generateKeys(plan.rotationSteps(), false); });
    measurement.galoisKeyBytes = adapter.galoisKeysSize();

    const auto encrypted = adapter.encrypt(adapter.encodeComplex(input));
    const auto contract = m2424::makeCoeffToSlotContract(candidate);
    const auto preflight = m2424::preflightCoeffToSlot(adapter, encrypted, contract, plan.requirements());

    measurement.prepareMs = elapsedMs([&] { plan.prepare(adapter, encrypted); });
    measurement.preparedPlaintextBytes = plan.preparedPlaintextBytes(adapter);

    const auto expected = directDft(input);
    const auto validationResult = plan.apply(adapter, encrypted);
    const auto decoded = adapter.decodeComplex(adapter.decrypt(validationResult));
    const m2424::ComplexVector actual(decoded.begin(), decoded.begin() + static_cast<std::ptrdiff_t>(input.size()));
    measurement.error = maxError(expected, actual);
    const auto resultInfo = adapter.info(validationResult);
    const auto inputInfo = adapter.info(encrypted);

    std::vector<double> samples;
    samples.reserve(kTrials);
    for (std::size_t trial = 0; trial < kTrials; ++trial) {
        samples.push_back(elapsedMs([&] {
            const auto result = plan.apply(adapter, encrypted);
            (void)result;
        }));
    }
    measurement.minApplyMs = *std::min_element(samples.begin(), samples.end());
    measurement.medianApplyMs = median(std::move(samples));
    measurement.passes = preflight.ready
        && measurement.error <= candidate.errorBudget.coeffToSlot
        && resultInfo.chainIndex + measurement.levels == inputInfo.chainIndex
        && std::abs(std::log2(resultInfo.scale) - std::log2(inputInfo.scale)) <= 0.25;
    return measurement;
}

void printRow(const char* strategy, std::size_t size, std::size_t babyStep, const Measurement& measurement) {
    std::printf("precision_8192_s59,%s,%zu,%zu,%zu,%.3f,%.3f,%.3f,%.3f,%.3e,%zu,%zu,%zu,%zu,%s\n",
                strategy, size, babyStep, kTrials, measurement.keySetupMs, measurement.prepareMs,
                measurement.minApplyMs, measurement.medianApplyMs, measurement.error, measurement.levels,
                measurement.rotationKeyCount, measurement.galoisKeyBytes, measurement.preparedPlaintextBytes,
                measurement.passes ? "PASS" : "FAIL");
}

} // namespace

int main() {
    std::printf("candidate,strategy,transform_slots,baby_step,trials,key_setup_ms,prepare_ms,min_apply_ms,"
                "median_apply_ms,max_abs_error,levels,rotation_key_count,galois_key_bytes,"
                "prepared_plaintext_bytes,status\n");
    bool ok = true;
    for (const std::size_t size : {std::size_t{4}, std::size_t{8}, std::size_t{16}}) {
        const auto input = makeInput(size);
        m2424::CoeffToSlotFftPlan fft(size, 8192);
        const auto fftMeasurement = measurePlan(fft, input);
        printRow("fft", size, 0, fftMeasurement);
        ok = fftMeasurement.passes && ok;

        for (std::size_t babyStep = 1; babyStep <= size; babyStep *= 2) {
            m2424::CoeffToSlotBsgsPlan bsgs(size, 8192, babyStep);
            const auto bsgsMeasurement = measurePlan(bsgs, input);
            printRow("bsgs", size, babyStep, bsgsMeasurement);
            ok = bsgsMeasurement.passes && ok;
        }
    }
    return ok ? 0 : 1;
}
