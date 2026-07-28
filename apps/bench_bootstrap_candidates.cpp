#include "m2424/accuracy.hpp"
#include "m2424/bootstrap_candidates.hpp"
#include "m2424/seal_adapter.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <exception>
#include <string>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;
constexpr std::size_t kCalibrationTrials = 5;

struct PrimitiveMeasurement {
    double setupMs{};
    double addError{};
    double multiplyError{};
    double rotationError{};
    double maxRotationMs{};
    std::size_t inputChainIndex{};
    std::size_t multiplyChainIndex{};
    std::size_t publicKeyBytes{};
    std::size_t relinKeyBytes{};
    std::size_t galoisKeyBytes{};
    bool keysWithinLimit{};
    bool precisionPasses{};
};

template <class Fn>
double elapsedMs(Fn&& fn) {
    const auto start = Clock::now();
    fn();
    const auto finish = Clock::now();
    return std::chrono::duration<double, std::milli>(finish - start).count();
}

m2424::CkksProfile makeCalibrationProfile(const m2424::BootstrapCandidate& candidate) {
    const int scaleBits = static_cast<int>(std::ceil(candidate.targetScaleLog2));
    const std::size_t workPrimes = candidate.polyModulusDegree >= 16384 ? 2 : 1;
    std::vector<int> coeffModulusBits{60};
    coeffModulusBits.insert(coeffModulusBits.end(), workPrimes, scaleBits);
    coeffModulusBits.push_back(60);
    return {candidate.polyModulusDegree, coeffModulusBits, std::exp2(candidate.targetScaleLog2), candidate.activeSlots};
}

std::vector<int> calibrationRotationSteps(std::size_t activeSlots) {
    std::vector<int> steps;
    for (std::size_t step = 1; step <= activeSlots / 2; step *= 2) {
        steps.push_back(static_cast<int>(step));
        if (step > activeSlots / 4) {
            break;
        }
    }
    return steps;
}

std::vector<double> makeInput(std::size_t activeSlots, std::size_t trial) {
    std::vector<double> input;
    input.reserve(activeSlots);
    for (std::size_t index = 0; index < activeSlots; ++index) {
        const double phase = static_cast<double>(index) * 0.6180339887498948
            + static_cast<double>(trial) * 0.4142135623730950;
        input.push_back(std::sin(phase));
    }
    if (!input.empty()) {
        input.front() = -1.0;
        input.back() = 1.0;
    }
    return input;
}

std::vector<double> head(const std::vector<double>& values, std::size_t size) {
    return {values.begin(), values.begin() + static_cast<std::ptrdiff_t>(size)};
}

std::vector<double> rotateReference(const std::vector<double>& input,
                                    std::size_t physicalSlots,
                                    int steps) {
    std::vector<double> padded(physicalSlots, 0.0);
    std::copy(input.begin(), input.end(), padded.begin());
    std::vector<double> expected;
    expected.reserve(input.size());
    for (std::size_t index = 0; index < input.size(); ++index) {
        expected.push_back(padded[(index + static_cast<std::size_t>(steps)) % physicalSlots]);
    }
    return expected;
}

PrimitiveMeasurement measure(const m2424::BootstrapCandidate& candidate) {
    const auto profile = makeCalibrationProfile(candidate);
    const auto steps = calibrationRotationSteps(candidate.activeSlots);

    m2424::SealAdapter adapter;
    PrimitiveMeasurement result;
    result.setupMs = elapsedMs([&] {
        adapter = m2424::SealAdapter::create(profile);
        adapter.generateKeys(steps, true);
    });

    for (std::size_t trial = 0; trial < kCalibrationTrials; ++trial) {
        const auto input = makeInput(candidate.activeSlots, trial);
        const auto encrypted = adapter.encrypt(adapter.encode(input));
        result.inputChainIndex = adapter.chainIndex(encrypted);
        const auto added = adapter.add(encrypted, encrypted);
        const auto multiplied = adapter.rescaleToNext(adapter.relinearize(adapter.multiply(encrypted, encrypted)));
        result.multiplyChainIndex = adapter.chainIndex(multiplied);

        std::vector<double> doubled;
        std::vector<double> squared;
        doubled.reserve(input.size());
        squared.reserve(input.size());
        for (double value : input) {
            doubled.push_back(2.0 * value);
            squared.push_back(value * value);
        }
        result.addError = std::max(result.addError, m2424::compare(
            doubled, head(adapter.decode(adapter.decrypt(added)), input.size()),
            m2424::kTargetAbsoluteError).max_abs_error);
        result.multiplyError = std::max(result.multiplyError, m2424::compare(
            squared, head(adapter.decode(adapter.decrypt(multiplied)), input.size()),
            m2424::kTargetAbsoluteError).max_abs_error);

        for (const int step : steps) {
            m2424::Cipher rotated;
            const double rotationMs = elapsedMs([&] { rotated = adapter.rotate(encrypted, step); });
            const auto expected = rotateReference(input, adapter.slotCount(), step);
            const auto actual = head(adapter.decode(adapter.decrypt(rotated)), input.size());
            result.rotationError = std::max(result.rotationError,
                                            m2424::compare(expected, actual, m2424::kTargetAbsoluteError).max_abs_error);
            result.maxRotationMs = std::max(result.maxRotationMs, rotationMs);
        }
    }

    result.publicKeyBytes = adapter.publicKeySize();
    result.relinKeyBytes = adapter.relinKeysSize();
    result.galoisKeyBytes = adapter.galoisKeysSize();
    const double evaluationKeyMiB = static_cast<double>(result.relinKeyBytes + result.galoisKeyBytes)
        / (1024.0 * 1024.0);
    result.keysWithinLimit = evaluationKeyMiB <= static_cast<double>(candidate.resources.maxEvaluationKeyMiB);
    result.precisionPasses = result.addError <= m2424::kTargetAbsoluteError
        && result.multiplyError <= m2424::kTargetAbsoluteError
        && result.rotationError <= m2424::kTargetAbsoluteError;
    return result;
}

void printMeasurement(const m2424::BootstrapCandidate& candidate, const PrimitiveMeasurement& measurement) {
    std::printf("%s,%zu,%zu,%.0f,%zu,%.3f,%.3e,%.3e,%.3e,%.3f,%zu,%zu,%zu,%zu,%zu,%s,%s\n",
                candidate.id.c_str(),
                candidate.activeSlots,
                candidate.polyModulusDegree,
                candidate.targetScaleLog2,
                kCalibrationTrials,
                measurement.setupMs,
                measurement.addError,
                measurement.multiplyError,
                measurement.rotationError,
                measurement.maxRotationMs,
                measurement.inputChainIndex,
                measurement.multiplyChainIndex,
                measurement.publicKeyBytes,
                measurement.relinKeyBytes,
                measurement.galoisKeyBytes,
                measurement.precisionPasses ? "PASS" : "FAIL",
                measurement.keysWithinLimit ? "PASS" : "FAIL");
}

bool measureCandidate(const m2424::BootstrapCandidate& candidate) {
    try {
        const auto measurement = measure(candidate);
        printMeasurement(candidate, measurement);
        return measurement.precisionPasses && measurement.keysWithinLimit;
    } catch (const std::exception& error) {
        std::printf("%s,%zu,%zu,%.0f,%zu,0,0,0,0,0,0,0,0,0,0,ERROR,%s\n",
                    candidate.id.c_str(), candidate.activeSlots, candidate.polyModulusDegree,
                    candidate.targetScaleLog2, kCalibrationTrials, error.what());
        return false;
    }
}

} // namespace

int main(int argc, char** argv) {
    std::printf("candidate,active_slots,poly_modulus_degree,scale_log2,trials,key_setup_ms,add_max_abs_error,"
                "multiply_max_abs_error,rotate_max_abs_error,max_rotate_ms,input_chain_index,"
                "multiply_chain_index,public_key_bytes,relin_key_bytes,galois_key_bytes,precision,status_keys\n");

    bool ok = true;
    if (argc == 2 && std::string(argv[1]) != "--all") {
        return measureCandidate(m2424::bootstrapCandidateById(argv[1])) ? 0 : 1;
    }
    if (argc > 2) {
        std::fprintf(stderr, "usage: bench_bootstrap_candidates [candidate_id|--all]\n");
        return 2;
    }
    for (const auto& candidate : m2424::bootstrapCandidates()) {
        ok = measureCandidate(candidate) && ok;
    }
    return ok ? 0 : 1;
}
