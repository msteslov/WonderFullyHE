#include "m2424/bootstrap_candidates.hpp"
#include "m2424/coeff_to_slot_bsgs_plan.hpp"
#include "m2424/coeff_to_slot_contract.hpp"
#include "m2424/coeff_to_slot_fft_plan.hpp"
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
constexpr std::size_t kRuntimeTrials = 7;
constexpr std::size_t kAccuracyTrials = 7;
constexpr double kPi = 3.141592653589793238462643383279502884;

struct ComparisonProfile {
    const char* candidateId;
    std::vector<std::size_t> transformSizes;
    bool fftSupported;
};

struct Selection {
    std::string candidateId;
    std::vector<std::size_t> transformSizes;
};

struct Measurement {
    double keySetupMs{};
    double prepareMs{};
    double minApplyMs{};
    double medianApplyMs{};
    double maxError{};
    double medianError{};
    std::size_t levels{};
    std::size_t rotationKeyCount{};
    std::size_t galoisKeyBytes{};
    std::size_t preparedPlaintextBytes{};
    bool keysWithinLimit{};
    bool passes{};
};

const std::vector<ComparisonProfile>& comparisonProfiles() {
    static const std::vector<ComparisonProfile> profiles{
        {"validation_4_s55", {4}, false},
        {"validation_16_s55", {16}, false},
        {"validation_64_s55", {64}, false},
        {"memory_1024_s55", {64}, false},
        {"throughput_1024_s55", {64}, false},
        {"memory_4096_s55", {64}, false},
        {"throughput_4096_s55", {64}, false},
        {"balanced_8192_s55", {4, 8, 16}, true},
        {"precision_8192_s59", {4, 8, 16}, true},
        {"memory_8192_s55", {4, 8, 16}, true},
        {"throughput_16384_s55", {4, 8, 16, 32, 64}, true},
        {"precision_16384_s59", {4, 8, 16, 32, 64}, true},
    };
    return profiles;
}

template <class Function>
double elapsedMs(Function&& function) {
    const auto start = Clock::now();
    function();
    const auto finish = Clock::now();
    return std::chrono::duration<double, std::milli>(finish - start).count();
}

std::size_t fftLevels(std::size_t transformSlots) {
    std::size_t levels = 1; // bit-reversal
    while (transformSlots > 1) {
        transformSlots >>= 1U;
        ++levels;
    }
    return levels;
}

m2424::CkksProfile makeComparisonProfile(const m2424::BootstrapCandidate& candidate,
                                          std::size_t transformSlots,
                                          bool fftSupported) {
    const int scaleBits = static_cast<int>(std::ceil(candidate.targetScaleLog2));
    std::vector<int> chain{60};
    chain.insert(chain.end(), fftSupported ? fftLevels(transformSlots) : 1, scaleBits);
    chain.push_back(60);
    // Диагонали CoeffToSlot занимают физическую CKKS-ёмкость, поэтому benchmark
    // не ограничивает encoder логическим числом активных слотов кандидата.
    return {candidate.polyModulusDegree, std::move(chain), std::exp2(candidate.targetScaleLog2), 0};
}

std::vector<std::size_t> bsgsBabySteps(std::size_t transformSlots) {
    std::vector<std::size_t> steps;
    for (std::size_t step = 1; step < transformSlots; step *= 2) {
        if (step * step >= transformSlots / 4 && step * step <= transformSlots * 4) {
            steps.push_back(step);
        }
    }
    if (steps.empty()) {
        steps.push_back(1);
    }
    return steps;
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
Measurement measurePlan(const m2424::BootstrapCandidate& candidate,
                        const m2424::CkksProfile& profile,
                        Plan& plan,
                        const m2424::ComplexVector& input) {
    auto adapter = m2424::SealAdapter::create(profile);
    Measurement measurement;
    measurement.rotationKeyCount = plan.rotationSteps().size();
    measurement.levels = plan.requiredLevels();
    measurement.keySetupMs = elapsedMs([&] { adapter.generateKeys(plan.rotationSteps(), false); });
    measurement.galoisKeyBytes = adapter.galoisKeysSize();
    measurement.keysWithinLimit = static_cast<double>(measurement.galoisKeyBytes) / (1024.0 * 1024.0)
        <= static_cast<double>(candidate.resources.maxEvaluationKeyMiB);

    const auto encodedInput = adapter.encodeComplex(input);
    const auto encrypted = adapter.encrypt(encodedInput);
    const auto preflight = m2424::preflightCoeffToSlot(
        adapter, encrypted, m2424::makeCoeffToSlotContract(candidate), plan.requirements());
    measurement.prepareMs = elapsedMs([&] { plan.prepare(adapter, encrypted); });
    measurement.preparedPlaintextBytes = plan.preparedPlaintextBytes(adapter);

    const auto expected = directDft(input);
    std::vector<double> accuracySamples;
    accuracySamples.reserve(kAccuracyTrials);
    m2424::CipherInfo resultInfo;
    for (std::size_t trial = 0; trial < kAccuracyTrials; ++trial) {
        const auto validationInput = adapter.encrypt(encodedInput);
        const auto validationResult = plan.apply(adapter, validationInput);
        const auto decoded = adapter.decodeComplex(adapter.decrypt(validationResult));
        const m2424::ComplexVector actual(
            decoded.begin(), decoded.begin() + static_cast<std::ptrdiff_t>(input.size()));
        accuracySamples.push_back(maxError(expected, actual));
        resultInfo = adapter.info(validationResult);
    }
    measurement.maxError = *std::max_element(accuracySamples.begin(), accuracySamples.end());
    measurement.medianError = median(std::move(accuracySamples));

    std::vector<double> runtimeSamples;
    runtimeSamples.reserve(kRuntimeTrials);
    for (std::size_t trial = 0; trial < kRuntimeTrials; ++trial) {
        runtimeSamples.push_back(elapsedMs([&] {
            const auto result = plan.apply(adapter, encrypted);
            (void)result;
        }));
    }
    measurement.minApplyMs = *std::min_element(runtimeSamples.begin(), runtimeSamples.end());
    measurement.medianApplyMs = median(std::move(runtimeSamples));
    const auto inputInfo = adapter.info(encrypted);
    measurement.passes = preflight.ready && measurement.keysWithinLimit
        && measurement.maxError <= candidate.errorBudget.coeffToSlot
        && resultInfo.chainIndex + measurement.levels == inputInfo.chainIndex
        && std::abs(std::log2(resultInfo.scale) - std::log2(inputInfo.scale)) <= 0.25;
    return measurement;
}

void printRow(const m2424::BootstrapCandidate& candidate,
              const char* strategy,
              std::size_t transformSlots,
              std::size_t babyStep,
              const Measurement& measurement) {
    std::printf("%s,%zu,%zu,%s,%zu,%zu,%zu,%zu,%.3f,%.3f,%.3f,%.3f,%.3e,%.3e,%zu,%zu,%zu,%zu,%s\n",
                candidate.id.c_str(), candidate.polyModulusDegree, candidate.activeSlots, strategy,
                transformSlots, babyStep, kRuntimeTrials, kAccuracyTrials, measurement.keySetupMs,
                measurement.prepareMs, measurement.minApplyMs, measurement.medianApplyMs,
                measurement.maxError, measurement.medianError, measurement.levels,
                measurement.rotationKeyCount, measurement.galoisKeyBytes, measurement.preparedPlaintextBytes,
                measurement.passes ? "PASS" : "FAIL");
}

Selection parseSelection(int argc, char** argv) {
    Selection selection;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument == "--profile" && index + 1 < argc) {
            selection.candidateId = argv[++index];
        } else if (argument == "--size" && index + 1 < argc) {
            selection.transformSizes = {std::stoull(argv[++index])};
        } else {
            throw std::invalid_argument("usage: bench_coeff_to_slot [--profile candidate_id] [--size power_of_two]");
        }
    }
    return selection;
}

bool isSelected(const Selection& selection, const ComparisonProfile& profile, std::size_t size) {
    return (selection.candidateId.empty() || selection.candidateId == profile.candidateId)
        && (selection.transformSizes.empty() || selection.transformSizes.front() == size);
}

} // namespace

int main(int argc, char** argv) {
    try {
        const auto selection = parseSelection(argc, argv);
        std::printf("candidate,poly_modulus_degree,active_slots,strategy,transform_slots,baby_step,runtime_trials,"
                    "accuracy_trials,key_setup_ms,prepare_ms,min_apply_ms,median_apply_ms,max_abs_error,"
                    "median_abs_error,levels,rotation_key_count,galois_key_bytes,prepared_plaintext_bytes,status\n");
        bool found = false;
        bool ok = true;
        for (const auto& profileSelection : comparisonProfiles()) {
            const auto& candidate = m2424::bootstrapCandidateById(profileSelection.candidateId);
            for (const std::size_t size : profileSelection.transformSizes) {
                if (!isSelected(selection, profileSelection, size)) {
                    continue;
                }
                found = true;
                const auto profile = makeComparisonProfile(candidate, size, profileSelection.fftSupported);
                const auto input = makeInput(size);
                const std::size_t physicalSlotCount = candidate.polyModulusDegree / 2;
                if (profileSelection.fftSupported) {
                    m2424::CoeffToSlotFftPlan fft(size, physicalSlotCount);
                    const auto fftMeasurement = measurePlan(candidate, profile, fft, input);
                    printRow(candidate, "fft", size, 0, fftMeasurement);
                    ok = fftMeasurement.passes && ok;
                }
                for (const std::size_t babyStep : bsgsBabySteps(size)) {
                    m2424::CoeffToSlotBsgsPlan bsgs(size, physicalSlotCount, babyStep);
                    const auto bsgsMeasurement = measurePlan(candidate, profile, bsgs, input);
                    printRow(candidate, "bsgs", size, babyStep, bsgsMeasurement);
                    ok = bsgsMeasurement.passes && ok;
                }
            }
        }
        if (!found) {
            std::fprintf(stderr, "no configured comparison case matches selection\n");
            return 2;
        }
        return ok ? 0 : 1;
    } catch (const std::exception& error) {
        std::fprintf(stderr, "%s\n", error.what());
        return 2;
    }
}
