#include "m2424/coeff_to_slot.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <optional>
#include <cstdlib>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

template <class Function>
double elapsedMs(Function&& function) {
    const auto start = Clock::now();
    function();
    return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}

m2424::RaisedCipher makeRaised(m2424::SealAdapter& adapter,
                               const std::vector<double>& values) {
    const auto encrypted = adapter.encrypt(adapter.encode(values));
    auto lowered = encrypted;
    const std::size_t drops = adapter.info(encrypted).chainIndex;
    for (std::size_t level = 0; level < drops; ++level) {
        lowered = adapter.rescaleToNext(adapter.multiplyPlain(
            lowered, adapter.encodeScalarAtScaleFor(1.0, std::exp2(60.0), lowered)));
    }
    return adapter.modRaiseToTop(lowered);
}

double maxHalfError(const std::vector<double>& expected,
                    std::size_t offset,
                    const m2424::ComplexVector& actual) {
    double error = 0.0;
    for (std::size_t index = 0; index < actual.size(); ++index) {
        error = std::max(error, std::abs(actual[index] - expected[index + offset]));
    }
    return error;
}

} // namespace

int main(int argc, char** argv) {
    constexpr std::size_t degree = 16384;
    constexpr std::size_t slots = degree / 2;
    const m2424::CkksProfile profile{
        degree, {60, 60, 60, 60, 60, 60, 60}, std::exp2(59.5), slots
    };
    const m2424::CoeffToSlotContract contract{
        "bench_full_s59_5", slots, degree, 59.5, 59.5, 0.25, 2e-10
    };

    m2424::SealAdapter adapter;
    double contextMs = elapsedMs([&] { adapter = m2424::SealAdapter::create(profile); });

    const auto ranked = m2424::CoeffToSlotPlan::rankFactorizations(degree, 4, 5);
    const int candidateIndex = argc > 1 ? std::atoi(argv[1]) : -1;
    if (candidateIndex >= static_cast<int>(ranked.size())) return 2;
    std::optional<m2424::CoeffToSlot> transform;
    double planMs = elapsedMs([&] {
        if (candidateIndex >= 0) {
            transform.emplace(degree, ranked[static_cast<std::size_t>(candidateIndex)].factorization);
        } else {
            transform.emplace(degree);
        }
    });
    const auto requirements = transform->requirements();
    auto keySteps = requirements.rotationSteps;
    keySteps.push_back(0);
    double keysMs = elapsedMs([&] { adapter.generateKeys(keySteps, false); });

    std::vector<double> values(slots);
    for (std::size_t index = 0; index < slots; ++index) {
        values[index] = static_cast<double>(static_cast<int>((index * 29) % 41) - 20) / 64.0;
    }

    std::optional<m2424::RaisedCipher> firstRaised;
    double firstInputMs = elapsedMs([&] { firstRaised.emplace(makeRaised(adapter, values)); });
    std::vector<double> oracle;
    double oracleMs = elapsedMs([&] {
        oracle = adapter.decryptRaisedCoefficientsAtRaisedModulus(*firstRaised);
    });

    std::optional<m2424::PreparedCoeffToSlotPlan> prepared;
    double prepareMs = elapsedMs([&] {
        prepared.emplace(transform->prepare(adapter, *firstRaised, contract));
    });
    std::size_t preparedBytes = 0;
    double preparedSizeMs = elapsedMs([&] {
        preparedBytes = prepared->serializedPlaintextBytes(adapter);
    });

    std::optional<m2424::CoeffToSlotResult> firstResult;
    double firstApplyMs = elapsedMs([&] {
        firstResult.emplace(transform->apply(
            adapter, std::move(*firstRaised), contract, *prepared));
    });

    std::optional<m2424::RaisedCipher> secondRaised;
    double secondInputMs = elapsedMs([&] { secondRaised.emplace(makeRaised(adapter, values)); });
    std::vector<double> secondOracle;
    double secondOracleMs = elapsedMs([&] {
        secondOracle = adapter.decryptRaisedCoefficientsAtRaisedModulus(*secondRaised);
    });
    std::optional<m2424::CoeffToSlotResult> secondResult;
    double warmApplyMs = elapsedMs([&] {
        secondResult.emplace(transform->apply(
            adapter, std::move(*secondRaised), contract, *prepared));
    });

    m2424::ComplexVector firstSlots;
    m2424::ComplexVector secondSlots;
    double decodeMs = elapsedMs([&] {
        firstSlots = adapter.decodeComplex(adapter.decrypt(secondResult->slotCipherFirst));
        secondSlots = adapter.decodeComplex(adapter.decrypt(secondResult->slotCipherSecond));
    });
    const double firstError = maxHalfError(secondOracle, 0, firstSlots);
    const double secondError = maxHalfError(secondOracle, slots, secondSlots);
    const auto metrics = transform->plan().metrics();

    std::printf("phase,ms\n");
    std::printf("context,%.3f\nplan,%.3f\nkeys,%.3f\ninput_first,%.3f\n",
                contextMs, planMs, keysMs, firstInputMs);
    std::printf("oracle,%.3f\nprepare,%.3f\nprepared_size,%.3f\n",
                oracleMs, prepareMs, preparedSizeMs);
    std::printf("apply_first,%.3f\ninput_second,%.3f\noracle_second,%.3f\napply_warm,%.3f\ndecode,%.3f\n",
                firstApplyMs, secondInputMs, secondOracleMs, warmApplyMs, decodeMs);
    std::printf("metric,value\nplaintexts,%zu\nprepared_bytes,%zu\n",
                prepared->plaintextCount(), preparedBytes);
    std::printf("multiplications,%zu\nrotations,%zu\nadditions,%zu\nrescales,%zu\nkeys,%zu\n",
                metrics.plaintextMultiplicationsPerApply, metrics.rotationsPerApply,
                metrics.additionsPerApply, metrics.rescalesPerApply,
                metrics.uniqueEvaluationKeys);
    std::printf("hoisted_decompositions,%zu\nhoisted_automorphisms,%zu\nordinary_giant_rotations,%zu\n",
                metrics.hoistedDecompositionsPerApply,
                metrics.hoistedAutomorphismsPerApply,
                metrics.ordinaryGiantRotationsPerApply);
    std::printf("inner_moddowns,%zu\nfinal_moddowns,%zu\n",
                metrics.innerModDownsPerApply, metrics.finalModDownsPerApply);
    std::printf("first_error,%.3e\nsecond_error,%.3e\n", firstError, secondError);
    for (std::size_t index = 0; index < ranked.size(); ++index) {
        m2424::CoeffToSlotPlan candidatePlan(degree, ranked[index].factorization);
        const auto candidateMetrics = candidatePlan.metrics();
        std::printf("candidate_%zu,", index);
        for (std::size_t radix : ranked[index].factorization.radices) {
            std::printf("%zu-", radix);
        }
        std::printf("diag=%zu,exec_rot=%zu,keys=%zu,actual_mul=%zu,actual_rot=%zu,actual_keys=%zu\n",
                    ranked[index].estimatedDiagonals,
                    ranked[index].estimatedRotationsPerOutput,
                    ranked[index].estimatedUniqueRotations,
                    candidateMetrics.plaintextMultiplicationsPerApply,
                    candidateMetrics.rotationsPerApply,
                    candidateMetrics.uniqueEvaluationKeys);
    }
    return firstError <= contract.maxAbsError && secondError <= contract.maxAbsError ? 0 : 1;
}
