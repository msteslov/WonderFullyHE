#pragma once

#include <cstddef>

namespace m2424::experimental {

struct EvalModCircuitCost {
    std::size_t degree{};
    std::size_t multiplicativeDepth{};
    std::size_t ciphertextMultiplications{};
    std::size_t ciphertextPlaintextMultiplications{};
    std::size_t additions{};
    std::size_t relinearizations{};
    std::size_t rescales{};
    std::size_t levelConsumption{};
    std::size_t peakLiveCiphertexts{};
    std::size_t preparedPlaintextBytes{};
    std::size_t backendScratchBytes{};
    std::size_t modulusSwitches{};
    std::size_t scaleAlignments{};
    std::size_t plaintextAdditions{};
};

struct EvalModBackendCostModel {
    double ciphertextMultiplyMs{};
    double ciphertextPlaintextMultiplyMs{};
    double additionMs{};
    double relinearizeMs{};
    double rescaleMs{};
    std::size_t bytesPerLiveCiphertext{};
    double modulusSwitchMs{};
    double scaleAlignmentMs{};
    double plaintextAdditionMs{};
};

struct EvalModCostEstimate {
    double latencyMs{};
    std::size_t peakWorkingSetBytes{};
    std::size_t dataModulusBitsLowerBound{};
};

EvalModCostEstimate estimateEvalModCost(const EvalModCircuitCost& circuit,
                                       const EvalModBackendCostModel& backend,
                                       std::size_t scaleBits);

} // namespace m2424::experimental
