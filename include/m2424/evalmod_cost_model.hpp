#pragma once

#include <cstddef>

namespace m2424 {

struct EvalModCircuitCost {
    std::size_t degree{};
    std::size_t multiplicativeDepth{};
    std::size_t ciphertextMultiplications{};
    std::size_t relinearizations{};
    std::size_t rescales{};
};

struct EvalModBackendCostModel {
    double ciphertextMultiplyMs{};
    double relinearizeMs{};
    double rescaleMs{};
    std::size_t bytesPerLiveCiphertext{};
};

struct EvalModCostEstimate {
    double latencyMs{};
    std::size_t peakWorkingSetBytes{};
    std::size_t requiredModulusBits{};
};

EvalModCostEstimate estimateEvalModCost(const EvalModCircuitCost& circuit,
                                       const EvalModBackendCostModel& backend,
                                       std::size_t scaleBits,
                                       std::size_t liveCiphertexts);

} // namespace m2424
