#include "m2424/experimental/evalmod_analysis/cost_model.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace m2424::experimental {

EvalModCostEstimate estimateEvalModCost(const EvalModCircuitCost& circuit,
                                       const EvalModBackendCostModel& backend,
                                       std::size_t scaleBits) {
    if (circuit.degree == 0 || scaleBits == 0 || circuit.peakLiveCiphertexts == 0
        || !std::isfinite(backend.ciphertextMultiplyMs)
        || backend.ciphertextMultiplyMs < 0.0 || !std::isfinite(backend.relinearizeMs)
        || backend.relinearizeMs < 0.0 || !std::isfinite(backend.rescaleMs)
        || backend.rescaleMs < 0.0 || backend.bytesPerLiveCiphertext == 0
        || !std::isfinite(backend.ciphertextPlaintextMultiplyMs)
        || backend.ciphertextPlaintextMultiplyMs < 0.0 || !std::isfinite(backend.additionMs)
        || backend.additionMs < 0.0
        || circuit.relinearizations > circuit.ciphertextMultiplications) {
        throw std::invalid_argument("invalid EvalMod cost input");
    }
    const auto maximum = std::numeric_limits<std::size_t>::max();
    if (circuit.peakLiveCiphertexts > maximum / backend.bytesPerLiveCiphertext
        || circuit.levelConsumption == maximum
        || circuit.levelConsumption + 1 > maximum / scaleBits
        || circuit.preparedPlaintextBytes > maximum - circuit.backendScratchBytes
        || circuit.peakLiveCiphertexts * backend.bytesPerLiveCiphertext
            > maximum - circuit.preparedPlaintextBytes - circuit.backendScratchBytes) {
        throw std::overflow_error("EvalMod cost estimate overflow");
    }
    const double latency = circuit.ciphertextMultiplications * backend.ciphertextMultiplyMs
            + circuit.ciphertextPlaintextMultiplications * backend.ciphertextPlaintextMultiplyMs
            + circuit.additions * backend.additionMs
            + circuit.relinearizations * backend.relinearizeMs
            + circuit.rescales * backend.rescaleMs;
    if (!std::isfinite(latency)) throw std::overflow_error("EvalMod latency estimate overflow");
    return {
        latency,
        circuit.peakLiveCiphertexts * backend.bytesPerLiveCiphertext
            + circuit.preparedPlaintextBytes + circuit.backendScratchBytes,
        (circuit.levelConsumption + 1) * scaleBits
    };
}

} // namespace m2424::experimental
