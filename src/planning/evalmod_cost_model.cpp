#include "m2424/evalmod_cost_model.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace m2424 {

EvalModCostEstimate estimateEvalModCost(const EvalModCircuitCost& circuit,
                                       const EvalModBackendCostModel& backend,
                                       std::size_t scaleBits,
                                       std::size_t liveCiphertexts) {
    if (scaleBits == 0 || liveCiphertexts == 0 || !std::isfinite(backend.ciphertextMultiplyMs)
        || backend.ciphertextMultiplyMs < 0.0 || !std::isfinite(backend.relinearizeMs)
        || backend.relinearizeMs < 0.0 || !std::isfinite(backend.rescaleMs)
        || backend.rescaleMs < 0.0 || backend.bytesPerLiveCiphertext == 0
        || circuit.relinearizations > circuit.ciphertextMultiplications
        || circuit.rescales > circuit.ciphertextMultiplications) {
        throw std::invalid_argument("invalid EvalMod cost input");
    }
    if (liveCiphertexts > std::numeric_limits<std::size_t>::max() / backend.bytesPerLiveCiphertext
        || circuit.multiplicativeDepth > std::numeric_limits<std::size_t>::max() / scaleBits) {
        throw std::overflow_error("EvalMod cost estimate overflow");
    }
    return {
        circuit.ciphertextMultiplications * backend.ciphertextMultiplyMs
            + circuit.relinearizations * backend.relinearizeMs
            + circuit.rescales * backend.rescaleMs,
        liveCiphertexts * backend.bytesPerLiveCiphertext,
        (circuit.multiplicativeDepth + 1) * scaleBits
    };
}

} // namespace m2424
