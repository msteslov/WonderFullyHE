#pragma once

#include "m2424/bootstrap_candidates.hpp"
#include "m2424/experimental/evalmod_analysis/approximation_lab.hpp"
#include "m2424/experimental/evalmod_analysis/cost_model.hpp"
#include "m2424/experimental/evalmod_analysis/exact_modular_oracle.hpp"

#include <optional>

namespace m2424::experimental {

struct EvalModProblem {
    ExactInteger qSource;
    ExactScale coeffToSlotScale;
    ExactScale outputScale;
    EvalModCiphertextModel ciphertextModel;
    std::size_t availableLevels{};
    std::size_t targetPrecisionBits{};
    double targetAbsoluteError{};
    EvalModBackendCostModel backendCost;
    double slotToCoeffOperatorNorm{1.0};
    double slotToCoeffAdditive{};
    double finalAdditive{};
};

enum class EvalModApproximationFamily {
    PeriodicSineBaseline,
    MultiIntervalLeastSquaresPrototype,
    MultiIntervalMinimax,
    MinimaxInverseSine
};

struct EvalModScaleStage {
    std::size_t inputScaleBits{};
    std::size_t multiplicationGrowthBits{};
    std::size_t rescalePrimeBits{};
    std::size_t requiredHeadroomBits{};
};

struct EvalModCandidate {
    EvalModApproximationFamily family{};
    EvalModDomain domain;
    EvalModPolynomial polynomial;
    EvalModCircuitCost circuit;
    std::vector<EvalModScaleStage> scaleSchedule;
    EvalModGridDiagnostic diagnostic;
    BootstrapPropagationBounds propagationBounds;
    double predictedBootstrapError{};
    EvalModCostEstimate cost;
    bool feasible{};
};

struct EvalModSynthesisResult {
    EvalModProblem problem;
    EvalModDomain domain;
    std::vector<EvalModCandidate> candidates;
    std::optional<std::size_t> selectedCandidate;
};

EvalModSynthesisResult synthesizeEvalMod(const EvalModProblem& problem);
double evaluateEvalModPolynomial(const EvalModPolynomial& polynomial, double input);
std::string evalModSynthesisJson(const EvalModSynthesisResult& result);
std::string evalModSynthesisCsv(const EvalModSynthesisResult& result);

} // namespace m2424::experimental
