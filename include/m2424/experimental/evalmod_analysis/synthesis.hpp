#pragma once

#include "m2424/bootstrap_candidates.hpp"
#include "m2424/experimental/evalmod_analysis/approximation_lab.hpp"
#include "m2424/experimental/evalmod_analysis/cost_model.hpp"
#include "m2424/experimental/evalmod_analysis/exact_modular_oracle.hpp"

#include <optional>
#include <limits>

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
    std::size_t maxDataModulusBits{std::numeric_limits<std::size_t>::max()};
    double maxLatencyMs{std::numeric_limits<double>::infinity()};
    std::size_t maxWorkingSetBytes{std::numeric_limits<std::size_t>::max()};
};

enum class EvalModApproximationFamily {
    PeriodicSineBaseline,
    MultiIntervalLeastSquaresPrototype,
    MultiIntervalMinimax,
    MinimaxInverseSine
};

enum class EvalModCandidateStage {
    Generated, GridDiagnosed, IntervalCertified, CircuitCompiled, ScaleScheduled, BackendValidated
};

enum class EvalModRejectionReason {
    None, ApproximationError, Uncertified, ArithmeticErrorUnknown, InsufficientLevels,
    ModulusBudget, SecurityBudget, ScaleScheduleFailure
};

enum class EvalModOperation {
    Input, EncodeConstant, Add, MultiplyPlain, MultiplyCipher, Relinearize, Rescale
};

struct EvalModDagNode {
    EvalModOperation operation{};
    std::vector<std::size_t> inputs;
    std::size_t chainIndex{};
    ExactScale inputScale;
    ExactScale outputScale;
    std::string constantDecimal;
};

struct CompiledEvalModCircuit {
    std::vector<EvalModDagNode> nodes;
    EvalModCircuitCost cost;
    std::size_t outputNode{};
    std::size_t babyStep{};
    ExactScale normalizationGain;
    ExactScale denormalizationGain;
};

struct EvalModScaleStage {
    std::size_t inputScaleBits{};
    std::size_t multiplicationGrowthBits{};
    std::size_t rescalePrimeBits{};
    std::size_t outputScaleBits{};
    std::size_t availableModulusBits{};
    std::size_t requiredHeadroomBits{};
};

struct EvalModCandidate {
    EvalModApproximationFamily family{};
    EvalModDomain domain;
    EvalModPolynomial polynomial;
    CompiledEvalModCircuit compiledCircuit;
    std::vector<EvalModScaleStage> scaleSchedule;
    EvalModGridDiagnostic diagnostic;
    EvalModIntervalCertificate intervalCertificate;
    BootstrapPropagationBounds propagationBounds;
    std::optional<double> polynomialArithmeticError;
    double predictedBootstrapError{};
    EvalModCostEstimate cost;
    EvalModCandidateStage stage{EvalModCandidateStage::Generated};
    EvalModRejectionReason rejectionReason{EvalModRejectionReason::None};
    bool satisfiesNumericalTarget{};
    bool satisfiesLevelBudget{};
    bool intervalCertified{};
    bool executable{};
};

struct EvalModSynthesisResult {
    EvalModProblem problem;
    EvalModDomain domain;
    std::vector<EvalModCandidate> candidates;
    std::optional<std::size_t> provisionalSelection;
};

EvalModSynthesisResult synthesizeEvalMod(const EvalModProblem& problem);
CompiledEvalModCircuit compileEvalModPolynomial(const EvalModPolynomial& polynomial,
                                               const EvalModProblem& problem,
                                               std::size_t babyStep);
double evaluateEvalModPolynomial(const EvalModPolynomial& polynomial, double input);
std::string evalModSynthesisJson(const EvalModSynthesisResult& result);
std::string evalModSynthesisCsv(const EvalModSynthesisResult& result);

} // namespace m2424::experimental
