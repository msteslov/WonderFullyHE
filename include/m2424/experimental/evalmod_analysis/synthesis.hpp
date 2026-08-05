#pragma once

#include "m2424/bootstrap_candidates.hpp"
#include "m2424/coeff_to_slot.hpp"
#include "m2424/seal_adapter.hpp"
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
    Generated, GridDiagnosed, IntervalCertified, CircuitCompiled, ScaleScheduled,
    BackendMeasured, BackendValidated
};

enum class EvalModRejectionReason {
    None, ApproximationError, ApproximationNotConverged, Uncertified,
    ArithmeticErrorUnknown, InsufficientLevels,
    ModulusBudget, SecurityBudget, ScaleScheduleFailure
};

enum class EvalModOperation {
    Input, EncodeConstant, Add, AddPlain, MultiplyPlain, MultiplyCipher,
    Relinearize, Rescale, ModSwitch, AlignScale
};

struct EvalModDagNode {
    EvalModOperation operation{};
    std::vector<std::size_t> inputs;
    std::size_t chainIndex{};
    ExactScale inputScale;
    ExactScale outputScale;
    std::string constantDecimal;
    double plannedScaleCorrectionLog2{};
    double expectedRelativeScaleError{};
};

struct CompiledEvalModCircuit {
    std::vector<EvalModDagNode> nodes;
    EvalModCircuitCost cost;
    std::size_t outputNode{};
    std::size_t babyStep{};
    ExactScale normalizationGain;
    ExactScale denormalizationGain;
    double maxMetadataScaleCorrectionLog2{1e-6};
};

struct EvalModScaleStage {
    std::size_t inputScaleBits{};
    std::size_t multiplicationGrowthBits{};
    std::size_t rescalePrimeBits{};
    std::size_t outputScaleBits{};
    std::size_t availableModulusBits{};
    std::size_t requiredHeadroomBits{};
};

struct EvalModNodeErrorState {
    double valueAbsBound{};
    double absoluteErrorBound{};
    double relativeScaleErrorBound{};
    double noiseBound{};
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
    bool arithmeticErrorRigorous{};
    bool approximationConverged{};
    double predictedBootstrapError{};
    EvalModCostEstimate cost;
    EvalModCandidateStage stage{EvalModCandidateStage::Generated};
    EvalModRejectionReason rejectionReason{EvalModRejectionReason::None};
    bool satisfiesNumericalTarget{};
    bool satisfiesLevelBudget{};
    bool intervalCertified{};
    bool executable{};
    std::optional<double> measuredBackendError;
    bool circuitValid{};
    bool backendRunnable{};
    bool backendMeasured{};
    bool approximationCertified{};
    bool arithmeticErrorCertified{};
    bool rigorouslyValidated{};
    std::optional<double> analyticalArithmeticBound;
    bool analyticalArithmeticBoundRigorous{};
    std::vector<EvalModNodeErrorState> nodeErrorStates;
};

struct EvalModExecutionTrace {
    std::size_t levelsConsumed{};
    double outputScale{};
    std::vector<CipherInfo> nodeStates;
};

struct EvalModCoeffToSlotResult {
    Cipher slotCipherFirst;
    Cipher slotCipherSecond;
    EvalModExecutionTrace firstTrace;
    EvalModExecutionTrace secondTrace;
};

struct EvalModSynthesisResult {
    EvalModProblem problem;
    EvalModDomain domain;
    std::vector<EvalModCandidate> candidates;
    std::optional<std::size_t> provisionalSelection;
};

struct EvalModBackendValidation {
    bool passed{};
    bool executionSucceeded{};
    bool matchesPolynomialReference{};
    bool matchesEvalModTarget{};
    bool predictionCoveredMeasurement{};
    double implementationError{};
    double approximationError{};
    double totalMeasuredError{};
    double maxAbsoluteError{};
    std::size_t executedNodes{};
    std::size_t outputChainIndex{};
    double outputScale{};
    std::string failure;
};

EvalModSynthesisResult synthesizeEvalMod(const EvalModProblem& problem);
CompiledEvalModCircuit compileEvalModPolynomial(const EvalModPolynomial& polynomial,
                                               const EvalModProblem& problem,
                                               std::size_t babyStep);
double evaluateEvalModPolynomial(const EvalModPolynomial& polynomial, double input);
std::vector<double> executeEvalModDagPlaintext(const CompiledEvalModCircuit& circuit,
                                               const std::vector<double>& input);
std::vector<double> evaluateEvalModReferenceMpfr(const EvalModPolynomial& polynomial,
                                                 const ExactScale& normalizationGain,
                                                 const ExactScale& denormalizationGain,
                                                 const std::vector<double>& rawInput,
                                                 std::size_t precisionBits = 384);
std::vector<double> evaluateExactEvalModTargetMpfr(
    const ExactScale& normalizationGain, const ExactScale& denormalizationGain,
    const std::vector<double>& rawInput, std::size_t precisionBits = 384);
bool isCompiledEvalModCircuitValid(const CompiledEvalModCircuit& circuit);
Cipher executeEvalModCircuit(SealAdapter& adapter, const CompiledEvalModCircuit& circuit,
                             const Cipher& coeffToSlotOutput,
                             EvalModExecutionTrace* trace = nullptr);
EvalModCoeffToSlotResult executeEvalModAfterCoeffToSlot(
    SealAdapter& adapter, const CompiledEvalModCircuit& circuit,
    CoeffToSlotResult coeffToSlotOutput);
EvalModBackendValidation validateEvalModCandidateBackend(EvalModCandidate& candidate,
                                                        const EvalModProblem& problem,
                                                        const CkksProfile& profile,
                                                        const std::vector<double>& rawInput);
std::string evalModSynthesisJson(const EvalModSynthesisResult& result);
std::string evalModSynthesisCsv(const EvalModSynthesisResult& result);

} // namespace m2424::experimental
