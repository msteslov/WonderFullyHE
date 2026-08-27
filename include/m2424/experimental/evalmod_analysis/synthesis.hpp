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

struct EvalModArithmeticErrorModel {
    double encodingAbsolute{};
    double additionAbsolute{};
    double multiplicationRelative{};
    double relinearizationAbsolute{};
    double rescaleAbsolute{};
    double modSwitchAbsolute{};
    double metadataScaleRelative{};
    bool calibrated{};
    bool rigorous{};
    std::string provenance;
};

struct EvalModPrecisionBudget {
    double implementation{};
    double approximation{};
};

/// Точные активные data-primes на входе EvalMod и special prime зафиксированного SEALContext.
/// Для входного ciphertext `dataPrimes` получают через `SealAdapter::coeffModulusValues`.
struct EvalModExactModulusContext {
    std::vector<std::uint64_t> dataPrimes;
    std::uint64_t specialPrime{};
};

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
    EvalModArithmeticErrorModel arithmeticErrorModel;
    EvalModPrecisionBudget precisionBudget;
    std::optional<EvalModExactModulusContext> exactModulusContext;
    bool requireCertifiedScaleSchedule{};
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
    ModulusBudget, SecurityBudget, ScaleScheduleFailure, HeadroomViolation
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
    double maxPlannedScaleDriftLog2{0.55};
};

struct EvalModScaleStage {
    std::size_t inputScaleBits{};
    std::size_t multiplicationGrowthBits{};
    std::size_t rescalePrimeBits{};
    std::size_t outputScaleBits{};
    std::size_t availableModulusBits{};
    std::size_t requiredHeadroomBits{};
};

/// Точное prime-aware состояние scale и modulus для каждого DAG-узла.
struct EvalModExactScaleSchedule {
    bool available{};
    bool valid{};
    bool rigorous{};
    std::vector<ExactScale> nodeScales;
    std::vector<ExactInteger> nodeModuli;
    std::vector<std::uint64_t> rescalePrimes;
    std::optional<std::size_t> failingNode;
    std::string failure;
};

struct EvalModModSwitchHeadroomGate {
    std::size_t node{};
    ExactInteger targetModulus;
    ExactInteger encodedMagnitudeUpperBound;
    double headroomBits{};
    bool proved{};
    bool rigorous{};
};

struct EvalModHeadroomCertificate {
    bool available{};
    bool valid{};
    bool rigorous{};
    std::vector<EvalModModSwitchHeadroomGate> modSwitchGates;
    std::optional<std::size_t> failingNode;
    std::string failure;
};

struct PreparedEvalModConstant {
    std::size_t node{};
    std::string decimal;
    ExactScale encodingScale;
    ExactInteger roundedScaledInteger;
    std::string encodingErrorUpperBoundDecimal;
    double encodingErrorUpperBound{};
    Plain plaintext;
    bool rigorous{};
};

struct PreparedEvalModConstants {
    std::vector<std::uint64_t> inputDataPrimes;
    std::vector<PreparedEvalModConstant> constants;
    bool rigorous{};
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
    double estimatedBootstrapError{};
    std::optional<double> certifiedBootstrapError;
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
    EvalModExactScaleSchedule exactScaleSchedule;
    EvalModHeadroomCertificate headroomCertificate;
};

struct EvalModExecutionTrace {
    std::size_t levelsConsumed{};
    double outputScale{};
    std::vector<CipherInfo> nodeStates;
};

struct EvalModNodeDifferential {
    std::size_t node{};
    EvalModOperation operation{};
    std::size_t chainIndex{};
    double actualScale{};
    double mpfrSemanticValue{};
    double decryptedCkksValue{};
    double absoluteError{};
    double errorIncrease{};
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
    std::vector<EvalModNodeDifferential> differentialTrace;
    std::optional<std::size_t> firstImplementationBudgetExceedingNode;
    bool preparedConstantsUsed{};
    double maxPreparedConstantEncodingError{};
};

EvalModSynthesisResult synthesizeEvalMod(const EvalModProblem& problem);
CompiledEvalModCircuit compileEvalModPolynomial(const EvalModPolynomial& polynomial,
                                               const EvalModProblem& problem,
                                               std::size_t babyStep);
EvalModExactScaleSchedule buildExactEvalModScaleSchedule(
    const CompiledEvalModCircuit& circuit,
    const EvalModExactModulusContext& context,
    bool rigorous);
EvalModHeadroomCertificate certifyEvalModModSwitchHeadroom(
    const CompiledEvalModCircuit& circuit,
    const EvalModExactScaleSchedule& schedule,
    const std::vector<EvalModNodeErrorState>& nodeErrors,
    bool errorBoundsRigorous);
PreparedEvalModConstants prepareEvalModConstants(
    SealAdapter& adapter,
    const CompiledEvalModCircuit& circuit,
    const Cipher& evalModInput);
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
                             EvalModExecutionTrace* trace = nullptr,
                             const std::vector<double>* semanticInput = nullptr,
                             std::vector<EvalModNodeDifferential>* differentialTrace = nullptr,
                             const PreparedEvalModConstants* preparedConstants = nullptr);
EvalModCoeffToSlotResult executeEvalModAfterCoeffToSlot(
    SealAdapter& adapter, const CompiledEvalModCircuit& circuit,
    CoeffToSlotResult coeffToSlotOutput);
EvalModBackendValidation validateEvalModCandidateBackend(EvalModCandidate& candidate,
                                                        const EvalModProblem& problem,
                                                        const CkksProfile& profile,
                                                        const std::vector<double>& rawInput);
EvalModArithmeticErrorModel calibrateEvalModArithmeticModel(
    const std::vector<EvalModBackendValidation>& measurements, double safetyFactor = 4.0);
std::string evalModSynthesisJson(const EvalModSynthesisResult& result);
std::string evalModSynthesisCsv(const EvalModSynthesisResult& result);

} // namespace m2424::experimental
