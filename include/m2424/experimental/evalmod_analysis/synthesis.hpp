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

/// Empirical model retained exclusively for diagnostics and benchmark prediction.
using EvalModEmpiricalArithmeticModel = EvalModArithmeticErrorModel;

enum class BoundKind {
    Deterministic,
    Probabilistic,
    Unknown
};

/// A semantic upper bound used by the certified path. `Unknown` is never
/// interpreted as zero and makes the enclosing certificate non-rigorous.
struct CertifiedBound {
    double upperBound{std::numeric_limits<double>::infinity()};
    std::string upperBoundDecimal{"inf"};
    BoundKind kind{BoundKind::Unknown};
    double log2FailureProbability{0.0};
    bool rigorous{};
    std::string provenance;
};

enum class EvalModCertificationStatus {
    Certified,
    InvalidInput,
    MissingExactModulusContext,
    InvalidScaleSchedule,
    PreparedConstantMismatch,
    DiscontinuityMarginViolation,
    NoEvalModErrorBudgetRemaining,
    InsufficientLevels,
    SecurityBudgetExceeded,
    ScaleScheduleInfeasible,
    ModSwitchHeadroomViolation,
    RigorousRescaleBoundUnavailable,
    RigorousKeySwitchBoundUnavailable,
    UnknownOperationBound,
    ApproximationInsufficient,
    ArithmeticNoiseTooLarge,
    FailureProbabilityExceeded,
    MissingEvaluationKeys,
    ContextMismatch,
    InputLevelMismatch,
    InputScaleMismatch,
    CiphertextRepresentationMismatch
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
    MultiIntervalChebyshev,
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

/// Ideal semantic scale and the exact binary64 metadata value produced by SEAL arithmetic.
struct EvalModScaleValue {
    ExactScale ideal;
    std::uint64_t runtimeBinary64Bits{};
    ExactScale runtimeExact;
};

/// Точное prime-aware состояние scale и modulus для каждого DAG-узла.
struct EvalModExactScaleSchedule {
    bool available{};
    bool valid{};
    bool rigorous{};
    std::vector<EvalModScaleValue> scaleValues;
    /// Legacy alias for ideal scales. Certified code uses `scaleValues` explicitly.
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
    std::uint64_t specialPrime{};
    std::size_t polyModulusDegree{};
    std::size_t secretCoefficientAbsSupport{1};
    std::size_t evaluationKeyNoiseCoefficientAbsSupport{};
    bool deterministicEvaluationKeyNoiseSupport{};
    std::vector<PreparedEvalModConstant> constants;
    bool rigorous{};
};

enum class EvaluationKeySamplerKind {
    CenteredBinomial,
    ClippedRoundedGaussian,
    Unknown
};

struct EvaluationKeyNoiseCertificateMetadata {
    EvaluationKeySamplerKind sampler{EvaluationKeySamplerKind::Unknown};
    double standardDeviation{3.2};
    std::size_t noiseCoefficientAbsSupport{};
    std::size_t secretCoefficientAbsSupport{};
    std::size_t polyModulusDegree{};
    std::uint64_t specialPrime{};
    std::vector<std::uint64_t> dataPrimes;
    std::size_t relevantEvaluationKeyComponents{};
    bool rigorous{};
    std::string provenance;
};

struct EvalModNodeErrorState {
    double valueAbsBound{};
    double absoluteErrorBound{};
    double relativeScaleErrorBound{};
    double noiseBound{};
};

struct EvalModNodeCertificate {
    CertifiedBound valueAbs;
    CertifiedBound semanticError;
    CertifiedBound localAddedError;
    EvalModScaleValue scale;
    ExactInteger modulus;
    double headroomBits{-std::numeric_limits<double>::infinity()};
};

struct EvalModArithmeticCertificate {
    EvalModCertificationStatus status{EvalModCertificationStatus::InvalidInput};
    std::vector<EvalModNodeCertificate> nodes;
    CertifiedBound outputError;
    EvaluationKeyNoiseCertificateMetadata keyNoiseMetadata;
    double log2FailureProbability{};
    std::optional<std::size_t> failingNode;
    bool rigorous{};
    std::string detail;
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
    EvalModArithmeticCertificate arithmeticCertificate;
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

struct PreparedEvalModPlan {
    std::array<std::uint64_t, 4> contextFingerprint{};
    std::size_t inputChainIndex{};
    std::uint64_t inputScaleBinary64Bits{};
    EvalModDomain domain;
    CertifiedBound inputError;
    CertifiedBound discontinuityMargin;
    ExactScale normalizationGain;
    ExactScale denormalizationGain;
    EvalModPolynomial approximation;
    CompiledEvalModCircuit circuit;
    std::vector<std::uint64_t> dataPrimes;
    std::uint64_t specialPrime{};
    PreparedEvalModConstants constants;
    EvalModExactScaleSchedule scaleSchedule;
    EvalModArithmeticCertificate arithmeticCertificate;
    CertifiedBound approximationError;
    CertifiedBound arithmeticError;
    CertifiedBound normalizedEvalModError;
    CertifiedBound denormalizedEvalModError;
    CertifiedBound bootstrapContribution;
    double normalizedEvalModBudget{};
    std::size_t levelsConsumed{};
    int securityBits{};
    double log2FailureProbability{};
    EvalModCertificationStatus status{EvalModCertificationStatus::InvalidInput};
    bool rigorous{};
    std::string detail;
};

struct EvalModPreflightResult {
    EvalModCertificationStatus status{EvalModCertificationStatus::InvalidInput};
    bool compatible{};
    std::string detail;
};

enum class EvalModPlanSearchStatus {
    Certified,
    NoCertifiedPlanInSearchSpace,
    ProfileResourceInfeasible,
    RigorousBoundUnavailable
};

struct EvalModPlanSearchResult {
    EvalModPlanSearchStatus status{EvalModPlanSearchStatus::NoCertifiedPlanInSearchSpace};
    std::optional<PreparedEvalModPlan> plan;
    std::vector<EvalModApproximationFamily> familiesSearched;
    std::size_t minimumDegree{};
    std::size_t maximumDegree{};
    std::size_t maximumDepth{};
    double bestApproximationBound{std::numeric_limits<double>::infinity()};
    double bestArithmeticBound{std::numeric_limits<double>::infinity()};
    EvalModCertificationStatus firstFailingGate{EvalModCertificationStatus::InvalidInput};
    bool globalImpossibilityProved{};
    std::string proofScope;
    std::string recommendedRigorousNextPath;
    std::string detail;
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
    bool runtimeScaleBitsMatch{};
    std::optional<std::size_t> firstRuntimeScaleMismatchNode;
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
EvalModArithmeticCertificate certifyEvalModDagArithmetic(
    const CompiledEvalModCircuit& circuit,
    const EvalModExactScaleSchedule& schedule,
    const PreparedEvalModConstants& preparedConstants,
    double inputValueAbsUpperBound,
    double inputSemanticErrorUpperBound);
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
Cipher executeEvalModCircuitDiagnostic(
    SealAdapter& adapter, const CompiledEvalModCircuit& circuit,
    const Cipher& coeffToSlotOutput,
    EvalModExecutionTrace* trace = nullptr,
    const std::vector<double>* semanticInput = nullptr,
    std::vector<EvalModNodeDifferential>* differentialTrace = nullptr);
Cipher executePreparedEvalMod(
    SealAdapter& adapter, const CompiledEvalModCircuit& circuit,
    const PreparedEvalModConstants& preparedConstants,
    const Cipher& coeffToSlotOutput,
    EvalModExecutionTrace* trace = nullptr,
    const std::vector<double>* semanticInput = nullptr,
    std::vector<EvalModNodeDifferential>* differentialTrace = nullptr);
EvalModCoeffToSlotResult executeEvalModAfterCoeffToSlot(
    SealAdapter& adapter, const CompiledEvalModCircuit& circuit,
    const PreparedEvalModConstants& preparedConstants,
    CoeffToSlotResult coeffToSlotOutput);
PreparedEvalModPlan prepareEvalMod(
    SealAdapter& adapter,
    const EvalModCandidate& candidate,
    const EvalModProblem& problem,
    const Cipher& evalModInput);
EvalModPlanSearchResult prepareEvalModSearch(
    SealAdapter& adapter,
    const EvalModSynthesisResult& synthesis,
    const Cipher& evalModInput);
EvalModPreflightResult preflightEvalMod(
    const SealAdapter& adapter,
    const PreparedEvalModPlan& plan,
    const Cipher& evalModInput);
Cipher applyEvalMod(
    SealAdapter& adapter,
    const PreparedEvalModPlan& plan,
    const Cipher& evalModInput,
    EvalModExecutionTrace* trace = nullptr);
EvalModBackendValidation validateEvalModCandidateBackend(EvalModCandidate& candidate,
                                                        const EvalModProblem& problem,
                                                        const CkksProfile& profile,
                                                        const std::vector<double>& rawInput);
EvalModArithmeticErrorModel calibrateEvalModArithmeticModel(
    const std::vector<EvalModBackendValidation>& measurements, double safetyFactor = 4.0);
std::string evalModSynthesisJson(const EvalModSynthesisResult& result);
std::string evalModSynthesisCsv(const EvalModSynthesisResult& result);

} // namespace m2424::experimental
