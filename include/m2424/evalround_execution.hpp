#pragma once
#include "m2424/evalround.hpp"
#include <functional>
#include <memory>
#include <optional>

namespace m2424 {
namespace experimental {
class EvalRoundExecutionCompiler;
struct EvalRoundBinaryDigitSearchResult;
}
// Deliberately no metadata scale rewrite operation.
enum class EvalRoundOperation { Input, Multiply, MultiplyPlain, Add, Subtract,
    AddPlain, Relinearize, Rescale, ModSwitch, Conjugate };
enum class EvalRoundEvaluationKey { None, Relinearization, Conjugation };
struct EvalRoundExactScale {
    std::string numerator, denominator; // exact positive rational, decimal integers
    std::uint64_t binary64Bits{};
};
struct EvalRoundExactBound {
    std::string numerator, denominator;
    BootstrapBoundKind kind{BootstrapBoundKind::Unknown};
    std::string provenance;
    // Present only when an outward binary64 projection is finite.
    std::optional<double> outwardBinary64;
};
struct EvalRoundExecutionNode {
    EvalRoundOperation operation{};
    std::vector<std::size_t> inputs;
    std::string stage;
    std::size_t chainIndex{};
    std::size_t ciphertextComponents{};
    std::vector<std::uint64_t> activePrimes;
    std::vector<EvalRoundExactScale> inputScales;
    EvalRoundExactScale outputScale;
    EvalRoundExactScale arithmeticScale; // exact scale before binary64 rounding
    // Constants are exact rationals; their encoding error is included below.
    std::string constantNumerator, constantDenominator;
    double constantScale{};
    EvalRoundExactScale constantEncodingScale;
    std::string encodedConstantInteger;
    std::string representedConstantNumerator, representedConstantDenominator;
    BootstrapBound idealMagnitude, propagatedSemanticError, localArithmeticError, semanticError;
    BootstrapBound constantEncodingError, scaleRepresentationError;
    EvalRoundExactBound exactIdealMagnitude, exactPropagatedSemanticError;
    EvalRoundExactBound exactLocalArithmeticError, exactSemanticError;
    EvalRoundExactBound exactConstantEncodingError, exactScaleRepresentationError;
    // Exact lower margin Q/2 - ceil(scale*(idealMagnitude+semanticError)).
    std::string centeredHeadroomNumerator; // denominator 2
    std::string centeredHeadroomProvenance;
    EvalRoundEvaluationKey requiredKey{EvalRoundEvaluationKey::None};
};
struct EvalRoundExecutionDigitDiagnostics {
    BootstrapBound approximationError;
    EvalRoundExactBound backendExtractionError;
    EvalRoundExactBound initialCleanerError;
    std::vector<EvalRoundExactBound> cleanerLocalErrors;
    std::vector<EvalRoundExactBound> cleanerErrorAfterRounds; // index zero is a_0
    EvalRoundExactBound reconstructionLocalError;
    std::size_t chebyshevNodeCount{};
};
struct EvalRoundBoundedCandidateDiagnostics {
    std::size_t digitIndex{};
    std::string family;
    std::size_t degree{};
    BootstrapBound approximationError;
    EvalRoundExactBound backendExtractionError;
    EvalRoundExactBound initialCleanerError;
    std::vector<EvalRoundExactBound> cleanerLocalErrors;
    std::vector<EvalRoundExactBound> cleanerErrorAfterRounds;
    std::vector<std::size_t> reachableOptionOutputNodes, reachableOptionLevels;
    std::vector<double> reachableOptionScales;
    std::vector<std::string> reachableOptionCenteredHeadroomNumerators;
    std::size_t maximumReachableCleaningRounds{};
    BootstrapCertificationStatus firstTrajectoryFailure{
        BootstrapCertificationStatus::Certified};
    std::string firstTrajectoryFailureProvenance;
    std::size_t firstTrajectoryFailureRound{};
    EvalRoundExactBound firstTrajectoryFailureInputError;
    EvalRoundExactBound firstTrajectoryFailureCleanerLocalError;
    EvalRoundExactBound firstTrajectoryFailureNextError;
    EvalRoundExactBound minimumReachableDigitError;
    std::size_t extractionOutputNode{}, outputLevel{};
    double outputScale{};
    std::string outputCenteredHeadroomNumerator;
    std::size_t reachableNodes{}, ciphertextMultiplications{}, relinearizations{};
    std::size_t rescales{}, modSwitches{}, plaintextMultiplications{};
    bool bestExecutableByTotalInitialError{};
};
struct EvalRoundExecutionDiagnostics {
    std::vector<EvalRoundExecutionDigitDiagnostics> digits;
    // Complete backend evaluation of the fixed bounded K=64 candidate space.
    std::vector<EvalRoundBoundedCandidateDiagnostics> boundedCandidates;
    std::vector<std::size_t> bestBoundedCandidateByDigit;
    EvalRoundExactBound minimumWeightedDigitError;
    std::string boundedSearchFailureDimension;
    std::size_t constructedNodes{};
    std::size_t ciphertextMultiplications{},relinearizations{},rescales{},modSwitches{};
    std::size_t plaintextMultiplications{},criticalMultiplicativeDepth{},criticalPathLevelConsumption{};
    std::optional<double> minimumRuntimeScale,maximumRuntimeScale;
    std::string minimumCenteredHeadroomNumerator;
};
class EvalRoundExecutionPlan {
public:
    EvalRoundExecutionPlan();
    const BootstrapContractResult& certification() const;
    const std::vector<EvalRoundExecutionNode>& nodes() const;
    const EvalRoundPlan& mathematicalPlan() const;
    double integerErrorUpper() const;
    std::size_t outputNode() const;
    // Available on fail-closed results as far as compilation reached. Nodes
    // themselves remain unpublished unless the complete plan is Certified.
    const EvalRoundExecutionDiagnostics& diagnostics() const;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
    friend class experimental::EvalRoundExecutionCompiler;
    friend BootstrapContractResult preflightEvalRound(const SealAdapter&, const Cipher&, const EvalRoundExecutionPlan&);
    friend Cipher executeEvalRound(SealAdapter&, const Cipher&, const EvalRoundExecutionPlan&,
        const std::function<void(std::size_t, const Cipher&)>&);
};
BootstrapContractResult preflightEvalRound(const SealAdapter&, const Cipher&, const EvalRoundExecutionPlan&);
// Observer only exposes ciphertexts. Decryption is not part of production computation.
Cipher executeEvalRound(SealAdapter&, const Cipher&, const EvalRoundExecutionPlan&,
    const std::function<void(std::size_t, const Cipher&)>& observer = {});
Cipher executeEvalRound(SealAdapter&, const Cipher&, const EvalRoundPlan&) = delete;
struct EvalRoundPairResult { Cipher first, second; double integerErrorUpper{}; };
EvalRoundPairResult executeEvalRoundPair(SealAdapter&, const Cipher&, const Cipher&, const EvalRoundExecutionPlan&);
} // namespace m2424
