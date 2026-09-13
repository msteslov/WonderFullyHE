#pragma once

#include "m2424/bootstrap_contract.hpp"
#include <cstdint>
#include "m2424/decimal_polynomial.hpp"

namespace m2424 {

enum class EvalRoundRadix : unsigned { Binary = 2, BalancedTernary = 3 };
enum class EvalRoundExtractionMethod {
    PiecewiseReference, BinaryQuadraticK1, TernaryPhaseReferenceK1,
    ExternalPolynomial, DigitExtract
};
enum class EvalRoundRejectionReason {
    None, InvalidProblem, DomainViolation, InvalidCandidate, DigitExtractDomainViolation,
    ExtractionNotCertified, BoundUnavailable, CleaningDomainViolation,
    ErrorBudgetExceeded, FailureProbabilityExceeded
};
enum class EvalRoundPlanStatus { Rejected, Diagnostic, Certified };

struct EvalRoundProblem {
    std::uint32_t K{};
    double rho{};
    double requiredIntegerError{}; // Required caller budget; no precision default.
    std::size_t maxCleaningRoundsPerDigit{16}; // Bounded search, at most 64.
    double maxFailureProbabilityLog2{-128};
};
struct EvalRoundCost {
    double extraction{};
    double cleaningPerDigitIteration{1};
    double reconstruction{};
};
enum class EvalRoundPolynomialProof { Unknown, GridDiagnostic, BinaryQuadraticIdentity, OutwardInterval };
enum class EvalRoundDigitTarget { BinaryOffsetDigit, TernaryRoot };
struct EvalRoundDigitPolynomial {
    EvalRoundRadix radix{EvalRoundRadix::Binary};
    std::size_t digitIndex{};
    experimental::EvalModPolynomial polynomial;
    EvalRoundDigitTarget target{EvalRoundDigitTarget::BinaryOffsetDigit};
    std::uint32_t certifiedK{};
    double certifiedRho{};
    BootstrapBound approximationError;
    EvalRoundPolynomialProof proof{EvalRoundPolynomialProof::Unknown};
    bool verified{};
    std::string provenance;
    std::vector<std::string> failureEventIds;
    std::size_t intervalSubdivisions{64};
};
struct EvalRoundExtractionDescription {
    EvalRoundExtractionMethod method{EvalRoundExtractionMethod::ExternalPolynomial};
    std::string description;
    std::vector<std::size_t> polynomialDegrees;
    std::optional<double> digitExtractEpsilon;
    // Evidence must cover ALL intervals of this domain, not a sampled grid.
    bool verified{};
    std::uint32_t certifiedK{};
    double certifiedRho{};
    std::string provenance;
    std::vector<EvalRoundDigitPolynomial> polynomials;
};
struct EvalRoundDigitBounds {
    BootstrapBound extractionError;
    std::vector<BootstrapBound> cleaningLocalErrors; // One entry per usable round.
    BootstrapBound reconstructionLocalError;
};
struct EvalRoundCandidate {
    std::string id;
    EvalRoundRadix radix{EvalRoundRadix::Binary};
    EvalRoundExtractionDescription extraction;
    std::vector<EvalRoundDigitBounds> digits;
    std::vector<BootstrapFailureEvent> failureEvents;
    EvalRoundCost cost;
};
struct EvalRoundDigitCertificate {
    std::size_t cleaningIterations{};
    std::vector<double> errorAfterRounds; // Includes extraction at index zero.
    double weightedReconstructionError{};
    BootstrapBound extractionError;
    std::vector<BootstrapBound> cleaningLocalErrors; // Only the selected reachable rounds.
    BootstrapBound reconstructionLocalError;
};
struct EvalRoundPlan {
    EvalRoundProblem problem;
    EvalRoundCost cost;
    std::vector<BootstrapFailureEvent> failureEvents;
    std::string candidateId;
    EvalRoundRadix radix{EvalRoundRadix::Binary};
    EvalRoundExtractionDescription extraction;
    EvalRoundPlanStatus status{EvalRoundPlanStatus::Rejected};
    EvalRoundRejectionReason rejection{EvalRoundRejectionReason::InvalidProblem};
    std::string provenance;
    std::vector<EvalRoundDigitCertificate> digits;
    std::size_t totalCleaningIterations{};
    double integerErrorUpper{std::numeric_limits<double>::infinity()};
    double log2FailureProbabilityUpper{std::numeric_limits<double>::infinity()};
    double configuredCost{std::numeric_limits<double>::infinity()};
    // PR-2 certifies only mathematical/reference bounds under the input contract.
    // It supplies no ciphertext DAG, scale schedule, headroom or key certificate.
};
struct EvalRoundPlanningResult {
    std::vector<EvalRoundPlan> candidates;
    std::optional<std::size_t> selected;
    EvalRoundRejectionReason rejection{EvalRoundRejectionReason::InvalidProblem};
    std::string provenance;
};

std::size_t evalRoundDigitCount(std::uint32_t K, EvalRoundRadix);
std::vector<int> evalRoundIntegerDigits(std::int64_t I, std::uint32_t K, EvalRoundRadix);
std::int64_t reconstructEvalRoundInteger(const std::vector<int>&, std::uint32_t K, EvalRoundRadix);
/// Outward arithmetic for the hard v9 bounds. Throws on invalid radix/domain.
double evalRoundCleaningErrorUpper(EvalRoundRadix, double a, double localError);
double evalRoundReconstructionErrorUpper(EvalRoundRadix,
    const std::vector<double>& errors, const std::vector<double>& localErrors);

/// Exact-arithmetic/reference candidates with analytic whole-domain extraction
/// bounds. Piecewise targets support any K; quadratic/phase methods require K=1.
/// Their zero LOCAL errors refer to exact mathematical operations, never CKKS.
EvalRoundCandidate makeEvalRoundReferenceCandidate(const EvalRoundProblem&, EvalRoundRadix,
    EvalRoundExtractionMethod, const EvalRoundCost& = {});
EvalRoundPlan planEvalRoundCandidate(const EvalRoundProblem&, const EvalRoundCandidate&);
EvalRoundPlanningResult planEvalRound(const EvalRoundProblem&, const std::vector<EvalRoundCandidate>&);

} // namespace m2424
