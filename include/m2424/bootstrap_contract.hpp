#pragma once

#include "m2424/seal_adapter.hpp"
#include <array>
#include <map>
#include <limits>
#include <optional>
#include <string>
#include <vector>

namespace m2424 {

struct BootstrapTarget {
    double targetAbsoluteError{1e-10};
    int targetSecurityBits{128};
    double maxFailureProbabilityLog2{-128};
    // Zero means no configured resource limit.
    std::size_t maxEvaluationKeyBytes{};
    double maxLatencyMs{};
};

enum class BootstrapCertificationStatus {
    Certified, InvalidInput, MissingExactModulusContext, InputScaleMismatch,
    DomainViolation, SparseSecretCertificateUnavailable, KeySwitchBoundUnavailable,
    ExtractionNotCertified, CleaningBoundUnavailable, SlotToCoeffGainUnavailable,
    ScaleScheduleInfeasible, HeadroomViolation, InsufficientLevels,
    MissingEvaluationKeys, ErrorBudgetExceeded, FailureProbabilityExceeded,
    SecurityBudgetExceeded, RequiredBoundUnavailable, ResourceBudgetExceeded,
    LiftBoundUnavailable, UnsupportedEvalRoundDomain, TestOnlyAssumption
};

enum class BootstrapBoundKind { Unknown, Deterministic, Probabilistic };
struct BootstrapBound {
    double upperBound{std::numeric_limits<double>::infinity()};
    BootstrapBoundKind kind{BootstrapBoundKind::Unknown};
    std::string provenance;
    // References to primary assumptions, not probabilities of derived bounds.
    std::vector<std::string> failureEventIds;
};

/// Exact factored integer qSource = product(sourcePrimes), never a bit-count proxy.
/// Resolve BEFORE ModRaise. The binary64 bit pattern represents Delta0 exactly.
struct BootstrapInputContext {
    std::vector<std::uint64_t> sourcePrimes;
    std::vector<std::uint64_t> raisedPrimes;
    std::uint64_t scaleBinary64Bits{};
    std::array<std::uint64_t, 4> contextFingerprint{};
    std::uint64_t specialPrime{};
    std::size_t chainIndex{};
};

/// Certified means this particular validation succeeded; input/target validation
/// alone is not a full-plan certificate. No executable Bootstrapper is exposed.
struct BootstrapContractResult {
    BootstrapCertificationStatus status{BootstrapCertificationStatus::InvalidInput};
    std::string gate;
    std::string provenance;
};
// Evidence is an auditable external estimator result, not SEAL's tc128 flag.
struct RlweSecurityEvidence {
    std::string familyId, statement, estimator, version, artifact, assumptions;
    std::optional<double> lowerSecurityBits;
};
struct PublicRlweFamily {
    std::string id;
    std::size_t degree{};
    std::vector<std::uint64_t> modulus;
    std::string secretDistribution,errorDistribution,relations,purpose,provenance;
    std::size_t samples{},components{};
    std::string statement;
    std::optional<double> searchSpaceCeilingBits;
    std::optional<RlweSecurityEvidence> evidence;
    BootstrapContractResult result;
};
struct BootstrapSecurityReport {
    std::vector<PublicRlweFamily> families;
    std::optional<double> minimumSecurityBits;
    BootstrapContractResult result;
};
BootstrapSecurityReport certifyPublicRlweFamilies(std::vector<PublicRlweFamily>,
    const std::vector<RlweSecurityEvidence>&,int targetBits);

struct BootstrapInputResolution {
    BootstrapContractResult result;
    std::optional<BootstrapInputContext> context;
};
BootstrapInputResolution resolveBootstrapInput(const SealAdapter&, const Cipher&,
    std::optional<std::uint64_t> expectedScaleBinary64Bits = std::nullopt);

/// Fixed required gates prevent an empty/partial list from certifying a plan.
/// Evidence is supplied by future stage verifiers; PR-0 does not generate proofs.
enum class BootstrapGate : std::size_t {
    Input, SparseSecret, KeySwitch, Domain, Extraction, Arithmetic,
    CoeffToSlotHP, CoeffToSlotLP, Reconstruction, Combination, SlotToCoeff,
    ScaleSchedule, Headroom, EvaluationKeys, Security, Count
};
struct BootstrapGateEvidence {
    bool verified{};
    std::string provenance;
    std::vector<BootstrapBound> requiredBounds;
};
/// ScaleSchedule and EvaluationKeys need only verified evidence; Security also
/// needs minimumSecurityBits. Other gates require at least one known bound.
bool bootstrapGateRequiresNumericalBound(BootstrapGate);

struct BootstrapFailureEvent {
    std::string id;
    double log2FailureProbability{std::numeric_limits<double>::quiet_NaN()};
    std::string provenance;
};
struct BootstrapPlanMetadata {
    std::string id;
    std::optional<BootstrapInputContext> input;
    std::size_t levelsUsed{};
    // Missing measurement is allowed only when its target limit is unset.
    std::optional<std::size_t> evaluationKeyBytes;
    std::optional<double> latencyMs;
};
struct BootstrapCertificate {
    std::array<BootstrapGateEvidence, static_cast<std::size_t>(BootstrapGate::Count)> gates;
    std::vector<BootstrapFailureEvent> failureEvents;
    BootstrapBound outputError;
    std::optional<int> minimumSecurityBits;
};
struct BootstrapTraceNode {
    std::string stage;
    std::size_t node{};
    BootstrapInputContext state;
    BootstrapBound valueAbs, semanticError, localError;
    std::optional<double> observedError;
    std::vector<std::uint64_t> activePrimes;
    std::size_t chainIndex{};
};
struct BootstrapTrace {
    std::optional<BootstrapSecurityReport> publicRlwe;
    BootstrapPlanMetadata plan;
    std::vector<BootstrapTraceNode> nodes;
    BootstrapCertificate certificate;
    BootstrapContractResult result;
    std::map<std::string,BootstrapBound> bounds;
    std::map<std::string,std::string> details;
    std::vector<BootstrapContractResult> gateResults;
};

/// Union bound over unique primary event IDs (v9 section 6.4), capped at 1.
/// Identical duplicate events count once; conflicting IDs or invalid events fail.
/// Uses outward-rounded arithmetic and bounded exp/log series, not libm accuracy
/// assumptions. Requires round-to-nearest, without fast-math reassociation.
/// Empty events return log2(0) = -infinity.
std::optional<double> bootstrapFailureLog2UpperBound(const std::vector<BootstrapFailureEvent>&);
BootstrapContractResult validateBootstrapTarget(const BootstrapTarget&);
BootstrapContractResult validateBootstrapCertificate(const BootstrapTarget&,
    const BootstrapPlanMetadata&, const BootstrapCertificate&);

} // namespace m2424
