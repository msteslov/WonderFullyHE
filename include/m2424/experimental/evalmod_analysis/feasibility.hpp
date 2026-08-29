#pragma once

#include "m2424/experimental/evalmod_analysis/synthesis.hpp"

#include <string>
#include <vector>

namespace m2424::experimental {

enum class EvalModFeasibilityStatus {
    CandidateFeasible,
    CertificateClassInfeasible,
    RigorousBoundUnavailable
};

struct EvalModErrorSourceBudget {
    std::string source;
    double currentNormalizedBound{};
    double bestPlausibleNormalizedBound{};
    double log2FailureProbability{};
    std::size_t occurrences{};
    double propagatedFinalContribution{};
    bool rigorous{};
    std::string note;
};

struct EvalModProfileFeasibility {
    std::size_t polyModulusDegree{};
    std::size_t tc128CoeffModulusBits{};
    std::size_t stableScaleBits{};
    std::size_t maximumDataPrimes{};
    std::size_t maximumLevels{};
    double mandatoryRescaleFloor{};
    double maximumNormalizedBudget{};
    double approximationLowerBound{};
    bool excludedByArithmetic{};
    bool excludedByApproximation{};
    std::string stopCondition;
};

struct EvalModFamilyFeasibility {
    std::string family;
    std::size_t diagnosticMinimumDegree{};
    std::size_t diagnosticMaximumDegree{};
    std::size_t resourceDegreeUpperBound{};
    std::size_t maximumDepth{};
    double bestDiagnosticApproximation{std::numeric_limits<double>::infinity()};
    double mandatoryArithmeticFloor{};
    std::string status;
};

/// A finite, explicitly scoped feasibility result for the stable-scale SEAL
/// certificate class used by the prepared EvalMod runtime.
struct EvalModFeasibilityReport {
    EvalModFeasibilityStatus status{EvalModFeasibilityStatus::RigorousBoundUnavailable};
    bool globalImpossibilityProved{};
    std::size_t polyModulusDegree{};
    std::size_t integerBound{};
    double rho{};
    std::size_t currentSecurityBits{};
    std::size_t tc128CoeffModulusBits{};
    std::size_t currentCoeffModulusBits{};
    std::size_t maximumSecurityValidDataPrimes{};
    std::size_t maximumSecurityValidLevels{};
    std::size_t maximumResourceFeasibleDegreeUpperBound{};
    double normalizedTarget{};
    double finalTarget{};
    double outputGain{};
    double mandatoryCurrentRescaleFloor{};
    double mandatoryRescaleImprovementFactor{};
    double optimisticFirstUseCbdKeySwitchBound{};
    double optimisticKeySwitchPlusModDown{};
    double dominantCurrentBound{};
    double requiredImprovementFactor{};
    double compositeBranchMargin{};
    double quadraticApproximationLowerBound{};
    double centralSensitivityLowerBound{};
    std::vector<EvalModErrorSourceBudget> sources;
    std::vector<EvalModProfileFeasibility> profiles;
    std::vector<EvalModFamilyFeasibility> families;
    std::string dominantSource;
    std::string proofScope;
    std::string probabilisticStopCondition;
    std::string directStopCondition;
    std::string compositeStopCondition;
    std::string remainingAlternatives;
};

EvalModFeasibilityReport analyzeEvalModFeasibility(
    const EvalModProblem& problem,
    const PreparedEvalModPlan& representativePlan,
    const EvalModSynthesisResult* diagnosticSynthesis = nullptr);

std::string evalModFeasibilityJson(const EvalModFeasibilityReport& report);

} // namespace m2424::experimental
