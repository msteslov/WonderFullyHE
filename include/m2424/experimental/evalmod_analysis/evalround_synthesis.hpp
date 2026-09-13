#pragma once

#include "m2424/experimental/evalmod_analysis/evalround_execution.hpp"
#include "m2424/experimental/evalmod_analysis/synthesis.hpp"

#include <limits>
#include <optional>

namespace m2424::experimental {

enum class EvalRoundDigitGeneratorStatus {
    Generated,
    NotApplicableToTarget,
    GenerationFailed
};

struct EvalRoundDigitSearchRecord {
    std::size_t digitIndex{};
    EvalModApproximationFamily family{};
    std::size_t requestedDegree{};
    EvalRoundDigitGeneratorStatus generatorStatus{EvalRoundDigitGeneratorStatus::GenerationFailed};
    bool generatorConverged{};
    double gridMaximumError{std::numeric_limits<double>::infinity()};
    double rigorousIntervalError{std::numeric_limits<double>::infinity()};
    double maximumCoefficientMagnitude{};
    bool cleanerInputDomainSatisfied{};
    std::string provenance;
    std::optional<EvalRoundDigitPolynomial> certificate;
};

struct EvalRoundBinaryDigitSearchConfig {
    std::vector<std::size_t> degrees{64, 128, 192, 256};
    std::vector<EvalModApproximationFamily> families{
        EvalModApproximationFamily::MultiIntervalMinimax,
        EvalModApproximationFamily::MultiIntervalChebyshev};
    std::size_t gridPointsPerInterval{3};
    std::size_t intervalSubdivisions{16};
};

struct EvalRoundBinaryDigitSearchResult {
    EvalRoundProblem problem;
    EvalRoundBinaryDigitSearchConfig config;
    std::vector<std::vector<int>> targetTable; // Rows I=-K,...,+K.
    std::vector<EvalRoundDigitSearchRecord> records;
    std::vector<std::optional<std::size_t>> selectedRecordByDigit;
    bool allDigitsCertified{}; // Eight finite outward interval certificates.
    bool allCleanerInputDomainsSatisfied{};
    bool extractorCertified{};
    std::string status;
};

EvalRoundBinaryDigitSearchResult searchEvalRoundBinaryDigitPolynomials(
    const EvalRoundProblem&, const EvalRoundBinaryDigitSearchConfig& = {});

/// Builds the mathematical/reference candidate only from all selected outward
/// certificates. A missing digit produces an unverified candidate and cannot
/// be promoted by the planner or backend compiler.
EvalRoundCandidate makeEvalRoundBinaryPolynomialCandidate(
    const EvalRoundBinaryDigitSearchResult&);

} // namespace m2424::experimental
