#include "m2424/experimental/evalmod_analysis/synthesis.hpp"

#include <cmath>
#include <cstdio>
#include <cstdint>

namespace em = m2424::experimental;

int main() {
    const em::ExactInteger qSource = em::ExactInteger(1) << 60;
    const em::ExactInteger ctsScale = em::ExactInteger(1) << 59;
    const em::ExactInteger outputScale = em::ExactInteger(1) << 57;
    em::EvalModProblem problem{
        qSource, em::ExactScale::rational(ctsScale, 1), em::ExactScale::rational(outputScale, 1),
        {em::TailModel::Deterministic, 32768, 1.0, 0.0, 0.08, 0.0, 0.0, 0.0, 0.0,
         -128.0, 384, {"non-unit backend deterministic validation", "CKKS arithmetic",
                       "bounded slots and exact binary scales"}},
        7, 59, 1e-10, {1.0, 0.2, 0.05, 0.2, 0.1, 4096}, 1.0, 0.0, 0.0
    };
    auto synthesis = em::synthesizeEvalMod(problem);
    auto selected = synthesis.candidates.end();
    for (auto it = synthesis.candidates.begin(); it != synthesis.candidates.end(); ++it) {
        if (it->family == em::EvalModApproximationFamily::PeriodicSineBaseline
            && it->approximationConverged && it->satisfiesLevelBudget) { selected = it; break; }
    }
    if (selected == synthesis.candidates.end()) return 1;
    auto& candidate = *selected;
    const double rho = synthesis.domain.normalizedResidualBound;
    const std::vector<double> zValues{-1.0 - rho, -1.0 + rho, -0.037, 0.0, 0.061,
                                      1.0 - rho, 1.0 + rho};
    std::vector<double> rawInput;
    for (double z : zValues) rawInput.push_back(2.0 * z); // normalization gain = 1/2.
    std::uint64_t randomState = 0x9e3779b97f4a7c15ULL;
    for (std::int64_t integer = -static_cast<std::int64_t>(synthesis.domain.integerBound);
         integer <= static_cast<std::int64_t>(synthesis.domain.integerBound); ++integer) {
        rawInput.push_back(2.0 * (integer - rho));
        rawInput.push_back(2.0 * (integer + rho));
        for (int sample = 0; sample < 3; ++sample) {
            randomState = randomState * 6364136223846793005ULL + 1;
            const double unit = static_cast<double>(randomState >> 11) * 0x1.0p-53;
            rawInput.push_back(2.0 * (integer - rho + 2.0 * rho * unit));
        }
    }
    const auto expected = em::evaluateEvalModReferenceMpfr(
        candidate.polynomial, candidate.compiledCircuit.normalizationGain,
        candidate.compiledCircuit.denormalizationGain, rawInput, 512);
    const auto smokeReference = em::executeEvalModDagPlaintext(candidate.compiledCircuit, rawInput);
    bool independentReferenceMatches = expected.size() == smokeReference.size();
    for (std::size_t index = 0; index < expected.size(); ++index)
        independentReferenceMatches = independentReferenceMatches
            && std::abs(expected[index] - smokeReference[index]) < 1e-10;

    const m2424::CkksProfile profile{32768, {60, 59, 59, 59, 59, 59, 59, 59, 60},
                                      std::exp2(59.0), rawInput.size()};
    const auto validation = em::validateEvalModCandidateBackend(candidate, problem, profile, rawInput);
    const m2424::CkksProfile shortProfile{32768, {60, 59, 60}, std::exp2(59.0), rawInput.size()};
    auto shortCandidate = candidate;
    const auto shortValidation = em::validateEvalModCandidateBackend(
        shortCandidate, problem, shortProfile, rawInput);
    auto prototype = synthesis.candidates[1];
    const auto prototypeValidation = em::validateEvalModCandidateBackend(
        prototype, problem, profile, rawInput);
    bool minimaxMatrixPassed = true;
    std::vector<em::EvalModBackendValidation> calibrationMeasurements{validation};
    const std::vector<std::pair<std::size_t, std::size_t>> minimaxCases{{7, 2}, {9, 3}, {11, 4}};
    for (const auto [degree, babyStep] : minimaxCases) {
        auto found = synthesis.candidates.end();
        for (auto it = synthesis.candidates.begin(); it != synthesis.candidates.end(); ++it) {
            if (it->family == em::EvalModApproximationFamily::MultiIntervalMinimax
                && it->compiledCircuit.cost.degree == degree
                && it->compiledCircuit.babyStep == babyStep) { found = it; break; }
        }
        minimaxMatrixPassed = minimaxMatrixPassed && found != synthesis.candidates.end();
        if (found == synthesis.candidates.end()) continue;
        for (int seedTrial = 0; seedTrial < 2; ++seedTrial) {
            auto minimax = *found;
            const auto measured = em::validateEvalModCandidateBackend(
                minimax, problem, profile, rawInput);
            calibrationMeasurements.push_back(measured);
            minimaxMatrixPassed = minimaxMatrixPassed && measured.executionSucceeded
                && minimax.backendRunnable && minimax.backendMeasured
                && measured.executedNodes == minimax.compiledCircuit.nodes.size();
        }
    }
    const auto calibratedModel = em::calibrateEvalModArithmeticModel(calibrationMeasurements, 4.0);

    const bool ok = validation.passed && independentReferenceMatches
        && candidate.stage == em::EvalModCandidateStage::BackendMeasured
        && candidate.measuredBackendError.has_value() && candidate.executable
        && candidate.circuitValid && candidate.backendRunnable && candidate.backendMeasured
        && !candidate.rigorouslyValidated
        && candidate.intervalCertified && candidate.approximationCertified
        && !candidate.arithmeticErrorRigorous && !candidate.arithmeticErrorCertified
        && validation.executedNodes == candidate.compiledCircuit.nodes.size()
        && validation.executionSucceeded && !validation.matchesPolynomialReference
        && !validation.matchesEvalModTarget
        && validation.implementationError == validation.maxAbsoluteError
        && std::abs(std::log2(validation.outputScale) - 59.0) < 0.2
        && !prototypeValidation.passed && !prototype.executable
        && minimaxMatrixPassed
        && calibratedModel.calibrated && !calibratedModel.rigorous
        && calibratedModel.encodingAbsolute >= validation.implementationError
        && !shortValidation.passed && shortValidation.executedNodes == 0
        && shortValidation.failure == "profile_chain_too_short";
    if (!validation.passed) {
        for (std::size_t index = 0; index < candidate.compiledCircuit.nodes.size(); ++index)
            std::printf("node[%zu]=%d chain=%zu\n", index,
                        static_cast<int>(candidate.compiledCircuit.nodes[index].operation),
                        candidate.compiledCircuit.nodes[index].chainIndex);
    }
    std::printf("[test_evalmod_backend] error=%.3e nodes=%zu failure=%s %s\n",
                validation.maxAbsoluteError, validation.executedNodes, validation.failure.c_str(),
                ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
