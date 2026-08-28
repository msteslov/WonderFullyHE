#include "m2424/experimental/evalmod_analysis/synthesis.hpp"

#include <algorithm>
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
    problem.precisionBudget = {1e-10, 1e-10};
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
    m2424::SealAdapter contextAdapter = m2424::SealAdapter::create(profile);
    const std::vector<std::uint64_t> expectedDataPrimes{
        1152921504598720513ULL,
        576460752291954689ULL,
        576460752293134337ULL,
        576460752298180609ULL,
        576460752298835969ULL,
        576460752300015617ULL,
        576460752301391873ULL,
        576460752301785089ULL
    };
    const std::uint64_t expectedSpecialPrime = 1152921504606584833ULL;
    const em::EvalModExactModulusContext exactContext{
        contextAdapter.dataModulusValues(), contextAdapter.specialKeyModulusValue()};
    auto certifiedProblem = problem;
    certifiedProblem.exactModulusContext = exactContext;
    certifiedProblem.requireCertifiedScaleSchedule = true;
    const auto certifiedCircuit = em::compileEvalModPolynomial(
        candidate.polynomial, certifiedProblem, candidate.compiledCircuit.babyStep);
    const auto certifiedScaleSchedule = em::buildExactEvalModScaleSchedule(
        certifiedCircuit, exactContext, true);
    const auto certifiedAlignments = std::count_if(
        certifiedCircuit.nodes.begin(), certifiedCircuit.nodes.end(), [](const auto& node) {
            return node.operation == em::EvalModOperation::AlignScale;
        });

    em::CompiledEvalModCircuit oneMultiply;
    const auto scale59 = em::ExactScale::rational(em::ExactInteger(1) << 59, 1);
    const auto scale118 = em::ExactScale::rational(em::ExactInteger(1) << 118, 1);
    oneMultiply.nodes = {
        {em::EvalModOperation::Input, {}, 0, scale59, scale59},
        {em::EvalModOperation::MultiplyCipher, {0, 0}, 0, scale59, scale118},
        {em::EvalModOperation::Relinearize, {1}, 0, scale118, scale118},
        {em::EvalModOperation::Rescale, {2}, 1, scale118, scale59}
    };
    oneMultiply.outputNode = 3;
    const auto exactSchedule = em::buildExactEvalModScaleSchedule(
        oneMultiply, exactContext, true);
    auto metadataAligned = oneMultiply;
    metadataAligned.nodes.push_back({em::EvalModOperation::AlignScale, {3, 3}, 1,
                                     scale59, scale59});
    metadataAligned.outputNode = 4;
    const auto prohibitedSchedule = em::buildExactEvalModScaleSchedule(
        metadataAligned, exactContext, true);
    em::CompiledEvalModCircuit oneModSwitch;
    oneModSwitch.nodes = {
        {em::EvalModOperation::Input, {}, 0, scale59, scale59},
        {em::EvalModOperation::Input, {}, 1, scale59, scale59},
        {em::EvalModOperation::ModSwitch, {0, 1}, 1, scale59, scale59}
    };
    oneModSwitch.outputNode = 2;
    const auto modSwitchSchedule = em::buildExactEvalModScaleSchedule(
        oneModSwitch, exactContext, true);
    std::vector<em::EvalModNodeErrorState> boundedStates(3);
    for (auto& state : boundedStates) state = {1.0, 1e-12, 0.0, 1e-12};
    const auto headroom = em::certifyEvalModModSwitchHeadroom(
        oneModSwitch, modSwitchSchedule, boundedStates, true);
    auto overflowingStates = boundedStates;
    overflowingStates[2].valueAbsBound = 1e150;
    const auto noHeadroom = em::certifyEvalModModSwitchHeadroom(
        oneModSwitch, modSwitchSchedule, overflowingStates, true);
    const std::vector<double> precisionInput{-0.12, -0.06, 0.0, 0.06, 0.12};
    const auto validation = em::validateEvalModCandidateBackend(
        candidate, problem, profile, precisionInput);
    auto certifiedCandidate = candidate;
    certifiedCandidate.compiledCircuit = certifiedCircuit;
    certifiedCandidate.satisfiesLevelBudget = certifiedCircuit.cost.levelConsumption
        <= certifiedProblem.availableLevels;
    const auto certifiedValidation = em::validateEvalModCandidateBackend(
        certifiedCandidate, certifiedProblem, profile, precisionInput);
    auto planAdapter = m2424::SealAdapter::create(profile);
    planAdapter.generateKeys(std::vector<int>{}, true);
    const auto planInputCipher = planAdapter.encrypt(planAdapter.encode(precisionInput));
    const auto nonlinearPlan = em::prepareEvalMod(
        planAdapter, certifiedCandidate, certifiedProblem, planInputCipher);
    const auto planSearch = em::prepareEvalModSearch(
        planAdapter, synthesis, planInputCipher);
    bool operationCertificatesCover = true;
    bool sawCertifiedRescale = false;
    bool sawCertifiedRelinearize = false;
    for (const auto& differential : certifiedValidation.differentialTrace) {
        if (differential.operation != em::EvalModOperation::Rescale
            && differential.operation != em::EvalModOperation::Relinearize) continue;
        const auto& nodeCertificate =
            certifiedCandidate.arithmeticCertificate.nodes[differential.node];
        operationCertificatesCover = operationCertificatesCover
            && nodeCertificate.semanticError.rigorous
            && differential.absoluteError <= nodeCertificate.semanticError.upperBound;
        sawCertifiedRescale = sawCertifiedRescale
            || differential.operation == em::EvalModOperation::Rescale;
        sawCertifiedRelinearize = sawCertifiedRelinearize
            || differential.operation == em::EvalModOperation::Relinearize;
    }
    const m2424::CkksProfile shortProfile{32768, {60, 59, 60}, std::exp2(59.0), rawInput.size()};
    auto shortCandidate = candidate;
    const auto shortValidation = em::validateEvalModCandidateBackend(
        shortCandidate, problem, shortProfile, rawInput);
    auto prototype = synthesis.candidates[1];
    const auto prototypeValidation = em::validateEvalModCandidateBackend(
        prototype, problem, profile, rawInput);
    bool minimaxMatrixPassed = true;
    std::vector<em::EvalModBackendValidation> calibrationMeasurements;
    if (validation.executionSucceeded) calibrationMeasurements.push_back(validation);
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
            if (!measured.executionSucceeded)
                std::printf("minimax degree=%zu baby=%zu trial=%d failure=%s\n",
                            degree, babyStep, seedTrial, measured.failure.c_str());
            if (measured.executionSucceeded) calibrationMeasurements.push_back(measured);
            minimaxMatrixPassed = minimaxMatrixPassed && measured.executionSucceeded
                && minimax.backendRunnable && minimax.backendMeasured
                && measured.executedNodes == minimax.compiledCircuit.nodes.size();
        }
    }
    const auto calibratedModel = calibrationMeasurements.empty()
        ? em::EvalModArithmeticErrorModel{}
        : em::calibrateEvalModArithmeticModel(calibrationMeasurements, 4.0);

    const bool ok = !validation.passed && independentReferenceMatches
        && exactContext.dataPrimes == expectedDataPrimes
        && exactContext.specialPrime == expectedSpecialPrime
        && certifiedScaleSchedule.valid && certifiedScaleSchedule.rigorous
        && certifiedAlignments == 0
        && certifiedCircuit.cost.degree == 9
        && certifiedValidation.executionSucceeded
        && certifiedValidation.preparedConstantsUsed
        && certifiedValidation.runtimeScaleBitsMatch
        && certifiedValidation.executedNodes == certifiedCircuit.nodes.size()
        && certifiedCandidate.arithmeticCertificate.rigorous
        && certifiedCandidate.arithmeticCertificate.status
            == em::EvalModCertificationStatus::Certified
        && certifiedCandidate.arithmeticCertificate.keyNoiseMetadata.rigorous
        && sawCertifiedRescale && sawCertifiedRelinearize
        && operationCertificatesCover
        && !nonlinearPlan.rigorous
        && nonlinearPlan.status == em::EvalModCertificationStatus::ApproximationInsufficient
        && nonlinearPlan.securityBits >= 128
        && !planSearch.plan.has_value()
        && planSearch.status
            == em::EvalModPlanSearchStatus::NoCertifiedPlanInSearchSpace
        && !planSearch.globalImpossibilityProved
        && !planSearch.familiesSearched.empty()
        && planSearch.minimumDegree <= 7
        && planSearch.maximumDegree >= 15
        && std::isfinite(planSearch.bestApproximationBound)
        && std::isfinite(planSearch.bestArithmeticBound)
        && exactSchedule.available && exactSchedule.valid && exactSchedule.rigorous
        && exactSchedule.rescalePrimes[3] == expectedDataPrimes.back()
        && exactSchedule.nodeScales[3].numerator == (em::ExactInteger(1) << 118)
        && exactSchedule.nodeScales[3].denominator
            == em::ExactInteger(std::to_string(expectedDataPrimes.back()))
        && !prohibitedSchedule.valid && !prohibitedSchedule.rigorous
        && prohibitedSchedule.failingNode == 4
        && prohibitedSchedule.failure == "metadata_scale_alignment_prohibited"
        && modSwitchSchedule.valid && modSwitchSchedule.rigorous
        && headroom.available && headroom.valid && headroom.rigorous
        && headroom.modSwitchGates.size() == 1
        && headroom.modSwitchGates[0].proved && headroom.modSwitchGates[0].rigorous
        && headroom.modSwitchGates[0].encodedMagnitudeUpperBound > 0
        && !noHeadroom.valid && !noHeadroom.rigorous
        && noHeadroom.failingNode == 2
        && noHeadroom.failure == "modswitch_headroom_violation"
        && candidate.stage == em::EvalModCandidateStage::BackendMeasured
        && candidate.measuredBackendError.has_value() && candidate.executable
        && candidate.circuitValid && candidate.backendRunnable && candidate.backendMeasured
        && !candidate.rigorouslyValidated
        && candidate.intervalCertified && candidate.approximationCertified
        && !candidate.arithmeticErrorRigorous && !candidate.arithmeticErrorCertified
        && validation.executedNodes == candidate.compiledCircuit.nodes.size()
        && validation.preparedConstantsUsed
        && validation.runtimeScaleBitsMatch
        && !validation.firstRuntimeScaleMismatchNode.has_value()
        && validation.maxPreparedConstantEncodingError > 0.0
        && validation.maxPreparedConstantEncodingError < 1e-15
        && validation.executionSucceeded && validation.matchesPolynomialReference
        && !validation.matchesEvalModTarget
        && validation.implementationError <= problem.precisionBudget.implementation
        && validation.approximationError > problem.precisionBudget.approximation
        && validation.totalMeasuredError > problem.targetAbsoluteError
        && !validation.firstImplementationBudgetExceedingNode.has_value()
        && !validation.differentialTrace.empty()
        && validation.failure == "approximation_budget_exceeded"
        && validation.implementationError == validation.maxAbsoluteError
        && std::abs(std::log2(validation.outputScale) - 59.0) < 0.2
        && !prototypeValidation.passed && !prototype.executable
        && minimaxMatrixPassed
        && calibratedModel.calibrated && !calibratedModel.rigorous
        && calibratedModel.encodingAbsolute >= validation.implementationError
        && !shortValidation.passed && shortValidation.executedNodes == 0
        && shortValidation.failure == "profile_chain_too_short";
    if (!validation.matchesPolynomialReference) {
        for (const auto& node : validation.differentialTrace)
            std::printf("node[%zu]=%d chain=%zu scale=2^%.6f error=%.3e increase=%.3e\n",
                        node.node, static_cast<int>(node.operation), node.chainIndex,
                        std::log2(node.actualScale), node.absoluteError, node.errorIncrease);
    }
    if (!ok) {
        std::printf("diagnostic prepared=%d minimax=%d calibration=%d execution=%d polynomial=%d "
                    "target=%d implementation_budget=%d approximation_over=%d short_failure=%s\n",
                    validation.preparedConstantsUsed, minimaxMatrixPassed,
                    calibratedModel.calibrated, validation.executionSucceeded,
                    validation.matchesPolynomialReference, validation.matchesEvalModTarget,
                    validation.implementationError <= problem.precisionBudget.implementation,
                    validation.approximationError > problem.precisionBudget.approximation,
                    shortValidation.failure.c_str());
        std::printf("certified schedule valid=%d rigorous=%d alignments=%zu failure=%s nodes=%zu levels=%zu\n",
                    certifiedScaleSchedule.valid, certifiedScaleSchedule.rigorous,
                    certifiedAlignments, certifiedScaleSchedule.failure.c_str(),
                    certifiedCircuit.nodes.size(), certifiedCircuit.cost.levelConsumption);
    }
    std::printf("[test_evalmod_backend] error=%.3e certified_bound=%.3e plan_status=%d "
                "search_best=(%.3e,%.3e) "
                "prepared_constant_error=%.3e nodes=%zu failure=%s %s\n",
                validation.maxAbsoluteError,
                certifiedCandidate.arithmeticCertificate.outputError.upperBound,
                static_cast<int>(nonlinearPlan.status),
                planSearch.bestApproximationBound,
                planSearch.bestArithmeticBound,
                validation.maxPreparedConstantEncodingError,
                validation.executedNodes, validation.failure.c_str(),
                ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
