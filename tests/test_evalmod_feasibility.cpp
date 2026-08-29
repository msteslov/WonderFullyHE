#include "m2424/experimental/evalmod_analysis/feasibility.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

int main() {
    namespace em = m2424::experimental;
    const em::ExactInteger qSource = em::ExactInteger(1) << 60;
    const em::ExactInteger ctsScale = em::ExactInteger(1) << 59;
    const em::ExactInteger outputScale = em::ExactInteger(1) << 57;
    em::EvalModProblem problem{
        qSource, em::ExactScale::rational(ctsScale, 1),
        em::ExactScale::rational(outputScale, 1),
        {em::TailModel::Deterministic, 32768, 1.0, 0.0, 0.08, 0.0, 0.0,
         0.0, 0.0, -128.0, 512,
         {"deterministic feasibility regression", "CoeffToSlot and CKKS arithmetic",
          "bounded slots and exact binary scales"}},
        7, 59, 1e-10, {1.0, 0.2, 0.05, 0.2, 0.1, 4096},
        1.0, 0.0, 0.0};
    problem.precisionBudget = {1e-10, 1e-10};
    const auto synthesis = em::synthesizeEvalMod(problem);
    const auto selected = std::find_if(
        synthesis.candidates.begin(), synthesis.candidates.end(), [](const auto& candidate) {
            return candidate.family == em::EvalModApproximationFamily::PeriodicSineBaseline
                && candidate.approximationConverged;
        });
    if (selected == synthesis.candidates.end()) return 1;

    const m2424::CkksProfile profile{
        32768, {60, 59, 59, 59, 59, 59, 59, 59, 60},
        std::exp2(59.0), 8};
    auto adapter = m2424::SealAdapter::create(profile);
    adapter.generateKeys({}, true);
    const auto input = adapter.encrypt(adapter.encode(
        std::vector<double>{-2.16, -1.84, -0.16, 0.0, 0.16, 1.84, 2.16}));
    problem.exactModulusContext = em::EvalModExactModulusContext{
        adapter.coeffModulusValues(input), adapter.specialKeyModulusValue()};
    problem.requireCertifiedScaleSchedule = true;
    auto candidate = *selected;
    candidate.compiledCircuit = em::compileEvalModPolynomial(
        candidate.polynomial, problem, candidate.compiledCircuit.babyStep);
    const auto plan = em::prepareEvalMod(adapter, candidate, problem, input);
    if (!plan.arithmeticCertificate.rigorous) return 1;
    const auto report = em::analyzeEvalModFeasibility(problem, plan, &synthesis);

    const auto source = [&](const char* name) {
        return std::find_if(report.sources.begin(), report.sources.end(),
                            [&](const auto& item) { return item.source == name; });
    };
    const auto n8192 = std::find_if(
        report.profiles.begin(), report.profiles.end(), [](const auto& item) {
            return item.polyModulusDegree == 8192;
        });
    const bool allProfilesExcluded = std::all_of(
        report.profiles.begin(), report.profiles.end(), [](const auto& item) {
            return item.excludedByArithmetic || item.excludedByApproximation;
        });
    const double classified = plan.arithmeticCertificate.outputBreakdown.total();
    const bool breakdownCovers = classified + 1e-18
        >= plan.arithmeticCertificate.outputError.upperBound;
    const bool ok = report.status
            == em::EvalModFeasibilityStatus::CertificateClassInfeasible
        && !report.globalImpossibilityProved
        && report.integerBound == 1 && std::abs(report.rho - 0.08) < 1e-15
        && report.currentSecurityBits >= 128
        && report.tc128CoeffModulusBits == 881
        && report.currentCoeffModulusBits == 533
        && report.maximumSecurityValidDataPrimes == 13
        && report.maximumSecurityValidLevels == 12
        && report.maximumResourceFeasibleDegreeUpperBound == 41
        && report.mandatoryCurrentRescaleFloor > report.normalizedTarget
        && report.mandatoryRescaleImprovementFactor > 74.0
        && report.dominantSource == "Rescale"
        && report.requiredImprovementFactor > 1e5
        && report.optimisticKeySwitchPlusModDown < report.normalizedTarget
        && report.compositeBranchMargin > 0.4
        && report.quadraticApproximationLowerBound > 0.03
        && source("Input/CtS error") != report.sources.end()
        && source("Approximation error") != report.sources.end()
        && source("Coefficient quantization") != report.sources.end()
        && source("MultiplyPlain propagation") != report.sources.end()
        && source("CipherxCipher propagation") != report.sources.end()
        && source("Rescale") != report.sources.end()
        && source("Relinearization/KeySwitch") != report.sources.end()
        && source("ModSwitch") != report.sources.end()
        && source("Scale representation error") != report.sources.end()
        && source("Period mismatch") != report.sources.end()
        && source("Denormalization") != report.sources.end()
        && source("SlotToCoeff contribution") != report.sources.end()
        && source("Final additive error") != report.sources.end()
        && n8192 != report.profiles.end()
        && !n8192->excludedByArithmetic && n8192->excludedByApproximation
        && n8192->maximumLevels == 1
        && allProfilesExcluded && breakdownCovers
        && report.families.size() == 3
        && report.families[0].diagnosticMinimumDegree == 7
        && report.families[0].diagnosticMaximumDegree == 15
        && report.families[0].resourceDegreeUpperBound == 41
        && report.families[2].status == "PrunedByFeasibilityBeforeImplementation"
        && report.proofScope.find("not a global impossibility") != std::string::npos
        && report.remainingAlternatives.find("non-stationary") != std::string::npos;
    std::printf(
        "[test_evalmod_feasibility] status=%d target=%.3e rescale_floor=%.3e "
        "dominant=%.3e factor=%.3e max_levels=%zu quadratic_lb=%.3e %s\n",
        static_cast<int>(report.status), report.normalizedTarget,
        report.mandatoryCurrentRescaleFloor, report.dominantCurrentBound,
        report.requiredImprovementFactor, report.maximumSecurityValidLevels,
        report.quadraticApproximationLowerBound, ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
