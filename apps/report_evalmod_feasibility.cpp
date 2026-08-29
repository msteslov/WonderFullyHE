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
         {"deterministic feasibility domain", "CoeffToSlot and CKKS arithmetic",
          "bounded slots and exact binary scales"}},
        7, 59, 1e-10, {1.0, 0.2, 0.05, 0.2, 0.1, 4096},
        1.0, 0.0, 0.0};
    problem.precisionBudget = {1e-10, 1e-10};

    const m2424::CkksProfile profile{
        32768, {60, 59, 59, 59, 59, 59, 59, 59, 60},
        std::exp2(59.0), 8};
    auto adapter = m2424::SealAdapter::create(profile);
    adapter.generateKeys({}, true);
    const auto input = adapter.encrypt(adapter.encode(
        std::vector<double>{-2.16, -1.84, -0.16, 0.0, 0.16, 1.84, 2.16}));
    const auto synthesis = em::synthesizeEvalMod(problem);
    const auto selected = std::find_if(
        synthesis.candidates.begin(), synthesis.candidates.end(), [](const auto& candidate) {
            return candidate.family == em::EvalModApproximationFamily::PeriodicSineBaseline
                && candidate.approximationConverged;
        });
    if (selected == synthesis.candidates.end()) return 2;
    auto candidate = *selected;
    problem.exactModulusContext = em::EvalModExactModulusContext{
        adapter.coeffModulusValues(input), adapter.specialKeyModulusValue()};
    problem.requireCertifiedScaleSchedule = true;
    candidate.compiledCircuit = em::compileEvalModPolynomial(
        candidate.polynomial, problem, candidate.compiledCircuit.babyStep);
    const auto plan = em::prepareEvalMod(adapter, candidate, problem, input);
    if (!plan.arithmeticCertificate.rigorous) return 3;
    const auto report = em::analyzeEvalModFeasibility(problem, plan, &synthesis);
    const auto json = em::evalModFeasibilityJson(report);
    std::fwrite(json.data(), 1, json.size(), stdout);
    return report.status == em::EvalModFeasibilityStatus::CertificateClassInfeasible ? 0 : 1;
}
