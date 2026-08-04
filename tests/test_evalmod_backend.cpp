#include "m2424/experimental/evalmod_analysis/synthesis.hpp"

#include <cmath>
#include <cstdio>

namespace em = m2424::experimental;

int main() {
    const em::ExactInteger scaleInteger = em::ExactInteger(1) << 40;
    em::EvalModProblem problem{
        scaleInteger, em::ExactScale::rational(scaleInteger, 1),
        em::ExactScale::rational(scaleInteger, 1),
        {em::TailModel::Deterministic, 16384, 1.0, 0.0, 0.08, 0.0, 0.0, 0.0, 0.0,
         -128.0, 256, {"backend deterministic validation", "CKKS arithmetic", "bounded slots"}},
        8, 40, 1.0, {1.0, 0.2, 0.05, 0.2, 0.1, 4096}, 1.0, 0.0, 0.0
    };
    auto synthesis = em::synthesizeEvalMod(problem);
    auto selected = synthesis.candidates.end();
    for (auto it = synthesis.candidates.begin(); it != synthesis.candidates.end(); ++it) {
        if (it->family == em::EvalModApproximationFamily::MultiIntervalMinimax
            && it->satisfiesNumericalTarget && it->satisfiesLevelBudget) { selected = it; break; }
    }
    if (selected == synthesis.candidates.end()) return 1;
    auto& candidate = *selected;
    const std::vector<double> input{-1.06, -0.04, 0.0, 0.07, 1.03};
    const auto expected = em::executeEvalModDagPlaintext(candidate.compiledCircuit, input);
    const m2424::CkksProfile profile{16384, {50, 40, 40, 40, 40, 40, 40, 40, 50},
                                      std::exp2(40.0), input.size()};
    const auto validation = em::validateEvalModCandidateBackend(candidate, problem, profile,
                                                                 input, expected);
    const m2424::CkksProfile shortProfile{16384, {50, 40, 50}, std::exp2(40.0), input.size()};
    auto shortCandidate = candidate;
    shortCandidate.executable = false;
    const auto shortValidation = em::validateEvalModCandidateBackend(
        shortCandidate, problem, shortProfile, input, expected);
    auto prototype = synthesis.candidates[1];
    const auto prototypeValidation = em::validateEvalModCandidateBackend(
        prototype, problem, profile, input, expected);
    const bool ok = validation.passed && candidate.executable
        && candidate.stage == em::EvalModCandidateStage::BackendValidated
        && validation.executedNodes == candidate.compiledCircuit.nodes.size()
        && validation.maxAbsoluteError <= candidate.predictedBootstrapError * 1.25
        && std::abs(std::log2(validation.outputScale) - 40.0) < 0.1
        && !prototypeValidation.passed && !prototype.executable
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
