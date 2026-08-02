#include "m2424/experimental/evalmod_analysis/synthesis.hpp"

#include <cmath>
#include <cstdio>

namespace em = m2424::experimental;

int main() {
    em::EvalModProblem problem{
        101, em::ExactScale::rational(16, 1), em::ExactScale::rational(16, 1),
        {em::TailModel::Deterministic, 16, 1.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0,
         -128.0, 256, {"deterministic validation profile", "none", "bounded coefficients"}},
        12, 40, 1.0, {1.0, 0.2, 0.05, 0.2, 0.1, 4096}, 1.0, 0.0, 0.0
    };
    const auto result = em::synthesizeEvalMod(problem);
    const auto csv = em::evalModSynthesisCsv(result);
    const auto json = em::evalModSynthesisJson(result);
    bool plaintextMatches = true;
    const auto& polynomial = result.candidates[1].polynomial;
    for (int integer = -1; integer <= 1; ++integer) {
        for (int message : {-8, -4, 0, 4, 8}) {
            const double residual = static_cast<double>(message) / 101.0;
            const std::uint64_t residue = message < 0 ? static_cast<std::uint64_t>(101 + message)
                                                      : static_cast<std::uint64_t>(message);
            const auto oracle = em::exactCoefficientOracle({residue}, {101}, 101,
                                                            em::ExactScale::rational(16, 1), 256);
            const double expected = std::stod(oracle.expectedValueDecimal);
            plaintextMatches = plaintextMatches
                && std::abs(em::evaluateEvalModPolynomial(polynomial, integer + residual)
                                * (101.0 / 16.0) - expected) < 0.01;
        }
    }
    const bool ok = result.domain.integerBound == 1 && result.candidates.size() == 2
        && result.selectedCandidate.has_value() && result.candidates[0].circuit.degree == 9
        && result.candidates[1].circuit.degree == 15
        && !result.candidates[0].scaleSchedule.empty()
        && result.candidates[1].diagnostic.evaluations > 0
        && csv.find("periodic_sine_baseline") != std::string::npos
        && json.find("multi_interval_least_squares_prototype") != std::string::npos
        && json.find("tail_model_provenance") != std::string::npos && plaintextMatches;
    std::printf("[test_evalmod_synthesis] %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
