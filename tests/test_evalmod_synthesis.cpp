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
    auto strictProblem = problem;
    strictProblem.targetAbsoluteError = 1e-10;
    const auto strict = em::synthesizeEvalMod(strictProblem);
    auto shallowProblem = problem;
    shallowProblem.availableLevels = 1;
    const auto shallow = em::synthesizeEvalMod(shallowProblem);
    auto changedCts = problem;
    changedCts.coeffToSlotScale = em::ExactScale::rational(32, 1);
    const auto changedCtsCircuit = em::compileEvalModPolynomial(
        result.candidates[0].polynomial, changedCts, 3);
    auto changedOutput = problem;
    changedOutput.outputScale = em::ExactScale::rational(32, 1);
    const auto changedOutputCircuit = em::compileEvalModPolynomial(
        result.candidates[0].polynomial, changedOutput, 3);
    auto profile8192 = problem;
    profile8192.coeffToSlotScale = em::ExactScale::fromBinaryDouble(std::exp2(59.5));
    profile8192.outputScale = profile8192.coeffToSlotScale;
    profile8192.ciphertextModel.coefficientCount = 16384;
    profile8192.ciphertextModel.provenance.derivation = "N=16384, slots=8192 validation model";
    const auto profileResult = em::synthesizeEvalMod(profile8192);
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
    const auto& compiled = result.candidates[0].compiledCircuit;
    auto invalidTypedCircuit = compiled;
    for (auto& node : invalidTypedCircuit.nodes) {
        if (node.operation == em::EvalModOperation::ModSwitch) {
            node.inputs[0] = 1; // EncodeConstant cannot feed ciphertext alignment.
            break;
        }
    }
    const std::vector<double> dagRawInput{-6.5, -3.0, 0.0, 3.0, 6.5};
    const auto dagOutput = em::executeEvalModDagPlaintext(compiled, dagRawInput);
    bool dagMatchesReference = dagOutput.size() == dagRawInput.size();
    for (std::size_t index = 0; index < dagRawInput.size(); ++index) {
        const double z = dagRawInput[index] * (16.0 / 101.0);
        const double expected = em::evaluateEvalModPolynomial(result.candidates[0].polynomial, z)
            * (101.0 / 16.0);
        dagMatchesReference = dagMatchesReference && std::abs(dagOutput[index] - expected) < 1e-10;
    }
    std::size_t ctCt = 0, ctPt = 0, relin = 0, rescales = 0, additions = 0;
    std::size_t modSwitches = 0, scaleAlignments = 0, plaintextAdditions = 0;
    for (const auto& node : compiled.nodes) {
        ctCt += node.operation == em::EvalModOperation::MultiplyCipher;
        ctPt += node.operation == em::EvalModOperation::MultiplyPlain;
        relin += node.operation == em::EvalModOperation::Relinearize;
        rescales += node.operation == em::EvalModOperation::Rescale;
        additions += node.operation == em::EvalModOperation::Add
            || node.operation == em::EvalModOperation::AddPlain;
        modSwitches += node.operation == em::EvalModOperation::ModSwitch;
        scaleAlignments += node.operation == em::EvalModOperation::AlignScale;
        plaintextAdditions += node.operation == em::EvalModOperation::AddPlain;
    }
    bool prototypeSelected = result.provisionalSelection
        && result.candidates[*result.provisionalSelection].family
            == em::EvalModApproximationFamily::MultiIntervalLeastSquaresPrototype;
    const bool ok = result.domain.integerBound == 1 && result.candidates.size() >= 22
        && result.provisionalSelection.has_value()
        && result.candidates[0].compiledCircuit.cost.degree == 9
        && result.candidates[1].compiledCircuit.cost.degree == 15
        && !result.candidates[0].scaleSchedule.empty()
        && result.candidates[1].diagnostic.evaluations > 0
        && result.candidates[0].intervalCertified
        && result.candidates[0].intervalCertificate.proved && result.candidates[0].executable
        && result.candidates[0].polynomialArithmeticError.has_value()
        && !result.candidates[0].arithmeticErrorRigorous
        && !prototypeSelected && !strict.provisionalSelection.has_value()
        && !shallow.candidates[0].satisfiesLevelBudget
        && compiled.cost.ciphertextMultiplications == ctCt
        && compiled.cost.ciphertextPlaintextMultiplications == ctPt
        && compiled.cost.relinearizations == relin && compiled.cost.rescales == rescales
        && compiled.cost.additions == additions
        && compiled.cost.modulusSwitches == modSwitches
        && compiled.cost.scaleAlignments == scaleAlignments
        && compiled.cost.plaintextAdditions == plaintextAdditions
        && em::isCompiledEvalModCircuitValid(compiled)
        && !em::isCompiledEvalModCircuitValid(invalidTypedCircuit)
        && dagMatchesReference
        && changedCtsCircuit.normalizationGain.numerator
            != compiled.normalizationGain.numerator
        && changedOutputCircuit.denormalizationGain.denominator
            != compiled.denormalizationGain.denominator
        && profileResult.domain.integerBound >= 1 && profileResult.candidates.size() >= 22
        && csv.find("periodic_sine_baseline") != std::string::npos
        && json.find("multi_interval_least_squares_prototype") != std::string::npos
        && json.find("\"arithmetic_error_rigorous\":false") != std::string::npos
        && json.find("\"measured_backend_error\":null") != std::string::npos
        && json.find("tail_model_provenance") != std::string::npos && plaintextMatches;
    std::printf("[test_evalmod_synthesis] %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
