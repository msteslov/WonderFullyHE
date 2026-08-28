#include "m2424/experimental/evalmod_analysis/synthesis.hpp"

#include <cmath>
#include <cstdio>
#include <limits>

namespace em = m2424::experimental;

namespace {

em::EvalModDagNode node(em::EvalModOperation operation,
                        std::vector<std::size_t> inputs,
                        std::size_t chain,
                        em::ExactScale scale,
                        std::string decimal = {}) {
    return {operation, std::move(inputs), chain, scale, scale, std::move(decimal), 0.0, 0.0};
}

em::CompiledEvalModCircuit linearCircuit() {
    const auto scale16 = em::ExactScale::rational(16, 1);
    const auto scale1 = em::ExactScale::rational(1, 1);
    em::CompiledEvalModCircuit circuit;
    circuit.nodes = {
        node(em::EvalModOperation::Input, {}, 0, scale16),
        node(em::EvalModOperation::EncodeConstant, {}, 0, scale1, "2"),
        node(em::EvalModOperation::MultiplyPlain, {0, 1}, 0, scale16),
        node(em::EvalModOperation::EncodeConstant, {}, 0, scale16, "3"),
        node(em::EvalModOperation::AddPlain, {2, 3}, 0, scale16)
    };
    circuit.outputNode = 4;
    circuit.normalizationGain = scale1;
    circuit.denormalizationGain = scale1;
    return circuit;
}

em::PreparedEvalModConstants constants(const std::vector<std::uint64_t>& primes) {
    em::PreparedEvalModConstants result;
    result.inputDataPrimes = primes;
    result.specialPrime = 113;
    result.polyModulusDegree = 8;
    result.secretCoefficientAbsSupport = 1;
    result.evaluationKeyNoiseCoefficientAbsSupport = 21;
    result.deterministicEvaluationKeyNoiseSupport = true;
    result.rigorous = true;
    em::PreparedEvalModConstant two;
    two.node = 1;
    two.decimal = "2";
    two.encodingScale = em::ExactScale::rational(1, 1);
    two.encodingErrorUpperBoundDecimal = "0.1";
    two.encodingErrorUpperBound = 0.1;
    two.rigorous = true;
    result.constants.push_back(std::move(two));
    em::PreparedEvalModConstant three;
    three.node = 3;
    three.decimal = "3";
    three.encodingScale = em::ExactScale::rational(16, 1);
    three.encodingErrorUpperBoundDecimal = "0.2";
    three.encodingErrorUpperBound = 0.2;
    three.rigorous = true;
    result.constants.push_back(std::move(three));
    return result;
}

} // namespace

int main() {
    const auto convertedChebyshev = em::convertChebyshevToMonomial(
        {em::PolynomialBasis::Chebyshev, {"1", "2", "3"}});
    const em::EvalModExactModulusContext context{{97, 193}, 113};
    const auto linear = linearCircuit();
    const auto linearSchedule = em::buildExactEvalModScaleSchedule(linear, context, true);
    const auto prepared = constants(context.dataPrimes);
    const auto linearCertificate = em::certifyEvalModDagArithmetic(
        linear, linearSchedule, prepared, 4.0, 0.5);

    em::CompiledEvalModCircuit multiply;
    const auto scale16 = em::ExactScale::rational(16, 1);
    const auto scale256 = em::ExactScale::rational(256, 1);
    multiply.nodes = {
        node(em::EvalModOperation::Input, {}, 0, scale16),
        node(em::EvalModOperation::MultiplyCipher, {0, 0}, 0, scale256),
        node(em::EvalModOperation::Relinearize, {1}, 0, scale256)
    };
    multiply.outputNode = 1;
    multiply.normalizationGain = em::ExactScale::rational(1, 1);
    multiply.denormalizationGain = em::ExactScale::rational(1, 1);
    const auto multiplySchedule = em::buildExactEvalModScaleSchedule(multiply, context, true);
    em::PreparedEvalModConstants noConstants;
    noConstants.inputDataPrimes = context.dataPrimes;
    noConstants.specialPrime = context.specialPrime;
    noConstants.polyModulusDegree = 8;
    noConstants.secretCoefficientAbsSupport = 1;
    noConstants.evaluationKeyNoiseCoefficientAbsSupport = 21;
    noConstants.deterministicEvaluationKeyNoiseSupport = true;
    noConstants.rigorous = true;
    const auto multiplyCertificate = em::certifyEvalModDagArithmetic(
        multiply, multiplySchedule, noConstants, 4.0, 0.5);

    multiply.outputNode = 2;
    const auto relinCertificate = em::certifyEvalModDagArithmetic(
        multiply, multiplySchedule, noConstants, 4.0, 0.5);
    auto unsupportedNoise = noConstants;
    unsupportedNoise.deterministicEvaluationKeyNoiseSupport = false;
    const auto unsupportedRelinCertificate = em::certifyEvalModDagArithmetic(
        multiply, multiplySchedule, unsupportedNoise, 4.0, 0.5);

    em::CompiledEvalModCircuit rescale;
    rescale.nodes = {
        node(em::EvalModOperation::Input, {}, 0, scale256),
        node(em::EvalModOperation::Rescale, {0}, 1,
             em::ExactScale::rational(256, 193))
    };
    rescale.outputNode = 1;
    rescale.normalizationGain = em::ExactScale::rational(1, 1);
    rescale.denormalizationGain = em::ExactScale::rational(1, 1);
    const auto rescaleSchedule = em::buildExactEvalModScaleSchedule(rescale, context, true);
    const auto rescaleCertificate = em::certifyEvalModDagArithmetic(
        rescale, rescaleSchedule, noConstants, 4.0, 0.5);

    const m2424::CkksProfile runtimeProfile{
        8192, {60, 40, 60}, std::exp2(40.0), 4};
    auto adapter = m2424::SealAdapter::create(runtimeProfile);
    adapter.generateKeys(std::vector<int>{}, false);
    em::EvalModProblem planProblem{
        em::ExactInteger(1) << 40,
        em::ExactScale::rational(em::ExactInteger(1) << 40, 1),
        em::ExactScale::rational(em::ExactInteger(1) << 40, 1),
        {em::TailModel::Deterministic, 8192, 0.0, 0.0, 0.08,
         0.0, 0.0, 0.0, 0.0, -128.0, 256,
         {"deterministic K=0 prepared-plan regression", "none", "bounded real slots"}},
        2, 40, 1.0, {1.0, 0.2, 0.05, 0.2, 0.1, 4096},
        1.0, 0.0, 0.0};
    planProblem.exactModulusContext = em::EvalModExactModulusContext{
        adapter.dataModulusValues(), adapter.specialKeyModulusValue()};
    planProblem.requireCertifiedScaleSchedule = true;
    em::EvalModCandidate planCandidate;
    planCandidate.domain = em::estimateEvalModDomain(planProblem.ciphertextModel);
    planCandidate.polynomial = {em::PolynomialBasis::Monomial, {"0", "1"}};
    planCandidate.compiledCircuit = em::compileEvalModPolynomial(
        planCandidate.polynomial, planProblem, 2);
    planCandidate.approximationConverged = true;
    planCandidate.intervalCertificate = em::certifyEvalModPolynomialIntervals(
        planCandidate.polynomial, planCandidate.domain, "1e-30", 256, 384);
    planCandidate.intervalCertified = planCandidate.intervalCertificate.proved;
    const std::vector<double> planInput{-0.08, -0.02, 0.03, 0.08};
    const auto planCipher = adapter.encrypt(adapter.encode(planInput));
    const auto preparedPlan = em::prepareEvalMod(
        adapter, planCandidate, planProblem, planCipher);
    auto noBudgetProblem = planProblem;
    noBudgetProblem.targetAbsoluteError = 1e-12;
    noBudgetProblem.slotToCoeffAdditive = 1e-6;
    const auto noBudgetPlan = em::prepareEvalMod(
        adapter, planCandidate, noBudgetProblem, planCipher);
    auto missingContextProblem = planProblem;
    missingContextProblem.exactModulusContext.reset();
    bool missingContextRejected = false;
    try {
        (void)em::compileEvalModPolynomial(
            planCandidate.polynomial, missingContextProblem, 2);
    } catch (const std::invalid_argument& error) {
        missingContextRejected = std::string(error.what()) == "MissingExactModulusContext";
    }
    const auto preflight = em::preflightEvalMod(adapter, preparedPlan, planCipher);
    em::EvalModExecutionTrace planTrace;
    const auto planOutput = em::applyEvalMod(adapter, preparedPlan, planCipher, &planTrace);
    const auto planDecoded = adapter.decode(adapter.decrypt(planOutput));
    auto wrongScaleCipher = adapter.normalizeScale(
        planCipher, std::nextafter(adapter.scale(planCipher),
                                   std::numeric_limits<double>::infinity()));
    const auto wrongScalePreflight = em::preflightEvalMod(
        adapter, preparedPlan, wrongScaleCipher);
    auto publicEvaluator = m2424::SealAdapter::create(runtimeProfile);
    publicEvaluator.loadPublicKey(adapter.savePublicKey());
    const auto publicInput = publicEvaluator.loadCipher(adapter.saveCipher(planCipher));
    const auto publicPreflight = em::preflightEvalMod(
        publicEvaluator, preparedPlan, publicInput);
    const auto publicOutput = em::applyEvalMod(
        publicEvaluator, preparedPlan, publicInput);
    bool publicEvaluatorHasNoSecret = false;
    try {
        (void)publicEvaluator.decrypt(publicOutput);
    } catch (const std::runtime_error&) {
        publicEvaluatorHasNoSecret = true;
    }
    double planError = 0.0;
    for (std::size_t index = 0; index < planInput.size(); ++index)
        planError = std::max(planError, std::abs(planDecoded[index] - planInput[index]));

    const double expectedMultiplyError = 4.0 * 0.5 + 4.0 * 0.5 + 0.5 * 0.5;
    const double expectedPlainError = 2.0 * 0.5 + 4.0 * 0.1 + 0.5 * 0.1;
    const bool ok = linearSchedule.valid && linearSchedule.rigorous
        && convertedChebyshev.basis == em::PolynomialBasis::Monomial
        && convertedChebyshev.decimalCoefficients
            == std::vector<std::string>({"-2", "2", "6"})
        && linearCertificate.rigorous
        && linearCertificate.status == em::EvalModCertificationStatus::Certified
        && linearCertificate.log2FailureProbability
            == -std::numeric_limits<double>::infinity()
        && linearCertificate.nodes[2].valueAbs.upperBound >= 8.0
        && linearCertificate.nodes[2].semanticError.upperBound >= expectedPlainError
        && linearCertificate.outputError.upperBound >= expectedPlainError + 0.2
        && linearCertificate.nodes[2].localAddedError.upperBound == 0.0
        && linearCertificate.nodes[4].localAddedError.upperBound == 0.0
        && multiplySchedule.valid && multiplySchedule.rigorous
        && multiplyCertificate.rigorous
        && multiplyCertificate.outputError.upperBound >= expectedMultiplyError
        && multiplyCertificate.nodes[1].localAddedError.upperBound == 0.0
        && relinCertificate.rigorous
        && relinCertificate.status == em::EvalModCertificationStatus::Certified
        && relinCertificate.nodes[2].localAddedError.upperBound > 0.0
        && relinCertificate.keyNoiseMetadata.rigorous
        && relinCertificate.keyNoiseMetadata.relevantEvaluationKeyComponents == 1
        && !unsupportedRelinCertificate.rigorous
        && unsupportedRelinCertificate.status
            == em::EvalModCertificationStatus::RigorousKeySwitchBoundUnavailable
        && unsupportedRelinCertificate.failingNode == 2
        && rescaleCertificate.rigorous
        && rescaleCertificate.status == em::EvalModCertificationStatus::Certified
        && rescaleCertificate.nodes[1].localAddedError.upperBound > 0.0
        && preparedPlan.rigorous
        && preparedPlan.status == em::EvalModCertificationStatus::Certified
        && preparedPlan.securityBits >= 128
        && !noBudgetPlan.rigorous
        && noBudgetPlan.status
            == em::EvalModCertificationStatus::NoEvalModErrorBudgetRemaining
        && missingContextRejected
        && preflight.compatible
        && planTrace.levelsConsumed == 0
        && planError < 1e-6
        && !wrongScalePreflight.compatible
        && wrongScalePreflight.status == em::EvalModCertificationStatus::InputScaleMismatch
        && publicPreflight.compatible
        && publicEvaluatorHasNoSecret;
    std::printf("[test_evalmod_certificates] linear=%.17g multiply=%.17g relin=%d "
                "rescale=%d plan=%.3e %s\n",
                linearCertificate.outputError.upperBound,
                multiplyCertificate.outputError.upperBound,
                static_cast<int>(relinCertificate.status),
                static_cast<int>(rescaleCertificate.status), planError,
                ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
