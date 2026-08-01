#include "m2424/experimental/evalmod_analysis/approximation_lab.hpp"
#include "m2424/experimental/evalmod_analysis/cost_model.hpp"
#include "m2424/experimental/evalmod_analysis/exact_modular_oracle.hpp"

#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>

int main() {
    const auto base = m2424::experimental::analyzeEvalModDomain({4096, 1.0, 0.25, 0.1, 0.0004, 0.00005, 0.00005, -128.0});
    const auto moreNoise = m2424::experimental::analyzeEvalModDomain({4096, 1.0, 0.5, 0.1, 0.0, 0.0, 0.0, -128.0});
    const auto moreCoefficients = m2424::experimental::analyzeEvalModDomain({8192, 1.0, 0.25, 0.1, 0.0, 0.0, 0.0, -128.0});
    const auto lowerFailure = m2424::experimental::analyzeEvalModDomain({4096, 1.0, 0.25, 0.1, 0.0, 0.0, 0.0, -192.0});
    const auto zeroSigma = m2424::experimental::analyzeEvalModDomain({1, 2.0, 0.0, 0.1, 0.0, 0.0, 0.0, -128.0});

    const auto oracle = m2424::experimental::exactCoefficientOracle({34, 23}, {101, 103}, 101,
                                                       "16", 512);
    const auto positiveHalf = m2424::experimental::exactCoefficientOracle({50}, {101}, 101, "3.25", 384);
    const auto negativeHalf = m2424::experimental::exactCoefficientOracle({51}, {101}, 101, "3.25", 384);

    const m2424::experimental::EvalModDomain multiInterval{1, 0.2, 0.3, -128.0};
    const auto diagnostic = m2424::experimental::diagnoseEvalModPolynomialOnGrid(
        {m2424::experimental::PolynomialBasis::Monomial, {"0", "1"}}, multiInterval,
        17, "0.01", 1e-30, 384, 1000);

    const auto cost = m2424::experimental::estimateEvalModCost(
        {31, 5, 8, 3, 4, 8, 7, 6, 4, 2048, 4096},
        {2.0, 0.4, 0.1, 0.5, 0.25, 1024}, 50);

    bool discontinuityRejected = false;
    try {
        (void)m2424::experimental::analyzeEvalModDomain({1, 0.0, 0.0, 0.49, 0.0, 0.02, 0.0, -128.0});
    } catch (const std::domain_error&) { discontinuityRejected = true; }

    bool nonDivisorRejected = false;
    try {
        (void)m2424::experimental::exactCoefficientOracle({1, 2}, {101, 103}, 97, "7", 256);
    } catch (const std::invalid_argument&) { nonDivisorRejected = true; }

    bool nonCoprimeRejected = false;
    try {
        (void)m2424::experimental::exactCoefficientOracle({1, 2}, {15, 21}, 15, "7", 256);
    } catch (const std::invalid_argument&) { nonCoprimeRejected = true; }

    bool evaluationLimitRejected = false;
    try {
        (void)m2424::experimental::diagnoseEvalModPolynomialOnGrid(
            {m2424::experimental::PolynomialBasis::Monomial, {"0", "1"}}, multiInterval,
            100, "0.01", 0.0, 256, 100);
    } catch (const std::length_error&) { evaluationLimitRejected = true; }

    bool oversizedDomainRejected = false;
    try {
        const m2424::experimental::EvalModDomain oversized{
            static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max()), 0.1, 0.4, -128.0
        };
        (void)m2424::experimental::diagnoseEvalModPolynomialOnGrid(
            {m2424::experimental::PolynomialBasis::Monomial, {"0", "1"}}, oversized,
            2, "0.01", 0.0, 256, 100);
    } catch (const std::invalid_argument&) { oversizedDomainRejected = true; }

    bool costOverflowRejected = false;
    try {
        auto overflowing = m2424::experimental::EvalModCircuitCost{1, 1, 1, 0, 0, 1, 1,
            std::numeric_limits<std::size_t>::max(), 1, 0, 0};
        (void)m2424::experimental::estimateEvalModCost(overflowing, {1, 1, 1, 1, 1, 1}, 2);
    } catch (const std::overflow_error&) { costOverflowRejected = true; }

    const bool ok = base.integerBound >= 1 && base.normalizedResidualBound > 0.1005
        && moreNoise.integerBound >= base.integerBound
        && moreCoefficients.integerBound >= base.integerBound
        && lowerFailure.integerBound >= base.integerBound && zeroSigma.integerBound == 2
        && oracle.crtModulus == 10403 && oracle.centeredSourceCoefficient == 34
        && oracle.precisionBits == 512 && oracle.expectedValueDecimal.find("2.125") == 0
        && positiveHalf.centeredSourceCoefficient == 50
        && negativeHalf.centeredSourceCoefficient == -50
        && diagnostic.approximationMaxError > 0.9
        && diagnostic.realDerivativeMax >= 1.0 && diagnostic.complexDerivativeMax >= 1.0
        && diagnostic.complexBoundaryErrorMax > 0.9 && diagnostic.evaluations == 255
        && std::abs(cost.latencyMs - 23.35) < 1e-12 && cost.peakWorkingSetBytes == 10240
        && cost.dataModulusBitsLowerBound == 350
        && discontinuityRejected && nonDivisorRejected && nonCoprimeRejected
        && evaluationLimitRejected && oversizedDomainRejected && costOverflowRejected;
    std::printf("[test_evalmod_contract] %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
