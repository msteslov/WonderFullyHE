#include "m2424/experimental/evalmod_analysis/approximation_lab.hpp"
#include "m2424/experimental/evalmod_analysis/cost_model.hpp"
#include "m2424/experimental/evalmod_analysis/exact_modular_oracle.hpp"

#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>

namespace em = m2424::experimental;

em::EvalModCiphertextModel model(std::size_t count, double sigma, double failure = -128.0) {
    return {em::TailModel::Subgaussian, count, 1.0, sigma, 0.1, 0.0004, 0.00005,
            0.00001, 0.00002, failure, 384,
            {"union bound from component subgaussian parameters",
             "encryption,key-switching,CtS", "independent centered component bounds"}};
}

int main() {
    const auto base = em::estimateEvalModDomain(model(4096, 0.25));
    const auto moreNoise = em::estimateEvalModDomain(model(4096, 0.5));
    const auto moreCoefficients = em::estimateEvalModDomain(model(8192, 0.25));
    const auto lowerFailure = em::estimateEvalModDomain(model(4096, 0.25, -192.0));
    const auto deterministic = em::estimateEvalModDomain({
        em::TailModel::Deterministic, 1, 2.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, -128.0, 256,
        {"deterministic bound", "none", "fixed integer offset"}
    });
    auto highKModel = model(4096, 0.75);
    const auto highK = em::estimateEvalModDomain(highKModel);
    auto periodLowModel = model(1, 0.0);
    periodLowModel.tailModel = em::TailModel::Deterministic;
    periodLowModel.deterministicIntegerOffset = 1.0;
    auto periodHighModel = periodLowModel;
    periodHighModel.deterministicIntegerOffset = 10.0;
    const auto periodLow = em::estimateEvalModDomain(periodLowModel);
    const auto periodHigh = em::estimateEvalModDomain(periodHighModel);

    const auto binaryScale = em::ExactScale::fromBinaryDouble(3.25);
    const auto oracle = em::exactCoefficientOracle({34, 23}, {101, 103}, 101,
                                                    em::ExactScale::rational(16, 1), 512);
    const auto positiveHalf = em::exactCoefficientOracle({50}, {101}, 101, binaryScale, 384);
    const auto negativeHalf = em::exactCoefficientOracle({51}, {101}, 101, binaryScale, 384);

    const em::EvalModDomain multiInterval{1, 0.2, "0.2", -128.0};
    const auto diagnostic = em::diagnoseEvalModPolynomialOnGrid(
        {em::PolynomialBasis::Monomial, {"0", "1"}}, multiInterval,
        17, "1e-400", 1e-30, 384, 1000);

    const auto cost = em::estimateEvalModCost(
        {31, 5, 8, 3, 4, 8, 7, 6, 4, 2048, 4096},
        {2.0, 0.4, 0.1, 0.5, 0.25, 1024}, 50);
    const auto lazyCost = em::estimateEvalModCost(
        {31, 9, 8, 10, 4, 8, 10, 3, 4, 0, 0},
        {2.0, 0.4, 0.1, 0.5, 0.25, 1024}, 50);

    bool discontinuityRejected = false;
    try {
        (void)em::estimateEvalModDomain({em::TailModel::Deterministic, 1, 0.0, 0.0,
            0.49, 0.0, 0.02, 0.0, 0.0, -128.0, 256,
            {"deterministic bound", "CtS", "absolute bounds"}});
    } catch (const std::domain_error&) { discontinuityRejected = true; }

    bool unprovedSigmaRejected = false;
    try {
        (void)em::estimateEvalModDomain({em::TailModel::Deterministic, 1, 0.0, 1.0,
            0.1, 0.0, 0.0, 0.0, 0.0, -128.0, 256,
            {"invalid deterministic sigma", "noise", "none"}});
    } catch (const std::invalid_argument&) { unprovedSigmaRejected = true; }

    bool nonDivisorRejected = false;
    try {
        (void)em::exactCoefficientOracle({1, 2}, {101, 103}, 97,
                                         em::ExactScale::rational(7, 1), 256);
    } catch (const std::invalid_argument&) { nonDivisorRejected = true; }

    bool nonCoprimeRejected = false;
    try {
        (void)em::exactCoefficientOracle({1, 2}, {15, 21}, 15,
                                         em::ExactScale::rational(7, 1), 256);
    } catch (const std::invalid_argument&) { nonCoprimeRejected = true; }

    bool nonFiniteCoefficientRejected = false;
    try {
        (void)em::diagnoseEvalModPolynomialOnGrid(
            {em::PolynomialBasis::Monomial, {"nan"}}, multiInterval, 2, "0.1", 0.0, 256, 100);
    } catch (const std::invalid_argument&) { nonFiniteCoefficientRejected = true; }
    bool nonFiniteRadiusRejected = false;
    try {
        (void)em::diagnoseEvalModPolynomialOnGrid(
            {em::PolynomialBasis::Monomial, {"0", "1"}}, multiInterval, 2, "inf", 0.0, 256, 100);
    } catch (const std::invalid_argument&) { nonFiniteRadiusRejected = true; }
    bool nonFiniteScaleRejected = false;
    try {
        (void)em::ExactScale::fromBinaryDouble(std::numeric_limits<double>::infinity());
    } catch (const std::invalid_argument&) { nonFiniteScaleRejected = true; }

    bool evaluationLimitRejected = false;
    try {
        (void)em::diagnoseEvalModPolynomialOnGrid(
            {em::PolynomialBasis::Monomial, {"0", "1"}}, multiInterval,
            100, "0.01", 0.0, 256, 100);
    } catch (const std::length_error&) { evaluationLimitRejected = true; }

    bool oversizedDomainRejected = false;
    try {
        const em::EvalModDomain oversized{
            static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max()), 0.1, "0.1", -128.0
        };
        (void)em::diagnoseEvalModPolynomialOnGrid(
            {em::PolynomialBasis::Monomial, {"0", "1"}}, oversized, 2, "0.01", 0.0, 256, 100);
    } catch (const std::invalid_argument&) { oversizedDomainRejected = true; }

    bool costOverflowRejected = false;
    try {
        auto overflowing = em::EvalModCircuitCost{1, 1, 1, 0, 0, 1, 1,
            std::numeric_limits<std::size_t>::max(), 1, 0, 0};
        (void)em::estimateEvalModCost(overflowing, {1, 1, 1, 1, 1, 1}, 2);
    } catch (const std::overflow_error&) { costOverflowRejected = true; }

    const double rawMargin = 0.5 - base.normalizedResidualBound;
    const double nearHalf = std::nextafter(0.5, 0.0);
    const em::EvalModDomain tinyMargin{0, nearHalf, "0.499999999999999944488848768742172978818416595458984375", -128.0};
    const bool ok = base.integerBound >= 1 && base.normalizedResidualBound > 0.1005
        && base.discontinuityMargin() < rawMargin && base.discontinuityMargin() > 0.0
        && moreNoise.integerBound >= base.integerBound
        && moreCoefficients.integerBound >= base.integerBound
        && lowerFailure.integerBound >= base.integerBound && deterministic.integerBound == 2
        && highK.integerBound > base.integerBound
        && highK.normalizedResidualBound > base.normalizedResidualBound
        && periodHigh.normalizedResidualBound > periodLow.normalizedResidualBound
        && em::isEvalModDomainValid(tinyMargin) && tinyMargin.discontinuityMargin() > 0.0
        && oracle.crtModulus == 10403 && oracle.centeredSourceCoefficient == 34
        && oracle.scaleNumerator == 16 && oracle.scaleDenominator == 1
        && oracle.precisionBits == 512 && oracle.expectedValueDecimal.find("2.125") == 0
        && oracle.roundingErrorAbsBound.mantissa == 0
        && positiveHalf.centeredSourceCoefficient == 50
        && positiveHalf.roundingErrorAbsBound.mantissa == 1
        && negativeHalf.centeredSourceCoefficient == -50
        && binaryScale.numerator == 13 && binaryScale.denominator == 4
        && diagnostic.approximationMaxError > 0.9
        && diagnostic.realDerivativeMax >= 1.0 && diagnostic.complexDerivativeMax >= 1.0
        && diagnostic.complexBoundaryErrorMax > 0.9 && diagnostic.evaluations == 255
        && std::abs(cost.latencyMs - 23.35) < 1e-12 && cost.peakWorkingSetBytes == 10240
        && cost.dataModulusBitsLowerBound == 350 && lazyCost.dataModulusBitsLowerBound == 500
        && discontinuityRejected && unprovedSigmaRejected && nonDivisorRejected
        && nonCoprimeRejected && nonFiniteCoefficientRejected && nonFiniteRadiusRejected
        && nonFiniteScaleRejected && evaluationLimitRejected
        && oversizedDomainRejected && costOverflowRejected;
    std::printf("[test_evalmod_contract] %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
