#include "m2424/experimental/evalmod_analysis/evalround_synthesis.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>

using namespace m2424;
using namespace m2424::experimental;

int main() {
    // Exact binary64 value emitted by test_sparse_coeff_to_slot for the
    // unchanged N=16384, seven-60-bit-prime sparse certificate.
    const EvalRoundProblem problem{64, 5.2020386213659767e-5, 1e-2, 16};
    const auto result = searchEvalRoundBinaryDigitPolynomials(problem);
    std::cout << std::setprecision(17);
    std::cout << "scope K=" << problem.K << " rho=" << problem.rho
              << " intervals=" << result.targetTable.size() << " degrees=";
    for (const auto degree : result.config.degrees) std::cout << degree << ',';
    std::cout << "\nstatus=" << result.status << '\n';
    std::cout << "all_digit_interval_certificates=" << result.allDigitsCertified
              << " cleaner_domains=" << result.allCleanerInputDomainsSatisfied
              << " extractor_certified=" << result.extractorCertified << '\n';
    for (std::size_t digit = 0; digit < 8; ++digit) {
        const EvalRoundDigitSearchRecord* best = nullptr;
        for (const auto& record : result.records) if (record.digitIndex == digit
            && record.certificate && std::isfinite(record.rigorousIntervalError)
            && (!best || record.rigorousIntervalError < best->rigorousIntervalError)) best = &record;
        std::cout << "digit=" << digit;
        if (best) std::cout << " family=multi_interval_chebyshev degree=" << best->requestedDegree
            << " generated=" << best->generatorConverged
            << " grid=" << best->gridMaximumError
            << " rigorous=" << best->rigorousIntervalError
            << " max_coefficient=" << best->maximumCoefficientMagnitude
            << " cleaner_a_le_1=" << best->cleanerInputDomainSatisfied;
        else std::cout << " no_finite_interval_candidate";
        std::cout << '\n';
    }
    std::cout << "minimax_status=not_applicable: existing remezOdd target is nearest-integer residual\n";
    for (const auto& record : result.records) {
        if (record.family != EvalModApproximationFamily::MultiIntervalChebyshev) continue;
        std::cout << "record digit=" << record.digitIndex
                  << " degree=" << record.requestedDegree
                  << " status=" << static_cast<int>(record.generatorStatus)
                  << " converged=" << record.generatorConverged
                  << " grid=" << record.gridMaximumError
                  << " rigorous=" << record.rigorousIntervalError
                  << " max_coefficient=" << record.maximumCoefficientMagnitude
                  << " cleaner_a_le_1=" << record.cleanerInputDomainSatisfied << '\n';
    }
    const auto candidate = makeEvalRoundBinaryPolynomialCandidate(result);
    for (const double budget : {1e-2, 1e-4, 1e-6}) {
        auto diagnostic = problem; diagnostic.requiredIntegerError = budget;
        const auto plan = planEvalRoundCandidate(diagnostic, candidate);
        std::cout << "cleaning budget=" << budget
                  << " status=" << static_cast<int>(plan.status)
                  << " rejection=" << static_cast<int>(plan.rejection)
                  << " iterations=" << plan.totalCleaningIterations
                  << " final_E_I=" << plan.integerErrorUpper
                  << " cleaner_domain=false contributions=";
        for (std::size_t digit = 0; digit < 8; ++digit) {
            const auto& selected = result.records[*result.selectedRecordByDigit[digit]];
            std::cout << std::ldexp(selected.rigorousIntervalError,
                                    static_cast<int>(digit)) << ',';
        }
        std::cout << '\n';
    }
    std::cout << "ciphertext_compile=attempted_in_test; first_blocker=ExtractionNotCertified because mathematical planning stops at CleaningDomainViolation; no ciphertext depth claim\n";
}
