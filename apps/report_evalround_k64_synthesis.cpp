#include "m2424/experimental/evalmod_analysis/evalround_synthesis.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>

using namespace m2424;
using namespace m2424::experimental;

const char* familyName(EvalModApproximationFamily family) {
    return family == EvalModApproximationFamily::MultiIntervalMinimax
        ? "multi_interval_minimax" : "multi_interval_chebyshev";
}

const char* proofName(EvalRoundIntervalProofMethod method) {
    return method == EvalRoundIntervalProofMethod::CenteredShiftHorner
        ? "centered_shift" : method == EvalRoundIntervalProofMethod::DirectXHorner
            ? "direct_x" : "unknown";
}

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
        if (best) std::cout << " family=" << familyName(best->family)
            << " degree=" << best->requestedDegree
            << " generated=" << best->generatorConverged
            << " grid=" << best->gridMaximumError
            << " direct_x=" << best->directXRigorousIntervalError
            << " centered=" << best->centeredRigorousIntervalError
            << " rigorous=" << best->rigorousIntervalError
            << " selected_proof=" << proofName(best->selectedIntervalProofMethod)
            << " max_coefficient=" << best->maximumCoefficientMagnitude
            << " cleaner_a_le_1=" << best->cleanerInputDomainSatisfied;
        else std::cout << " no_finite_interval_candidate";
        std::cout << '\n';
    }
    for (const auto& record : result.records) {
        std::cout << "record family=" << familyName(record.family)
                  << " digit=" << record.digitIndex
                  << " degree=" << record.requestedDegree
                  << " status=" << static_cast<int>(record.generatorStatus)
                  << " converged=" << record.generatorConverged
                  << " iterations=" << record.exchangeIterations
                  << " exchange_inside_domain=" << record.exchangePointsInsideDomain
                  << " grid=" << record.gridMaximumError
                  << " direct_x=" << record.directXRigorousIntervalError
                  << " centered=" << record.centeredRigorousIntervalError
                  << " rigorous=" << record.rigorousIntervalError
                  << " direct_over_centered="
                  << record.directXRigorousIntervalError/record.centeredRigorousIntervalError
                  << " selected_proof=" << proofName(record.selectedIntervalProofMethod)
                  << " max_coefficient=" << record.maximumCoefficientMagnitude
                  << " cleaner_a_le_1=" << record.cleanerInputDomainSatisfied << '\n';
    }
    const auto candidate = makeEvalRoundBinaryPolynomialCandidate(result);
    bool mathematicalPlanCertified = false;
    for (const double budget : {1e-2, 1e-4, 1e-6}) {
        auto diagnostic = problem; diagnostic.requiredIntegerError = budget;
        const auto plan = planEvalRoundCandidate(diagnostic, candidate);
        mathematicalPlanCertified = mathematicalPlanCertified
            || plan.status == EvalRoundPlanStatus::Certified;
        std::cout << "cleaning budget=" << budget
                  << " status=" << static_cast<int>(plan.status)
                  << " rejection=" << static_cast<int>(plan.rejection)
                  << " rounds=";
        for (const auto& digit : plan.digits) std::cout << digit.cleaningIterations << ',';
        std::cout << " trajectories=";
        for (std::size_t digit = 0; digit < plan.digits.size(); ++digit) {
            std::cout << digit << ':';
            for (const auto error : plan.digits[digit].errorAfterRounds)
                std::cout << error << '/';
            std::cout << " weighted=" << plan.digits[digit].weightedReconstructionError << ';';
        }
        std::cout << " final_E_I=" << plan.integerErrorUpper << '\n';
    }
    std::cout << "ciphertext_compile="
              << (mathematicalPlanCertified ? "mathematically_permitted; see analysis test for backend gate" : "not_attempted")
              << "; mathematical_plan_certified=" << mathematicalPlanCertified << '\n';
}
