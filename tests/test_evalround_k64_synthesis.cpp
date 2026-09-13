#include "m2424/experimental/evalmod_analysis/evalround_synthesis.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <stdexcept>

using namespace m2424;
using namespace m2424::experimental;

void check(bool value, const char* why) { if (!value) throw std::runtime_error(why); }

int main() { try {
    // Exact binary64 value emitted by the selected sparse CtS certificate.
    const EvalRoundProblem problem{64, 5.2020386213659767e-5, 1e-2, 16};
    std::uint64_t rhoBits{};
    std::memcpy(&rhoBits, &problem.rho, sizeof(rhoBits));
    check(rhoBits == std::uint64_t{0x3f0b460edc2fc0f3ULL},
          "K=64 search uses the fixed sparse-CtS rho bits");
    check(evalRoundDigitCount(64, EvalRoundRadix::Binary) == 8,
          "K=64 must have eight offset-binary digits");
    const auto search = searchEvalRoundBinaryDigitPolynomials(problem);
    check(search.targetTable.size() == 129, "target table must contain every integer center");
    for (std::size_t row = 0; row < search.targetTable.size(); ++row) {
        check(search.targetTable[row].size() == 8, "target row has eight digits");
        check(search.targetTable[row] == evalRoundIntegerDigits(
            static_cast<std::int64_t>(row) - 64, 64, EvalRoundRadix::Binary),
            "target row is exactly bit_j(I+64)");
    }
    check(search.records.size() == 8 * 2 * 4, "every digit/family/degree is recorded");
    check(search.config.degrees == std::vector<std::size_t>({64, 128, 192, 256}),
          "bounded degrees remain exactly 64,128,192,256");
    std::vector<std::vector<bool>> remezAttempted(8, std::vector<bool>(4));
    bool sawConvergedRemez = false;
    for (const auto& record : search.records) {
        if (record.family == EvalModApproximationFamily::MultiIntervalMinimax) {
            const auto degree = std::find(search.config.degrees.begin(),
                search.config.degrees.end(), record.requestedDegree);
            check(degree != search.config.degrees.end(), "Remez attempted only configured degrees");
            remezAttempted[record.digitIndex][degree - search.config.degrees.begin()] = true;
            check(record.generatorStatus == EvalRoundDigitGeneratorStatus::Generated
                    && record.certificate && record.exchangeIterations > 0
                    && std::isfinite(record.rigorousIntervalError),
                  "direct target-agnostic Remez produces a separately certified candidate");
            check(record.exchangePointsInsideDomain,
                  "Remez exchange points never enter excluded gaps");
            check(record.cleanerInputDomainSatisfied
                    == (record.rigorousIntervalError <= 1),
                  "cleaner-domain flag is derived from the rigorous certificate");
            sawConvergedRemez = sawConvergedRemez || record.generatorConverged;
        }
        if (record.certificate) {
            const auto& proof = *record.certificate;
            check(proof.proof == EvalRoundPolynomialProof::OutwardInterval && proof.verified,
                  "accepted evidence is outward whole-domain evidence");
            check(proof.intervalSubdivisions == search.config.intervalSubdivisions,
                  "certificate records complete-cell subdivision coverage");
            check(proof.digitIndex == record.digitIndex && proof.certifiedK == 64,
                  "certificate is bound to digit and K");
            check(std::memcmp(&proof.certifiedRho, &problem.rho, sizeof(double)) == 0,
                  "certificate is bound to the actual rho bits");
        }
    }
    for (const auto& digits : remezAttempted) for (const bool attempted : digits)
        check(attempted, "all eight digits and four degrees are attempted by Remez");
    check(sawConvergedRemez && !search.extractorCertified,
          "Remez convergence alone cannot certify the extractor");
    for (std::size_t digit = 0; digit < 8; ++digit) {
        const auto& selected = search.records[*search.selectedRecordByDigit[digit]];
        for (const auto& record : search.records)
            if (record.digitIndex == digit && record.certificate
                && std::isfinite(record.rigorousIntervalError))
                check(selected.rigorousIntervalError <= record.rigorousIntervalError,
                      "selection uses the rigorous bound across Remez and diagnostics");
    }
    auto candidate = makeEvalRoundBinaryPolynomialCandidate(search);
    if (search.allDigitsCertified) {
        check(!search.allCleanerInputDomainsSatisfied && !search.extractorCertified,
              "finite outward bounds do not imply a cleaner-feasible extractor");
        check(!candidate.extraction.verified && candidate.extraction.polynomials.empty(),
              "cleaner-domain failure prevents candidate construction and planner input");
    } else {
        check(!candidate.extraction.verified && candidate.extraction.polynomials.empty()
                && search.status.find("not a global impossibility") != std::string::npos,
              "bounded incomplete search has explicit non-Certified provenance");
    }
    std::cout << search.status << '\n';
    std::cout << "PASS bounded K=64 target generation and certificate provenance\n";
} catch (const std::exception& error) {
    std::cerr << error.what() << '\n'; return 1;
}}
