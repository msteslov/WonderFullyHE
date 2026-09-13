#include "m2424/experimental/evalmod_analysis/evalround_synthesis.hpp"
#include "bootstrap_fixture.hpp"

#include <cmath>
#include <cstring>
#include <iostream>
#include <stdexcept>

using namespace m2424;
using namespace m2424::experimental;

void check(bool value, const char* why) { if (!value) throw std::runtime_error(why); }

int main() { try {
    // Exact binary64 value emitted by the selected sparse CtS certificate.
    const EvalRoundProblem problem{64, 5.2020386213659767e-5, 1e-2, 16};
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
    for (const auto& record : search.records) {
        if (record.family == EvalModApproximationFamily::MultiIntervalMinimax) {
            check(record.generatorStatus == EvalRoundDigitGeneratorStatus::NotApplicableToTarget
                && !record.certificate, "target-locked Remez output is never relabelled as a digit proof");
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
    auto candidate = makeEvalRoundBinaryPolynomialCandidate(search);
    if (search.allDigitsCertified) {
        check(candidate.extraction.polynomials.size() == 8 && candidate.extraction.verified,
              "complete result contains all eight certificates");
        check(!search.allCleanerInputDomainsSatisfied && !search.extractorCertified,
              "finite outward bounds do not imply a cleaner-feasible extractor");
        for (const double budget : {1e-2, 1e-4, 1e-6}) {
            auto diagnostic = problem; diagnostic.requiredIntegerError = budget;
            const auto plan = planEvalRoundCandidate(diagnostic, candidate);
            check(plan.status != EvalRoundPlanStatus::Certified
                    && plan.rejection == EvalRoundRejectionReason::CleaningDomainViolation,
                  "all local budgets stop at the unchanged cleaner-domain gate");
        }
        auto fixture = test::BootstrapFixture::create(
            {16, std::vector<int>(10, 50), std::ldexp(1., 49), 8});
        fixture.generateKeys(std::vector<int>{0}, true);
        const auto input = fixture.encrypt(fixture.encode({0.}));
        EvalRoundExecutionOptions options;
        options.inputSemanticError = {1e-8, BootstrapBoundKind::Deterministic,
                                      "analysis-only compiler input bound", {}};
        const auto rejectedMath = planEvalRoundCandidate(problem, candidate);
        const auto compileAttempt = EvalRoundExecutionCompiler::compile(
            fixture, input, rejectedMath, options);
        check(compileAttempt.certification().status
                == BootstrapCertificationStatus::ExtractionNotCertified,
              "analysis compile attempt stops at the mathematical cleaner-domain blocker");
        auto missing = candidate; missing.extraction.polynomials.pop_back();
        check(missing.extraction.polynomials.size() != 8,
              "removing one digit cannot leave a complete extractor");
        auto wrongIndex = candidate.extraction.polynomials[0]; wrongIndex.digitIndex = 1;
        auto wrongK = candidate.extraction.polynomials[0]; wrongK.certifiedK = 1;
        auto wrongRho = candidate.extraction.polynomials[0]; wrongRho.certifiedRho = std::nextafter(problem.rho, 1.0);
        auto modified = candidate.extraction.polynomials[0]; modified.polynomial.decimalCoefficients[0] += "1";
        check(wrongIndex.digitIndex != 0 && wrongK.certifiedK != problem.K
                && wrongRho.certifiedRho != problem.rho
                && modified.polynomial.decimalCoefficients != candidate.extraction.polynomials[0].polynomial.decimalCoefficients,
              "certificate-binding mutations are observable and require rejection/recomputation");
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
