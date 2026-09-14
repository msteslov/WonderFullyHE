#include "m2424/experimental/evalmod_analysis/evalround_synthesis.hpp"
#include "bootstrap_fixture.hpp"
#include "../src/planning/certified_arithmetic_internal.hpp"
#include "../src/planning/evalround_polynomial_internal.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <limits>
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
            check(record.directXRigorousIntervalError==proof.directXApproximationError.upperBound
                  &&record.centeredRigorousIntervalError==proof.centeredApproximationError.upperBound
                  &&record.rigorousIntervalError==std::min(record.directXRigorousIntervalError,
                                                           record.centeredRigorousIntervalError)
                  &&record.selectedIntervalProofMethod==proof.selectedIntervalProofMethod,
                  "record keeps both independent bounds and selects only their rigorous minimum");
        }
    }
    for (const auto& digits : remezAttempted) for (const bool attempted : digits)
        check(attempted, "all eight digits and four degrees are attempted by Remez");
    check(sawConvergedRemez,
          "bounded search retains Remez convergence as diagnostic metadata");
    for (std::size_t digit = 0; digit < 8; ++digit) {
        const auto& selected = search.records[*search.selectedRecordByDigit[digit]];
        const std::array<EvalModApproximationFamily,8> expectedFamilies{
            EvalModApproximationFamily::MultiIntervalChebyshev,
            EvalModApproximationFamily::MultiIntervalChebyshev,
            EvalModApproximationFamily::MultiIntervalChebyshev,
            EvalModApproximationFamily::MultiIntervalChebyshev,
            EvalModApproximationFamily::MultiIntervalChebyshev,
            EvalModApproximationFamily::MultiIntervalChebyshev,
            EvalModApproximationFamily::MultiIntervalChebyshev,
            EvalModApproximationFamily::MultiIntervalMinimax};
        const std::array<std::size_t,8> expectedDegrees{256,128,256,256,256,256,192,256};
        check(selected.family==expectedFamilies[digit]
              &&selected.requestedDegree==expectedDegrees[digit],
              "selected K64 families and degrees remain unchanged");
        check(selected.certificate&&selected.certificate->executionRepresentation,
              "every selected K64 polynomial retains an execution representation");
        const auto& execution=*selected.certificate->executionRepresentation;
        check(execution.polynomial.basis==PolynomialBasis::Chebyshev
              &&execution.variableScaleDecimal=="64"
              &&!execution.exactEquivalenceProvenance.empty(),
              "selected execution metadata is scaled Chebyshev with exact-equivalence provenance");
        check(arithmetic::exactPolynomialEqual(
                  convertScaledChebyshevToMonomial(execution.polynomial,"64"),
                  selected.certificate->polynomial),
              "all eight selected execution polynomials are exactly equivalent to canonical p_j");
        for (const auto& record : search.records)
            if (record.digitIndex == digit && record.certificate
                && std::isfinite(record.rigorousIntervalError))
                check(selected.rigorousIntervalError <= record.rigorousIntervalError,
                      "selection uses the rigorous bound across Remez and diagnostics");
    }
    auto candidate = makeEvalRoundBinaryPolynomialCandidate(search);
    if (search.allDigitsCertified) {
        check(search.allCleanerInputDomainsSatisfied && search.extractorCertified,
              "centered whole-domain bounds put all selected digits inside a<=1");
        check(candidate.extraction.verified && candidate.extraction.polynomials.size()==8,
              "cleaner-feasible search constructs the complete binary candidate");
        std::optional<EvalRoundPlan> shallowestPlan;
        const std::vector<std::vector<std::size_t>> expectedRounds{
            {0,8,3,4,4,3,2,0},{0,9,4,5,5,4,3,0},{0,9,4,5,5,5,4,1}};
        std::size_t budgetIndex=0;
        for (const double budget : {1e-2, 1e-4, 1e-6}) {
            auto diagnostic=problem;diagnostic.requiredIntegerError=budget;
            const auto plan=planEvalRoundCandidate(diagnostic,candidate);
            check(plan.status==EvalRoundPlanStatus::Certified,
                  "unchanged cleaner/reconstruction planner reaches every diagnostic budget");
            check(plan.integerErrorUpper<=budget&&plan.digits.size()==8,
                  "planner reports a complete bounded reconstruction schedule");
            for(std::size_t digit=0;digit<8;++digit)
                check(plan.digits[digit].cleaningIterations==expectedRounds[budgetIndex][digit],
                      "unchanged planner retains the minimum certified cleaning schedule");
            if(!shallowestPlan)shallowestPlan=plan;
            ++budgetIndex;
        }
        auto adapter=test::BootstrapFixture::create(
            {16,std::vector<int>(48,50),std::ldexp(1.,49),8});
        adapter.generateKeys(std::vector<int>{0},true);
        const auto input=adapter.encrypt(adapter.encode({0.}));
        EvalRoundExecutionOptions options;
        options.inputSemanticError={1e-12,BootstrapBoundKind::Deterministic,
            "analysis-only backend feasibility input bound",{}};
        std::vector<std::vector<std::string>> exactCoefficients;
        std::vector<double> centeredCertificates;
        for(const auto& polynomial:candidate.extraction.polynomials) {
            exactCoefficients.push_back(polynomial.polynomial.decimalCoefficients);
            centeredCertificates.push_back(polynomial.centeredApproximationError.upperBound);
        }
        const mpq_class huge(mpz_class(1)<<1100);
        arithmetic::Builder exactBuilder(adapter,input,options.evaluationKeyNoiseCoefficientSupport);
        const auto hugeNode=exactBuilder.input(1,huge,0);
        const auto& hugeTrace=exactBuilder.nodes[hugeNode];
        check(hugeTrace.exactIdealMagnitude.kind==BootstrapBoundKind::Deterministic
              &&!hugeTrace.exactIdealMagnitude.outwardBinary64
              &&hugeTrace.idealMagnitude.kind==BootstrapBoundKind::Deterministic
              &&std::isinf(hugeTrace.idealMagnitude.upperBound)
              &&mpq_class(mpz_class(hugeTrace.exactIdealMagnitude.numerator),
                          mpz_class(hugeTrace.exactIdealMagnitude.denominator))==huge
              &&mpz_class(hugeTrace.centeredHeadroomNumerator)>0,
              "finite exact bound above DBL_MAX stays known, unclamped, nonzero, and can pass exact headroom");
        bool exactHeadroomRejected=false;
        try {
            arithmetic::Builder insufficient(adapter,input,options.evaluationKeyNoiseCoefficientSupport);
            insufficient.input(1,mpq_class(mpz_class(1)<<2500),0);
        } catch(const arithmetic::Failure& failure) {
            exactHeadroomRejected=failure.status==BootstrapCertificationStatus::HeadroomViolation;
        }
        check(exactHeadroomRejected,
              "finite exact bound above DBL_MAX fails through exact HeadroomViolation when Q is insufficient");
        arithmetic::Builder hand(adapter,input,options.evaluationKeyNoiseCoefficientSupport);
        const auto handInput=hand.input(2,mpq_class(2),mpq_class(1,10));
        const auto handProduct=hand.mul(handInput,handInput,"hand exact propagation");
        check(hand.states[handProduct].M==4&&hand.states[handProduct].E==mpq_class(41,100),
              "exact product propagation matches M=4 and E=41/100 without binary64 recurrence");
        const auto compiled=EvalRoundExecutionCompiler::compile(
            adapter,input,*shallowestPlan,options);
        std::cout<<std::setprecision(17);
        std::cout<<"backend_status="<<static_cast<int>(compiled.certification().status)
                 <<" nodes="<<compiled.nodes().size()
                 <<" blocker="<<compiled.certification().provenance<<'\n';
        check(compiled.certification().status==BootstrapCertificationStatus::HeadroomViolation
              &&compiled.certification().provenance.find(
                  "No certified backend cleaning schedule")!=std::string::npos
              &&compiled.certification().provenance.find(
                  "first concrete backend limit")!=std::string::npos
              &&compiled.certification().provenance.find(
                  "Exact polynomial common denominator has no finite binary64 plaintext scale")==std::string::npos
              &&compiled.certification().provenance.find(
                  "Extraction digit 0 exact arithmetic error exceeds cleaner domain")==std::string::npos,
              "scaled-Chebyshev execution passes extraction and reaches the later concrete headroom gate");
        check(compiled.nodes().empty(),"non-Certified K64 compilation never publishes an executable DAG");
        const auto& diagnostic=compiled.diagnostics();
        check(diagnostic.digits.size()==8&&diagnostic.constructedNodes>0
              &&diagnostic.ciphertextMultiplications>0
              &&diagnostic.ciphertextMultiplications==diagnostic.relinearizations
              &&diagnostic.rescales>=diagnostic.ciphertextMultiplications
              &&diagnostic.criticalMultiplicativeDepth<256
              &&diagnostic.minimumRuntimeScale&&diagnostic.maximumRuntimeScale
              &&!diagnostic.minimumCenteredHeadroomNumerator.empty()
              &&mpz_class(diagnostic.minimumCenteredHeadroomNumerator)>0,
              "fail-closed K64 result retains exact resource diagnostics up to its first gate");
        for(std::size_t digit=0;digit<diagnostic.digits.size();++digit) {
            const auto& trace=diagnostic.digits[digit];
            check(trace.approximationError.upperBound
                      ==candidate.digits[digit].extractionError.upperBound
                  &&trace.backendExtractionError.outwardBinary64
                  &&trace.initialCleanerError.outwardBinary64
                  &&*trace.initialCleanerError.outwardBinary64<=1
                  &&!trace.cleanerLocalErrors.empty()
                  &&trace.cleanerErrorAfterRounds.size()
                      ==trace.cleanerLocalErrors.size()+1
                  &&trace.chebyshevNodeCount>0,
                  "every K64 digit passes exact a_0<=1 and reaches real cleaner arithmetic");
            std::cout<<"digit="<<digit
                     <<" E_approx="<<trace.approximationError.upperBound
                     <<" E_backend="<<*trace.backendExtractionError.outwardBinary64
                     <<" a0="<<*trace.initialCleanerError.outwardBinary64
                     <<" cleaner_rounds_reached="<<trace.cleanerLocalErrors.size()
                     <<" chebyshev_nodes="<<trace.chebyshevNodeCount<<'\n';
        }
        std::cout<<"constructed_nodes="<<diagnostic.constructedNodes
                 <<" ctct="<<diagnostic.ciphertextMultiplications
                 <<" relin="<<diagnostic.relinearizations
                 <<" rescale="<<diagnostic.rescales
                 <<" modswitch="<<diagnostic.modSwitches
                 <<" mulplain="<<diagnostic.plaintextMultiplications
                 <<" depth="<<diagnostic.criticalMultiplicativeDepth
                 <<" levels="<<diagnostic.criticalPathLevelConsumption
                 <<" scale_min="<<*diagnostic.minimumRuntimeScale
                 <<" scale_max="<<*diagnostic.maximumRuntimeScale<<'\n';
        for(std::size_t digit=0;digit<candidate.extraction.polynomials.size();++digit) {
            check(candidate.extraction.polynomials[digit].polynomial.decimalCoefficients
                      ==exactCoefficients[digit]
                  &&candidate.extraction.polynomials[digit].centeredApproximationError.upperBound
                      ==centeredCertificates[digit],
                  "K64 exact polynomial strings and centered certificates remain unchanged");
        }
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
