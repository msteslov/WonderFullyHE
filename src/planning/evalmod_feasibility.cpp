#include "m2424/experimental/evalmod_analysis/feasibility.hpp"

#include "m2424/security_report.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace m2424::experimental {
namespace {

double scaleValue(const ExactScale& scale) {
    return mpz_get_d(scale.numerator.get_mpz_t())
        / mpz_get_d(scale.denominator.get_mpz_t());
}

std::size_t bitCount(std::uint64_t value) {
    std::size_t bits = 0;
    while (value) { ++bits; value >>= 1; }
    return bits;
}

double divideRoundFloor(std::size_t degree, double scale) {
    const long double n = static_cast<long double>(degree);
    return static_cast<double>((0.5L * n * (1.0L + n)) / scale);
}

double quadraticApproximationLowerBound(double rho) {
    // At -1, 0, 1 the target is zero. Lagrange interpolation bounds any
    // quadratic at rho by eps*(1+rho-rho^2). Comparing with target(rho)=rho
    // gives eps >= rho/(2+rho-rho^2).
    return rho / (2.0 + rho - rho * rho);
}

std::size_t countOperation(const CompiledEvalModCircuit& circuit, EvalModOperation operation) {
    return static_cast<std::size_t>(std::count_if(
        circuit.nodes.begin(), circuit.nodes.end(), [&](const auto& node) {
            return node.operation == operation;
        }));
}

std::size_t maximumCompilerDegreeUpperBound(
    const EvalModProblem& source, std::size_t levels) {
    auto problem = source;
    problem.exactModulusContext.reset();
    problem.requireCertifiedScaleSchedule = false;
    std::size_t maximum = 0;
    // Odd dense polynomials are the relevant direct EvalMod class. The legacy
    // scale aligner can only consume no more levels than the exact compiler,
    // so this enumeration is a resource upper bound, not an optimistic plan.
    // If babyStep > L+1, baby powers alone exceed L. Otherwise a degree above
    // (L+1)^2 has more than L+1 giant blocks and also exceeds L. Thus this is
    // a finite exhaustive bound for this compiler, not a search heuristic.
    const std::size_t exhaustiveDegreeLimit = (levels + 1) * (levels + 1);
    for (std::size_t degree = 1; degree <= exhaustiveDegreeLimit; degree += 2) {
        bool feasible = false;
        EvalModPolynomial polynomial;
        polynomial.basis = PolynomialBasis::Monomial;
        polynomial.decimalCoefficients.assign(degree + 1, "0");
        for (std::size_t exponent = 1; exponent <= degree; exponent += 2)
            polynomial.decimalCoefficients[exponent] = "1";
        for (std::size_t babyStep = 2;
             babyStep <= std::min<std::size_t>(degree + 1, levels + 1); ++babyStep) {
            try {
                const auto circuit = compileEvalModPolynomial(polynomial, problem, babyStep);
                feasible = feasible || circuit.cost.levelConsumption <= levels;
            } catch (const std::exception&) {
            }
        }
        if (feasible) maximum = degree;
    }
    return maximum;
}

double optimisticFirstUseCbdKeySwitch(
    const PreparedEvalModPlan& plan, std::size_t relinearizations) {
    if (!relinearizations || plan.dataPrimes.empty() || plan.specialPrime < 2
        || !plan.constants.polyModulusDegree) return 0.0;
    const long double n = static_cast<long double>(plan.constants.polyModulusDegree);
    long double squaredWeights = 0.0L;
    for (const auto prime : plan.dataPrimes) {
        const long double weight = n * static_cast<long double>(prime - 1);
        squaredWeights += weight * weight;
    }
    // CBD(eta=21) is a sum of 42 centered Bernoulli variables and has the
    // exact Hoeffding MGF bound E exp(tX) <= exp(10.5 t^2 / 2).
    // We additionally (and optimistically) grant non-adaptive first-use
    // decomposition weights. Union bound covers real/imaginary coordinates,
    // all canonical roots, both ModDown outputs, and every relinearization.
    const long double events = 4.0L * n * 2.0L
        * static_cast<long double>(relinearizations);
    const long double logTail = std::log(events) + 128.0L * std::log(2.0L);
    const long double varianceProxy = 10.5L * n * squaredWeights;
    const long double embedding = std::sqrt(4.0L * varianceProxy * logTail);
    double minimumRelinearizationScale = std::numeric_limits<double>::infinity();
    for (std::size_t index = 0; index < plan.circuit.nodes.size(); ++index)
        if (plan.circuit.nodes[index].operation == EvalModOperation::Relinearize)
            minimumRelinearizationScale = std::min(
                minimumRelinearizationScale,
                scaleValue(plan.arithmeticCertificate.nodes[index].scale.runtimeExact));
    const long double perRelinearization = embedding
        / static_cast<long double>(plan.specialPrime)
        / static_cast<long double>(minimumRelinearizationScale);
    return static_cast<double>(
        perRelinearization * static_cast<long double>(relinearizations));
}

std::string number(double value) {
    if (!std::isfinite(value)) return "null";
    std::ostringstream out;
    out << std::setprecision(17) << value;
    return out.str();
}

std::string escapeJson(const std::string& input) {
    std::string output;
    output.reserve(input.size());
    for (const char value : input) {
        if (value == '\\' || value == '"') output.push_back('\\');
        if (value == '\n') output += "\\n";
        else output.push_back(value);
    }
    return output;
}

} // namespace

EvalModFeasibilityReport analyzeEvalModFeasibility(
    const EvalModProblem& problem,
    const PreparedEvalModPlan& plan,
    const EvalModSynthesisResult* diagnosticSynthesis) {
    if (plan.dataPrimes.empty() || !plan.constants.polyModulusDegree
        || !plan.arithmeticCertificate.outputError.rigorous
        || plan.domain.integerBound == 0 || !(problem.targetAbsoluteError > 0.0))
        throw std::invalid_argument("representative nonlinear prepared plan is incomplete");
    if (plan.circuit.normalizationGain.numerator * 2
            != plan.circuit.normalizationGain.denominator
        || plan.circuit.nodes.size() < 4
        || plan.circuit.nodes[0].operation != EvalModOperation::Input
        || plan.circuit.nodes[1].operation != EvalModOperation::EncodeConstant
        || plan.circuit.nodes[2].operation != EvalModOperation::MultiplyPlain
        || plan.circuit.nodes[3].operation != EvalModOperation::Rescale)
        throw std::invalid_argument(
            "feasibility proof requires the explicit alpha=1/2 normalization prefix");

    EvalModFeasibilityReport report;
    report.polyModulusDegree = plan.constants.polyModulusDegree;
    report.integerBound = plan.domain.integerBound;
    report.rho = plan.domain.normalizedResidualBound;
    report.currentSecurityBits = plan.securityBits;
    report.tc128CoeffModulusBits = static_cast<std::size_t>(
        coeff_modulus_max_bit_count(report.polyModulusDegree, SecurityLevel::TC128));
    report.currentCoeffModulusBits = bitCount(plan.specialPrime);
    for (const auto prime : plan.dataPrimes)
        report.currentCoeffModulusBits += bitCount(prime);
    const std::size_t workingBits = static_cast<std::size_t>(
        std::ceil(std::log2(scaleValue(problem.coeffToSlotScale))));
    const std::size_t specialBits = bitCount(plan.specialPrime);
    report.maximumSecurityValidDataPrimes = report.tc128CoeffModulusBits > specialBits
        ? (report.tc128CoeffModulusBits - specialBits) / workingBits : 0;
    report.maximumSecurityValidLevels = report.maximumSecurityValidDataPrimes
        ? report.maximumSecurityValidDataPrimes - 1 : 0;
    report.maximumResourceFeasibleDegreeUpperBound = maximumCompilerDegreeUpperBound(
        problem, report.maximumSecurityValidLevels);
    report.outputGain = scaleValue(plan.denormalizationGain);
    report.normalizedTarget = plan.normalizedEvalModBudget;
    report.finalTarget = problem.targetAbsoluteError;

    const double normalizationRescaleFloor =
        plan.arithmeticCertificate.nodes[3].localAddedError.upperBound;
    report.mandatoryCurrentRescaleFloor = normalizationRescaleFloor;
    report.mandatoryRescaleImprovementFactor =
        normalizationRescaleFloor / report.normalizedTarget;
    const std::size_t relins = countOperation(plan.circuit, EvalModOperation::Relinearize);
    report.optimisticFirstUseCbdKeySwitchBound =
        optimisticFirstUseCbdKeySwitch(plan, relins);
    double minimumRelinearizationScale = std::numeric_limits<double>::infinity();
    for (std::size_t index = 0; index < plan.circuit.nodes.size(); ++index)
        if (plan.circuit.nodes[index].operation == EvalModOperation::Relinearize)
            minimumRelinearizationScale = std::min(
                minimumRelinearizationScale,
                scaleValue(plan.arithmeticCertificate.nodes[index].scale.runtimeExact));
    report.optimisticKeySwitchPlusModDown =
        report.optimisticFirstUseCbdKeySwitchBound
        + relins * divideRoundFloor(report.polyModulusDegree,
                                    minimumRelinearizationScale);

    const auto& breakdown = plan.arithmeticCertificate.outputBreakdown;
    const double finalFactor = report.outputGain * problem.slotToCoeffOperatorNorm;
    const auto addSource = [&](std::string source, double current, double best,
                               double failure, std::size_t occurrences, bool rigorous,
                               std::string note) {
        report.sources.push_back({std::move(source), current, best, failure, occurrences,
                                  current * finalFactor, rigorous, std::move(note)});
    };
    addSource("Input/CtS error", breakdown.inputCoeffToSlot,
              breakdown.inputCoeffToSlot, -std::numeric_limits<double>::infinity(), 1, true,
              "exact DAG amplification; upstream error is also reserved by the bootstrap gate");
    const double currentApproximation = plan.approximationError.rigorous
        ? plan.approximationError.upperBound : std::numeric_limits<double>::infinity();
    addSource("Approximation error", currentApproximation,
              currentApproximation, -std::numeric_limits<double>::infinity(), 1,
              plan.approximationError.rigorous, "MPFR outward-rounded interval certificate");
    addSource("Coefficient quantization", breakdown.coefficientQuantization,
              breakdown.coefficientQuantization,
              -std::numeric_limits<double>::infinity(), plan.constants.constants.size(), true,
              "prepared exact-RNS constants with DAG amplification");
    addSource("MultiplyPlain propagation", 0.0, 0.0,
              -std::numeric_limits<double>::infinity(),
              countOperation(plan.circuit, EvalModOperation::MultiplyPlain), true,
              "adds no local error; amplification is retained on each originating source");
    addSource("CipherxCipher propagation", breakdown.multiplicationInteraction,
              breakdown.multiplicationInteraction,
              -std::numeric_limits<double>::infinity(),
              countOperation(plan.circuit, EvalModOperation::MultiplyCipher), true,
              "all E_a*E_b cross terms after exact value-bound amplification");
    addSource("Rescale", breakdown.rescale, normalizationRescaleFloor,
              -std::numeric_limits<double>::infinity(),
              countOperation(plan.circuit, EvalModOperation::Rescale), true,
              "SEAL divide-and-round coefficient support and canonical embedding; best column is the mandatory single-operation floor before sensitivity propagation");
    addSource("Relinearization/KeySwitch", breakdown.keySwitch, breakdown.keySwitch,
              -std::numeric_limits<double>::infinity(),
              countOperation(plan.circuit, EvalModOperation::Relinearize), true,
              "CBD support, exact decomposition primes, special-prime ModDown");
    addSource("ModSwitch", 0.0, 0.0, -std::numeric_limits<double>::infinity(),
              countOperation(plan.circuit, EvalModOperation::ModSwitch), true,
              "zero semantic error after the certified no-wrap headroom gate");
    addSource("Scale representation error", 0.0, 0.0,
              -std::numeric_limits<double>::infinity(), 0, true,
              "actual binary64 scale schedule; metadata-only alignment is prohibited");
    addSource("Period mismatch", plan.domain.errors.periodMismatch,
              plan.domain.errors.periodMismatch,
              problem.ciphertextModel.tailModel == TailModel::Deterministic
                  ? -std::numeric_limits<double>::infinity()
                  : plan.domain.failureProbabilityLog2,
              1, true,
              "outward-rounded domain model");
    addSource("Denormalization", 0.0, 0.0,
              -std::numeric_limits<double>::infinity(), 1, true,
              "adds no local error; all normalized rows are multiplied by the exact output gain");
    report.sources.push_back({"SlotToCoeff contribution", problem.slotToCoeffAdditive,
                              problem.slotToCoeffAdditive,
                              -std::numeric_limits<double>::infinity(), 1,
                              problem.slotToCoeffAdditive, true,
                              "operator norm is included in every propagated row"});
    report.sources.push_back({"Final additive error", problem.finalAdditive,
                              problem.finalAdditive,
                              -std::numeric_limits<double>::infinity(), 1,
                              problem.finalAdditive, true, "final bootstrap additive term"});

    const std::vector<std::pair<std::string, double>> dominantCandidates{
        {"Input/CtS error", breakdown.inputCoeffToSlot},
        {"Coefficient quantization", breakdown.coefficientQuantization},
        {"CipherxCipher propagation", breakdown.multiplicationInteraction},
        {"Rescale", breakdown.rescale}, {"Relinearization/KeySwitch", breakdown.keySwitch}};
    const auto dominant = std::max_element(
        dominantCandidates.begin(), dominantCandidates.end(),
        [](const auto& left, const auto& right) { return left.second < right.second; });
    report.dominantSource = dominant->first;
    report.dominantCurrentBound = dominant->second;
    report.requiredImprovementFactor = report.dominantCurrentBound / report.normalizedTarget;
    report.compositeBranchMargin = 0.5 - report.rho;
    report.quadraticApproximationLowerBound = quadraticApproximationLowerBound(report.rho);
    report.centralSensitivityLowerBound = std::max(
        0.0, 1.0 - report.normalizedTarget / report.rho);

    const auto addFamily = [&](EvalModApproximationFamily family, const char* name) {
        EvalModFamilyFeasibility item;
        item.family = name;
        item.resourceDegreeUpperBound = report.maximumResourceFeasibleDegreeUpperBound;
        item.mandatoryArithmeticFloor = report.mandatoryCurrentRescaleFloor
            * report.centralSensitivityLowerBound;
        item.status = "PrunedByMandatoryNormalizationRescale";
        if (diagnosticSynthesis) {
            item.diagnosticMinimumDegree = std::numeric_limits<std::size_t>::max();
            for (const auto& candidate : diagnosticSynthesis->candidates) {
                if (candidate.family != family) continue;
                const auto degree = candidate.compiledCircuit.cost.degree;
                item.diagnosticMinimumDegree = std::min(item.diagnosticMinimumDegree, degree);
                item.diagnosticMaximumDegree = std::max(item.diagnosticMaximumDegree, degree);
                item.maximumDepth = std::max(
                    item.maximumDepth, candidate.compiledCircuit.cost.multiplicativeDepth);
                if (candidate.intervalCertificate.proved)
                    item.bestDiagnosticApproximation = std::min(
                        item.bestDiagnosticApproximation,
                        candidate.intervalCertificate.approximationErrorUpperBoundDouble);
            }
            if (item.diagnosticMinimumDegree == std::numeric_limits<std::size_t>::max())
                item.diagnosticMinimumDegree = 0;
        }
        report.families.push_back(std::move(item));
    };
    addFamily(EvalModApproximationFamily::MultiIntervalMinimax, "DirectMinimax");
    addFamily(EvalModApproximationFamily::MultiIntervalChebyshev, "DirectChebyshev");
    EvalModFamilyFeasibility composite;
    composite.family = "CompositeArcsinSine";
    composite.resourceDegreeUpperBound = report.maximumResourceFeasibleDegreeUpperBound;
    composite.mandatoryArithmeticFloor = report.mandatoryCurrentRescaleFloor
        * report.centralSensitivityLowerBound;
    composite.status = "PrunedByFeasibilityBeforeImplementation";
    report.families.push_back(std::move(composite));

    for (const std::size_t degree : {1024U, 2048U, 4096U, 8192U, 16384U, 32768U}) {
        EvalModProfileFeasibility profile;
        profile.polyModulusDegree = degree;
        profile.tc128CoeffModulusBits = static_cast<std::size_t>(
            coeff_modulus_max_bit_count(degree, SecurityLevel::TC128));
        profile.stableScaleBits = std::min<std::size_t>(
            60, profile.tc128CoeffModulusBits / 3);
        profile.maximumDataPrimes = profile.stableScaleBits
            ? (profile.tc128CoeffModulusBits - profile.stableScaleBits)
                / profile.stableScaleBits : 0;
        profile.maximumLevels = profile.maximumDataPrimes
            ? profile.maximumDataPrimes - 1 : 0;
        profile.mandatoryRescaleFloor = profile.stableScaleBits
            ? divideRoundFloor(degree, std::ldexp(1.0, profile.stableScaleBits))
            : std::numeric_limits<double>::infinity();
        // Within scope S_out <= qSource, so the normalized budget cannot exceed
        // the final target (output gain >= 1).
        profile.maximumNormalizedBudget = problem.targetAbsoluteError;
        profile.excludedByArithmetic = profile.mandatoryRescaleFloor
                * std::max(0.0, 1.0 - profile.maximumNormalizedBudget / report.rho)
            > profile.maximumNormalizedBudget;
        profile.approximationLowerBound = profile.maximumLevels <= 1
            ? report.quadraticApproximationLowerBound : 0.0;
        profile.excludedByApproximation = !profile.excludedByArithmetic
            && profile.maximumLevels <= 1
            && profile.approximationLowerBound > profile.maximumNormalizedBudget;
        if (profile.excludedByArithmetic)
            profile.stopCondition = "the explicit alpha-normalization rescale floor times the central-interval sensitivity lower bound exceeds the maximum normalized budget";
        else if (profile.excludedByApproximation)
            profile.stopCondition = "one available level limits the nonlinear polynomial to degree <=2; interpolation lower bound exceeds target";
        else
            profile.stopCondition = "not excluded by this certificate-class proof";
        report.profiles.push_back(std::move(profile));
    }

    const bool profilesExhausted = std::all_of(
        report.profiles.begin(), report.profiles.end(), [](const auto& profile) {
            return profile.excludedByArithmetic || profile.excludedByApproximation;
        });
    report.status = profilesExhausted
        ? EvalModFeasibilityStatus::CertificateClassInfeasible
        : EvalModFeasibilityStatus::CandidateFeasible;
    report.probabilisticStopCondition =
        "KeySwitch is not the dominant source at its actual pre-rescale scale. A CBD tail can only improve that non-bottleneck, while the required Rescale improvement is at least the reported mandatory-floor factor and 1.7e5 for the representative DAG. Rescale residual independence/zero mean is not proved, so concentration is invalid.";
    report.directStopCondition =
        "The current N=32768 compiler has the reported resource-only dense odd degree upper bound, but every degree is pruned first by the explicit alpha-normalization rescale floor and central-interval sensitivity 1-epsilon/rho. Across profiles, N>=16384 has the same arithmetic stop; the only arithmetically non-excluded N=8192 profile has one level and degree<=2, whose rigorous multi-interval approximation lower bound is rho/(2+rho-rho^2).";
    report.compositeStopCondition =
        "The arcsin(sin(2*pi*z))/(2*pi) branch margin is positive, but under the same explicit alpha-normalization compiler every composition enters the normalization-rescale/sensitivity stop at N>=16384 or requires degree>2 after the one available level at N=8192. Feasibility therefore rejects implementation.";
    report.proofScope =
        "All plans emitted by the current explicit-alpha stable-scale polynomial EvalMod compiler (direct or composite) on SEAL standard tc128 degrees N=1024..32768, coefficient primes <=60 bits, a same-scale special prime, S_out<=qSource, and the implemented deterministic coefficient-support/canonical-embedding arithmetic certificate. The proof uses the mandatory alpha=1/2 normalization rescale and the sensitivity any target-accurate polynomial must have on [-rho,rho]. This is not a global impossibility result for EvalMod.";
    report.remainingAlternatives =
        "A compiler that folds alpha into polynomial coefficients and proves a different error layout; a backend with a proved stochastic rounding law independent of adaptive ciphertexts; a non-stationary/no-rescale exact scale layout with a new headroom proof; or a non-polynomial bootstrapping construction/backend.";
    return report;
}

std::string evalModFeasibilityJson(const EvalModFeasibilityReport& report) {
    std::ostringstream out;
    out << std::setprecision(17)
        << "{\n  \"status\":\""
        << (report.status == EvalModFeasibilityStatus::CertificateClassInfeasible
                ? "CertificateClassInfeasible"
                : report.status == EvalModFeasibilityStatus::CandidateFeasible
                    ? "CandidateFeasible" : "RigorousBoundUnavailable")
        << "\",\n  \"global_impossibility_proved\":"
        << (report.globalImpossibilityProved ? "true" : "false")
        << ",\n  \"N\":" << report.polyModulusDegree
        << ",\n  \"K\":" << report.integerBound
        << ",\n  \"rho\":" << number(report.rho)
        << ",\n  \"security_bits\":" << report.currentSecurityBits
        << ",\n  \"tc128_coeff_modulus_bits\":" << report.tc128CoeffModulusBits
        << ",\n  \"current_coeff_modulus_bits\":" << report.currentCoeffModulusBits
        << ",\n  \"max_security_valid_data_primes\":"
        << report.maximumSecurityValidDataPrimes
        << ",\n  \"max_security_valid_levels\":"
        << report.maximumSecurityValidLevels
        << ",\n  \"max_resource_feasible_degree_upper_bound\":"
        << report.maximumResourceFeasibleDegreeUpperBound
        << ",\n  \"normalized_target\":" << number(report.normalizedTarget)
        << ",\n  \"final_target\":" << number(report.finalTarget)
        << ",\n  \"output_gain\":" << number(report.outputGain)
        << ",\n  \"mandatory_current_rescale_floor\":"
        << number(report.mandatoryCurrentRescaleFloor)
        << ",\n  \"mandatory_rescale_improvement_factor\":"
        << number(report.mandatoryRescaleImprovementFactor)
        << ",\n  \"optimistic_first_use_cbd_keyswitch\":"
        << number(report.optimisticFirstUseCbdKeySwitchBound)
        << ",\n  \"optimistic_keyswitch_plus_moddown\":"
        << number(report.optimisticKeySwitchPlusModDown)
        << ",\n  \"dominant_source\":\"" << escapeJson(report.dominantSource)
        << "\",\n  \"dominant_current_bound\":" << number(report.dominantCurrentBound)
        << ",\n  \"required_improvement_factor\":"
        << number(report.requiredImprovementFactor)
        << ",\n  \"composite_branch_margin\":" << number(report.compositeBranchMargin)
        << ",\n  \"quadratic_approximation_lower_bound\":"
        << number(report.quadraticApproximationLowerBound)
        << ",\n  \"central_sensitivity_lower_bound\":"
        << number(report.centralSensitivityLowerBound)
        << ",\n  \"sources\":[";
    for (std::size_t index = 0; index < report.sources.size(); ++index) {
        const auto& source = report.sources[index];
        if (index) out << ',';
        out << "\n    {\"source\":\"" << escapeJson(source.source)
            << "\",\"current_normalized\":" << number(source.currentNormalizedBound)
            << ",\"best_plausible_normalized\":"
            << number(source.bestPlausibleNormalizedBound)
            << ",\"log2_failure\":" << number(source.log2FailureProbability)
            << ",\"occurrences\":" << source.occurrences
            << ",\"propagated_final\":" << number(source.propagatedFinalContribution)
            << ",\"rigorous\":" << (source.rigorous ? "true" : "false")
            << ",\"note\":\"" << escapeJson(source.note) << "\"}";
    }
    out << "\n  ],\n  \"profiles\":[";
    for (std::size_t index = 0; index < report.profiles.size(); ++index) {
        const auto& profile = report.profiles[index];
        if (index) out << ',';
        out << "\n    {\"N\":" << profile.polyModulusDegree
            << ",\"tc128_bits\":" << profile.tc128CoeffModulusBits
            << ",\"stable_scale_bits\":" << profile.stableScaleBits
            << ",\"max_data_primes\":" << profile.maximumDataPrimes
            << ",\"max_levels\":" << profile.maximumLevels
            << ",\"rescale_floor\":" << number(profile.mandatoryRescaleFloor)
            << ",\"max_normalized_budget\":" << number(profile.maximumNormalizedBudget)
            << ",\"approximation_lower_bound\":"
            << number(profile.approximationLowerBound)
            << ",\"excluded_by_arithmetic\":"
            << (profile.excludedByArithmetic ? "true" : "false")
            << ",\"excluded_by_approximation\":"
            << (profile.excludedByApproximation ? "true" : "false")
            << ",\"stop_condition\":\"" << escapeJson(profile.stopCondition) << "\"}";
    }
    out << "\n  ],\n  \"families\":[";
    for (std::size_t index = 0; index < report.families.size(); ++index) {
        const auto& family = report.families[index];
        if (index) out << ',';
        out << "\n    {\"family\":\"" << escapeJson(family.family)
            << "\",\"diagnostic_min_degree\":" << family.diagnosticMinimumDegree
            << ",\"diagnostic_max_degree\":" << family.diagnosticMaximumDegree
            << ",\"resource_degree_upper_bound\":"
            << family.resourceDegreeUpperBound
            << ",\"maximum_diagnostic_depth\":" << family.maximumDepth
            << ",\"best_diagnostic_approximation\":"
            << number(family.bestDiagnosticApproximation)
            << ",\"mandatory_arithmetic_floor\":"
            << number(family.mandatoryArithmeticFloor)
            << ",\"status\":\"" << escapeJson(family.status) << "\"}";
    }
    out << "\n  ],\n  \"probabilistic_stop\":\""
        << escapeJson(report.probabilisticStopCondition)
        << "\",\n  \"direct_stop\":\"" << escapeJson(report.directStopCondition)
        << "\",\n  \"composite_stop\":\"" << escapeJson(report.compositeStopCondition)
        << "\",\n  \"proof_scope\":\"" << escapeJson(report.proofScope)
        << "\",\n  \"remaining_alternatives\":\""
        << escapeJson(report.remainingAlternatives) << "\"\n}\n";
    return out.str();
}

} // namespace m2424::experimental
