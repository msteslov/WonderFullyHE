#include "m2424/experimental/evalmod_analysis/evalround_synthesis.hpp"

#include "m2424/experimental/evalmod_analysis/exact_decimal.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace m2424::experimental {
namespace {

long double digitAtInteger(std::int64_t integer, std::size_t digit) {
    // The generator's continuation is explicitly periodic. Certification does
    // not use it: the proof target remains bit_j(I+K) on the requested cells.
    const std::int64_t period = std::int64_t(1) << (digit + 1);
    std::int64_t residue = (integer + 64) % period;
    if (residue < 0) residue += period;
    return static_cast<long double>((residue >> digit) & 1);
}

long double centerPreservingBridge(long double x, std::size_t digit) {
    const auto left = static_cast<std::int64_t>(std::floor(x));
    const long double fraction = x - static_cast<long double>(left);
    const long double a = digitAtInteger(left, digit);
    const long double b = digitAtInteger(left + 1, digit);
    const long double blend = (1 - std::cos(std::acos(-1.0L) * fraction)) / 2;
    return a + (b - a) * blend;
}

std::string decimal(long double value) {
    if (std::abs(value) < 1e-30L) value = 0;
    std::ostringstream out;
    out << std::scientific << std::setprecision(std::numeric_limits<long double>::max_digits10)
        << value;
    return out.str();
}

EvalModPolynomial chebyshevCandidate(std::size_t digit, std::size_t degree) {
    const std::size_t count = degree + 1;
    std::vector<long double> samples(count);
    const long double pi = std::acos(-1.0L);
    for (std::size_t k = 0; k < count; ++k) {
        const long double theta = pi * (static_cast<long double>(k) + .5L) / count;
        samples[k] = centerPreservingBridge(64 * std::cos(theta), digit);
    }
    EvalModPolynomial polynomial;
    polynomial.basis = PolynomialBasis::Chebyshev;
    polynomial.decimalCoefficients.resize(count);
    for (std::size_t order = 0; order < count; ++order) {
        long double sum = 0;
        for (std::size_t k = 0; k < count; ++k) {
            const long double theta = pi * (static_cast<long double>(k) + .5L) / count;
            sum += samples[k] * std::cos(order * theta);
        }
        const long double coefficient = sum * (order == 0 ? 1.0L : 2.0L) / count;
        polynomial.decimalCoefficients[order] = decimal(coefficient);
    }
    return polynomial;
}

long double clenshaw(const EvalModPolynomial& polynomial, long double x) {
    long double next = 0, nextNext = 0;
    for (std::size_t k = polynomial.decimalCoefficients.size(); k-- > 1;) {
        const long double coefficient = std::stold(polynomial.decimalCoefficients[k]);
        const long double current = 2 * x * next - nextNext + coefficient;
        nextNext = next;
        next = current;
    }
    return x * next - nextNext + std::stold(polynomial.decimalCoefficients[0]);
}

double gridError(const EvalModPolynomial& chebyshev, const EvalRoundProblem& problem,
                 std::size_t digit, std::size_t points) {
    long double maximum = 0;
    for (std::int64_t integer = -static_cast<std::int64_t>(problem.K);
         integer <= static_cast<std::int64_t>(problem.K); ++integer) {
        const auto target = evalRoundIntegerDigits(integer, problem.K,
                                                   EvalRoundRadix::Binary)[digit];
        for (std::size_t sample = 0; sample < points; ++sample) {
            const long double fraction = points == 1 ? 0.0L
                : -1.0L + 2.0L * sample / static_cast<long double>(points - 1);
            const long double x = static_cast<long double>(integer)
                + static_cast<long double>(problem.rho) * fraction;
            maximum = std::max(maximum,
                std::abs(clenshaw(chebyshev, x / 64) - target));
        }
    }
    return static_cast<double>(maximum);
}

double maximumCoefficient(const EvalModPolynomial& polynomial) {
    double maximum = 0;
    for (const auto& coefficient : polynomial.decimalCoefficients) {
        const double magnitude = std::abs(parseExactDecimal(coefficient).get_d());
        maximum = std::max(maximum, magnitude);
    }
    return maximum;
}

bool supportedFamily(EvalModApproximationFamily family) {
    return family == EvalModApproximationFamily::MultiIntervalMinimax
        || family == EvalModApproximationFamily::MultiIntervalChebyshev;
}

} // namespace

EvalRoundBinaryDigitSearchResult searchEvalRoundBinaryDigitPolynomials(
    const EvalRoundProblem& problem, const EvalRoundBinaryDigitSearchConfig& config) {
    if (problem.K != 64 || evalRoundDigitCount(problem.K, EvalRoundRadix::Binary) != 8
        || !std::isfinite(problem.rho) || problem.rho < 0 || problem.rho >= .5
        || config.degrees.empty() || config.families.empty()
        || config.gridPointsPerInterval == 0 || config.intervalSubdivisions < 2
        || config.intervalSubdivisions > 4096 || config.remezSamplesPerInterval < 2
        || config.remezMaximumIterations == 0 || config.remezMaximumIterations > 24)
        throw std::invalid_argument("invalid bounded K=64 binary digit search");
    for (const auto degree : config.degrees)
        if (degree == 0 || degree > 256)
            throw std::invalid_argument("digit search degree is outside generic compiler limit");
    for (const auto family : config.families)
        if (!supportedFamily(family))
            throw std::invalid_argument("unsupported digit approximation family");

    EvalRoundBinaryDigitSearchResult result;
    result.problem = problem;
    result.config = config;
    for (std::int64_t integer = -64; integer <= 64; ++integer)
        result.targetTable.push_back(evalRoundIntegerDigits(
            integer, 64, EvalRoundRadix::Binary));
    result.selectedRecordByDigit.resize(8);

    for (std::size_t digit = 0; digit < 8; ++digit) {
        for (const auto family : config.families) for (const auto degree : config.degrees) {
            EvalRoundDigitSearchRecord record;
            record.digitIndex = digit;
            record.family = family;
            record.requestedDegree = degree;
            if (family == EvalModApproximationFamily::MultiIntervalMinimax) {
                try {
                    MultiIntervalRemezRequest request;
                    request.degree = degree;
                    request.basis = PolynomialBasis::Chebyshev;
                    request.variableScaleDecimal = "64";
                    request.samplesPerInterval = config.remezSamplesPerInterval;
                    request.maximumIterations = config.remezMaximumIterations;
                    for (std::int64_t integer = -64; integer <= 64; ++integer) {
                        const mpq_class center(static_cast<long>(integer));
                        const mpq_class rho(problem.rho);
                        request.intervals.push_back({
                            exactRationalTerminatingDecimal(center - rho),
                            exactRationalTerminatingDecimal(center + rho),
                            std::to_string(evalRoundIntegerDigits(
                                integer, 64, EvalRoundRadix::Binary)[digit])});
                    }
                    const auto remez = generateMultiIntervalRemez(request);
                    auto certificate = certifyEvalRoundDigitPolynomial(
                        problem, digit, remez.polynomial, config.intervalSubdivisions);
                    certificate.executionRepresentation = EvalRoundPolynomialExecutionRepresentation{
                        remez.executionPolynomial,remez.executionVariableScaleDecimal,
                        "Exact rational scaled-Chebyshev conversion equals the canonical certified monomial polynomial"};
                    record.generatorStatus = EvalRoundDigitGeneratorStatus::Generated;
                    record.generatorConverged = remez.converged;
                    record.exchangeIterations = remez.exchangeIterations;
                    record.exchangePointsInsideDomain = remez.exchangePointsInsideDomain;
                    record.gridMaximumError = remez.sampledMaximumError;
                    record.directXRigorousIntervalError = certificate.directXApproximationError.upperBound;
                    record.centeredRigorousIntervalError = certificate.centeredApproximationError.upperBound;
                    record.rigorousIntervalError = certificate.approximationError.upperBound;
                    record.selectedIntervalProofMethod = certificate.selectedIntervalProofMethod;
                    record.maximumCoefficientMagnitude = maximumCoefficient(remez.polynomial);
                    record.cleanerInputDomainSatisfied = std::isfinite(record.rigorousIntervalError)
                        && record.rigorousIntervalError <= 1;
                    record.provenance = "Generalized existing MPFR-384 Remez exchange over only the 129 exact closed digit intervals; full Chebyshev basis in x/64; convergence and sampled error are diagnostic only";
                    certificate.provenance += "; candidate generator: " + record.provenance;
                    record.certificate = std::move(certificate);
                } catch (const std::exception& error) {
                    record.generatorStatus = EvalRoundDigitGeneratorStatus::GenerationFailed;
                    record.provenance = std::string("Remez generation/certification failed: ")
                        + error.what();
                }
                result.records.push_back(std::move(record));
                continue;
            }
            try {
                const auto chebyshev = chebyshevCandidate(digit, degree);
                record.gridMaximumError = gridError(
                    chebyshev, problem, digit, config.gridPointsPerInterval);
                const auto monomial = convertScaledChebyshevToMonomial(chebyshev, "64");
                auto certificate = certifyEvalRoundDigitPolynomial(
                    problem, digit, monomial, config.intervalSubdivisions);
                certificate.executionRepresentation = EvalRoundPolynomialExecutionRepresentation{
                    chebyshev,"64",
                    "Exact rational scaled-Chebyshev conversion equals the canonical certified monomial polynomial"};
                record.generatorStatus = EvalRoundDigitGeneratorStatus::Generated;
                record.generatorConverged = true;
                record.exchangeIterations = 0;
                record.exchangePointsInsideDomain = false;
                record.directXRigorousIntervalError = certificate.directXApproximationError.upperBound;
                record.centeredRigorousIntervalError = certificate.centeredApproximationError.upperBound;
                record.rigorousIntervalError = certificate.approximationError.upperBound;
                record.selectedIntervalProofMethod = certificate.selectedIntervalProofMethod;
                record.maximumCoefficientMagnitude = maximumCoefficient(monomial);
                record.cleanerInputDomainSatisfied = std::isfinite(record.rigorousIntervalError)
                    && record.rigorousIntervalError <= 1;
                record.provenance = "MultiIntervalChebyshev DCT of an explicit C1 cosine bridge that equals bit_j(I+64) at every integer center; grid is diagnostic only; acceptance uses the outward digit verifier";
                certificate.provenance += "; candidate generator: " + record.provenance;
                record.certificate = std::move(certificate);
            } catch (const std::exception& error) {
                record.generatorStatus = EvalRoundDigitGeneratorStatus::GenerationFailed;
                record.provenance = std::string("candidate generation/certification failed: ") + error.what();
            }
            result.records.push_back(std::move(record));
        }
    }
    for (std::size_t index = 0; index < result.records.size(); ++index) {
        const auto& record = result.records[index];
        if (!record.certificate || !std::isfinite(record.rigorousIntervalError)) continue;
        auto& selected = result.selectedRecordByDigit[record.digitIndex];
        if (!selected || record.rigorousIntervalError
            < result.records[*selected].rigorousIntervalError) selected = index;
    }
    result.allDigitsCertified = std::all_of(
        result.selectedRecordByDigit.begin(), result.selectedRecordByDigit.end(),
        [](const auto& selected) { return selected.has_value(); });
    result.allCleanerInputDomainsSatisfied = result.allDigitsCertified
        && std::all_of(result.selectedRecordByDigit.begin(), result.selectedRecordByDigit.end(),
            [&](const auto& selected) { return result.records[*selected].cleanerInputDomainSatisfied; });
    result.extractorCertified = result.allDigitsCertified
        && result.allCleanerInputDomainsSatisfied;
    result.status = result.extractorCertified
        ? "All eight digits have deterministic whole-domain certificates inside a<=1"
        : result.allDigitsCertified
            ? "No Certified direct multi-interval Remez extractor was found in this bounded degree search."
            : "One or more digits lack a finite outward bound; no complete certified K=64 extractor in this bounded search; this is not a global impossibility result";
    return result;
}

EvalRoundCandidate makeEvalRoundBinaryPolynomialCandidate(
    const EvalRoundBinaryDigitSearchResult& search) {
    EvalRoundCandidate candidate;
    candidate.id = "k64_bounded_binary_polynomial_search";
    candidate.radix = EvalRoundRadix::Binary;
    candidate.extraction.method = EvalRoundExtractionMethod::ExternalPolynomial;
    candidate.extraction.description = "Bounded K=64 binary digit polynomial search";
    candidate.extraction.certifiedK = search.problem.K;
    candidate.extraction.certifiedRho = search.problem.rho;
    candidate.extraction.verified = search.extractorCertified;
    candidate.extraction.provenance = search.status;
    candidate.digits.resize(8);
    if (!search.extractorCertified) return candidate;
    for (std::size_t digit = 0; digit < 8; ++digit) {
        const auto& record = search.records[*search.selectedRecordByDigit[digit]];
        candidate.extraction.polynomials.push_back(*record.certificate);
        candidate.extraction.polynomialDegrees.push_back(record.requestedDegree);
        candidate.digits[digit].extractionError = record.certificate->approximationError;
        candidate.digits[digit].cleaningLocalErrors.assign(
            search.problem.maxCleaningRoundsPerDigit,
            BootstrapBound{0, BootstrapBoundKind::Deterministic,
                "Exact mathematical cleaner diagnostic; not a ciphertext arithmetic bound", {}});
        candidate.digits[digit].reconstructionLocalError = {
            0, BootstrapBoundKind::Deterministic,
            "Exact mathematical binary reconstruction diagnostic; not a ciphertext arithmetic bound", {}};
    }
    return candidate;
}

} // namespace m2424::experimental
