#include "m2424/experimental/evalmod_analysis/synthesis.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace m2424::experimental {
namespace {

std::vector<long double> fit(const EvalModDomain& domain, std::size_t degree, bool sineTarget) {
    const std::size_t count = std::max<std::size_t>(256, (degree + 1) * 16);
    std::vector<std::vector<long double>> matrix(degree + 1, std::vector<long double>(degree + 2));
    for (std::int64_t k = -static_cast<std::int64_t>(domain.integerBound);
         k <= static_cast<std::int64_t>(domain.integerBound); ++k) {
        for (std::size_t sample = 0; sample < count; ++sample) {
            const long double r = -domain.normalizedResidualBound
                + 2.0L * domain.normalizedResidualBound * sample / (count - 1);
            const long double x = k + r;
            const long double y = sineTarget ? std::sin(2.0L * std::acos(-1.0L) * x)
                                                   / (2.0L * std::acos(-1.0L)) : r;
            std::vector<long double> powers(2 * degree + 1, 1.0L);
            for (std::size_t i = 1; i < powers.size(); ++i) powers[i] = powers[i - 1] * x;
            for (std::size_t row = 0; row <= degree; ++row) {
                for (std::size_t col = 0; col <= degree; ++col) matrix[row][col] += powers[row + col];
                matrix[row][degree + 1] += y * powers[row];
            }
        }
    }
    for (std::size_t pivot = 0; pivot <= degree; ++pivot) {
        std::size_t best = pivot;
        for (std::size_t row = pivot + 1; row <= degree; ++row)
            if (std::abs(matrix[row][pivot]) > std::abs(matrix[best][pivot])) best = row;
        if (std::abs(matrix[best][pivot]) < 1e-30L) throw std::runtime_error("singular approximation fit");
        std::swap(matrix[pivot], matrix[best]);
        const long double divisor = matrix[pivot][pivot];
        for (std::size_t col = pivot; col <= degree + 1; ++col) matrix[pivot][col] /= divisor;
        for (std::size_t row = 0; row <= degree; ++row) if (row != pivot) {
            const long double factor = matrix[row][pivot];
            for (std::size_t col = pivot; col <= degree + 1; ++col)
                matrix[row][col] -= factor * matrix[pivot][col];
        }
    }
    std::vector<long double> result(degree + 1);
    for (std::size_t i = 0; i <= degree; ++i) result[i] = matrix[i][degree + 1];
    return result;
}

EvalModPolynomial polynomial(const std::vector<long double>& coefficients) {
    EvalModPolynomial result;
    result.basis = PolynomialBasis::Monomial;
    for (long double value : coefficients) {
        std::ostringstream text;
        text << std::setprecision(std::numeric_limits<long double>::max_digits10) << value;
        result.decimalCoefficients.push_back(text.str());
    }
    return result;
}

std::size_t depthForDegree(std::size_t degree) {
    std::size_t depth = 0, power = 1;
    while (power < degree) { power <<= 1; ++depth; }
    return depth + 1;
}

const char* familyName(EvalModApproximationFamily family) {
    switch (family) {
    case EvalModApproximationFamily::PeriodicSineBaseline: return "periodic_sine_baseline";
    case EvalModApproximationFamily::MultiIntervalLeastSquaresPrototype:
        return "multi_interval_least_squares_prototype";
    case EvalModApproximationFamily::MultiIntervalMinimax: return "multi_interval_minimax";
    case EvalModApproximationFamily::MinimaxInverseSine: return "minimax_inverse_sine";
    }
    return "unknown";
}

std::string jsonEscape(const std::string& value) {
    std::string result;
    for (char character : value) {
        if (character == '\\' || character == '"') result.push_back('\\');
        if (character == '\n') result += "\\n";
        else if (character != '\r') result.push_back(character);
    }
    return result;
}

} // namespace

EvalModSynthesisResult synthesizeEvalMod(const EvalModProblem& problem) {
    if (problem.qSource <= 1 || problem.availableLevels == 0 || problem.targetPrecisionBits < 32
        || !std::isfinite(problem.targetAbsoluteError) || problem.targetAbsoluteError <= 0.0) {
        throw std::invalid_argument("invalid EvalMod synthesis problem");
    }
    EvalModSynthesisResult result{problem, estimateEvalModDomain(problem.ciphertextModel), {}, {}};
    const std::size_t degrees[]{9, 15};
    const EvalModApproximationFamily families[]{EvalModApproximationFamily::PeriodicSineBaseline,
                                                EvalModApproximationFamily::MultiIntervalLeastSquaresPrototype};
    for (std::size_t index = 0; index < 2; ++index) {
        EvalModCandidate candidate;
        candidate.family = families[index];
        candidate.domain = result.domain;
        candidate.polynomial = polynomial(fit(result.domain, degrees[index], index == 0));
        const std::size_t depth = depthForDegree(degrees[index]);
        candidate.circuit = {degrees[index], depth, depth, degrees[index], degrees[index], depth,
                             depth, depth, 3, (degrees[index] + 1) * sizeof(double), 4096};
        const std::size_t scaleBits = problem.targetPrecisionBits + 8;
        for (std::size_t level = 0; level < depth; ++level)
            candidate.scaleSchedule.push_back({scaleBits, scaleBits, scaleBits, 8});
        candidate.diagnostic = diagnoseEvalModPolynomialOnGrid(candidate.polynomial, result.domain,
                                                               129, "1e-12", 0.0,
                                                               std::max<std::size_t>(256, problem.targetPrecisionBits * 4));
        const ExactInteger gainNumerator = problem.qSource * problem.outputScale.denominator;
        const double gain = gainNumerator.get_d() / problem.outputScale.numerator.get_d();
        candidate.propagationBounds = {problem.ciphertextModel.normalizedCoeffToSlotErrorAbsBound,
            candidate.diagnostic.approximationMaxError, 0.0, 0.0,
            problem.ciphertextModel.relativePeriodMismatchAbsBound,
            candidate.diagnostic.complexDerivativeMax, gain, problem.slotToCoeffOperatorNorm,
            problem.slotToCoeffAdditive, problem.finalAdditive};
        candidate.predictedBootstrapError = propagatedBootstrapError(candidate.propagationBounds);
        candidate.cost = estimateEvalModCost(candidate.circuit, problem.backendCost, scaleBits);
        candidate.feasible = candidate.circuit.levelConsumption <= problem.availableLevels
            && candidate.predictedBootstrapError <= problem.targetAbsoluteError;
        result.candidates.push_back(std::move(candidate));
    }
    for (std::size_t i = 0; i < result.candidates.size(); ++i) {
        if (result.candidates[i].feasible && (!result.selectedCandidate
            || result.candidates[i].cost.latencyMs
                < result.candidates[*result.selectedCandidate].cost.latencyMs)) result.selectedCandidate = i;
    }
    return result;
}

double evaluateEvalModPolynomial(const EvalModPolynomial& polynomial, double input) {
    if (polynomial.basis != PolynomialBasis::Monomial) throw std::invalid_argument("unsupported basis");
    long double value = 0.0L;
    for (auto it = polynomial.decimalCoefficients.rbegin(); it != polynomial.decimalCoefficients.rend(); ++it)
        value = value * input + std::stold(*it);
    return static_cast<double>(value);
}

std::string evalModSynthesisCsv(const EvalModSynthesisResult& result) {
    std::ostringstream out;
    out << "family,K,rho,degree,depth,ct_ct,ct_pt,rescales,modulus_bits,real_error,complex_error,"
           "real_derivative,complex_derivative,predicted_error,latency_ms,feasible,selected\n";
    for (std::size_t i = 0; i < result.candidates.size(); ++i) {
        const auto& c = result.candidates[i];
        out << familyName(c.family) << ',' << c.domain.integerBound << ',' << std::setprecision(17)
            << c.domain.normalizedResidualBound << ',' << c.circuit.degree << ','
            << c.circuit.multiplicativeDepth << ',' << c.circuit.ciphertextMultiplications << ','
            << c.circuit.ciphertextPlaintextMultiplications << ',' << c.circuit.rescales << ','
            << c.cost.dataModulusBitsLowerBound << ',' << c.diagnostic.approximationMaxError << ','
            << c.diagnostic.complexBoundaryErrorMax << ',' << c.diagnostic.realDerivativeMax << ','
            << c.diagnostic.complexDerivativeMax << ',' << c.predictedBootstrapError << ','
            << c.cost.latencyMs << ',' << c.feasible << ','
            << (result.selectedCandidate && *result.selectedCandidate == i) << '\n';
    }
    return out.str();
}

std::string evalModSynthesisJson(const EvalModSynthesisResult& result) {
    std::ostringstream out;
    out << "{\"K\":" << result.domain.integerBound << ",\"rho\":"
        << std::setprecision(17) << result.domain.normalizedResidualBound << ",\"selected\":";
    if (result.selectedCandidate) out << *result.selectedCandidate; else out << "null";
    out << ",\"tail_model_provenance\":{\"derivation\":\""
        << jsonEscape(result.problem.ciphertextModel.provenance.derivation)
        << "\",\"included_noise_sources\":\""
        << jsonEscape(result.problem.ciphertextModel.provenance.includedNoiseSources)
        << "\",\"assumptions\":\"" << jsonEscape(result.problem.ciphertextModel.provenance.assumptions)
        << "\"},\"candidates\":[";
    for (std::size_t i = 0; i < result.candidates.size(); ++i) {
        if (i) out << ',';
        const auto& c = result.candidates[i];
        out << "{\"family\":\"" << familyName(c.family) << "\",\"degree\":" << c.circuit.degree
            << ",\"depth\":" << c.circuit.multiplicativeDepth << ",\"real_error\":"
            << c.diagnostic.approximationMaxError << ",\"complex_error\":"
            << c.diagnostic.complexBoundaryErrorMax << ",\"complex_derivative\":"
            << c.diagnostic.complexDerivativeMax << ",\"predicted_error\":"
            << c.predictedBootstrapError << ",\"feasible\":" << (c.feasible ? "true" : "false") << '}';
    }
    out << "]}";
    return out.str();
}

} // namespace m2424::experimental
