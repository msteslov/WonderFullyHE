#include "m2424/experimental/evalmod_analysis/synthesis.hpp"

#include <mpfr.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace m2424::experimental {
namespace {

class MpReal {
public:
    explicit MpReal(mpfr_prec_t precision = 384) { mpfr_init2(value_, precision); }
    MpReal(const MpReal& other) { mpfr_init2(value_, mpfr_get_prec(other.value_)); mpfr_set(value_, other.value_, MPFR_RNDN); }
    MpReal(MpReal&& other) noexcept { mpfr_init2(value_, mpfr_get_prec(other.value_)); mpfr_swap(value_, other.value_); }
    MpReal& operator=(const MpReal& other) { mpfr_set(value_, other.value_, MPFR_RNDN); return *this; }
    ~MpReal() { mpfr_clear(value_); }
    mpfr_ptr get() { return value_; } mpfr_srcptr get() const { return value_; }
private: mpfr_t value_;
};

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

EvalModPolynomial remezOdd(const EvalModDomain& domain, std::size_t degree) {
    const std::size_t terms = (degree + 1) / 2;
    const std::size_t unknowns = terms + 1;
    std::vector<long double> grid;
    for (std::size_t integer = 0; integer <= domain.integerBound; ++integer) {
        for (std::size_t sample = 0; sample < 2048; ++sample) {
            long double residual = -domain.normalizedResidualBound
                + 2.0L * domain.normalizedResidualBound * sample / 2047.0L;
            const long double x = integer + residual;
            if (x >= 0.0L) grid.push_back(x);
        }
    }
    std::vector<std::size_t> extrema(unknowns);
    for (std::size_t i = 0; i < unknowns; ++i) extrema[i] = i * (grid.size() - 1) / (unknowns - 1);
    std::vector<MpReal> solution(unknowns);
    bool converged = false;
    long double previousMaximum = std::numeric_limits<long double>::infinity();
    for (std::size_t iteration = 0; iteration < 24; ++iteration) {
        std::vector<std::vector<MpReal>> a;
        a.reserve(unknowns);
        for (std::size_t row = 0; row < unknowns; ++row) {
            a.emplace_back(unknowns + 1);
            const long double x = grid[extrema[row]];
            const long double residual = x - std::floor(x + 0.5L);
            MpReal base, power;
            mpfr_set_ld(base.get(), x, MPFR_RNDN); mpfr_set(base.get(), base.get(), MPFR_RNDN);
            for (std::size_t col = 0; col < terms; ++col) {
                mpfr_pow_ui(a[row][col].get(), base.get(), 2 * col + 1, MPFR_RNDN);
            }
            mpfr_set_si(a[row][terms].get(), row % 2 ? -1 : 1, MPFR_RNDN);
            mpfr_set_ld(a[row][unknowns].get(), residual, MPFR_RNDN);
        }
        for (std::size_t pivot = 0; pivot < unknowns; ++pivot) {
            std::size_t best = pivot;
            for (std::size_t row = pivot + 1; row < unknowns; ++row)
                if (mpfr_cmpabs(a[row][pivot].get(), a[best][pivot].get()) > 0) best = row;
            std::swap(a[pivot], a[best]);
            if (mpfr_zero_p(a[pivot][pivot].get())) throw std::runtime_error("singular Remez exchange");
            for (std::size_t col = pivot + 1; col <= unknowns; ++col)
                mpfr_div(a[pivot][col].get(), a[pivot][col].get(), a[pivot][pivot].get(), MPFR_RNDN);
            mpfr_set_ui(a[pivot][pivot].get(), 1, MPFR_RNDN);
            for (std::size_t row = 0; row < unknowns; ++row) if (row != pivot) {
                MpReal factor; mpfr_set(factor.get(), a[row][pivot].get(), MPFR_RNDN);
                for (std::size_t col = pivot + 1; col <= unknowns; ++col) {
                    MpReal product; mpfr_mul(product.get(), factor.get(), a[pivot][col].get(), MPFR_RNDN);
                    mpfr_sub(a[row][col].get(), a[row][col].get(), product.get(), MPFR_RNDN);
                }
                mpfr_set_zero(a[row][pivot].get(), 0);
            }
        }
        for (std::size_t i = 0; i < unknowns; ++i) mpfr_set(solution[i].get(), a[i][unknowns].get(), MPFR_RNDN);
        long double globalMaximum = 0.0L;
        for (std::size_t segment = 0; segment < unknowns; ++segment) {
            const std::size_t begin = segment * grid.size() / unknowns;
            const std::size_t end = (segment + 1) * grid.size() / unknowns;
            long double maximum = -1.0L;
            for (std::size_t point = begin; point < end; ++point) {
                long double value = 0.0L;
                for (std::size_t col = terms; col-- > 0;)
                    value = value * grid[point] * grid[point] + mpfr_get_ld(solution[col].get(), MPFR_RNDN);
                value *= grid[point];
                const long double target = grid[point] - std::floor(grid[point] + 0.5L);
                if (std::abs(value - target) > maximum) {
                    maximum = std::abs(value - target); extrema[segment] = point;
                }
            }
            globalMaximum = std::max(globalMaximum, maximum);
        }
        if (std::isfinite(previousMaximum)
            && std::abs(globalMaximum - previousMaximum)
                <= 5e-2L * std::max(1e-12L, globalMaximum)) { converged = true; break; }
        previousMaximum = globalMaximum;
    }
    // При достижении iteration budget принимаем только конечную bounded exchange sequence;
    // независимый interval certificate ниже остаётся окончательным gate для кандидата.
    if (!converged) converged = std::isfinite(previousMaximum) && previousMaximum >= 0.0L;
    if (!converged) throw std::runtime_error("MPFR Remez exchange did not converge");
    EvalModPolynomial result;
    result.basis = PolynomialBasis::Monomial;
    result.decimalCoefficients.assign(degree + 1, "0");
    for (std::size_t index = 0; index < terms; ++index) {
        char* text = nullptr;
        mpfr_asprintf(&text, "%.120Rg", solution[index].get());
        result.decimalCoefficients[2 * index + 1] = text ? text : "0";
        mpfr_free_str(text);
    }
    return result;
}

std::vector<long double> multiplyPolynomial(const std::vector<long double>& left,
                                            const std::vector<long double>& right) {
    std::vector<long double> result(left.size() + right.size() - 1, 0.0L);
    for (std::size_t i = 0; i < left.size(); ++i)
        for (std::size_t j = 0; j < right.size(); ++j) result[i + j] += left[i] * right[j];
    return result;
}

std::vector<long double> inverseSineCorrection(const EvalModDomain& domain) {
    const auto sine = fit(domain, 7, true);
    const auto square = multiplyPolynomial(sine, sine);
    const auto cube = multiplyPolynomial(square, sine);
    const auto fifth = multiplyPolynomial(cube, square);
    std::vector<long double> result(fifth.size(), 0.0L);
    const long double pi = std::acos(-1.0L);
    for (std::size_t i = 0; i < sine.size(); ++i) result[i] += sine[i];
    for (std::size_t i = 0; i < cube.size(); ++i) result[i] += (4.0L * pi * pi / 6.0L) * cube[i];
    for (std::size_t i = 0; i < fifth.size(); ++i)
        result[i] += (3.0L * 16.0L * pi * pi * pi * pi / 40.0L) * fifth[i];
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

ExactScale multiplyScale(const ExactScale& left, const ExactScale& right) {
    return ExactScale::rational(left.numerator * right.numerator,
                                left.denominator * right.denominator);
}

double exactScaleUp(const ExactScale& scale) {
    mpfr_t value;
    mpfr_init2(value, 256);
    mpq_class rational(scale.numerator, scale.denominator);
    rational.canonicalize();
    mpfr_set_q(value, rational.get_mpq_t(), MPFR_RNDU);
    const double result = mpfr_get_d(value, MPFR_RNDU);
    mpfr_clear(value);
    if (!std::isfinite(result)) throw std::overflow_error("EvalMod gain exceeds planner range");
    return result;
}

std::size_t exactScaleBitsUp(const ExactScale& scale) {
    const double value = exactScaleUp(scale);
    if (value <= 0.0) throw std::invalid_argument("non-positive EvalMod scale");
    return static_cast<std::size_t>(std::max(0.0, std::ceil(std::log2(value))));
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

const char* stageName(EvalModCandidateStage stage) {
    switch (stage) {
    case EvalModCandidateStage::Generated: return "generated";
    case EvalModCandidateStage::GridDiagnosed: return "grid_diagnosed";
    case EvalModCandidateStage::IntervalCertified: return "interval_certified";
    case EvalModCandidateStage::CircuitCompiled: return "circuit_compiled";
    case EvalModCandidateStage::ScaleScheduled: return "scale_scheduled";
    case EvalModCandidateStage::BackendValidated: return "backend_validated";
    }
    return "unknown";
}

const char* rejectionName(EvalModRejectionReason reason) {
    switch (reason) {
    case EvalModRejectionReason::None: return "none";
    case EvalModRejectionReason::ApproximationError: return "approximation_error";
    case EvalModRejectionReason::Uncertified: return "uncertified";
    case EvalModRejectionReason::ArithmeticErrorUnknown: return "arithmetic_error_unknown";
    case EvalModRejectionReason::InsufficientLevels: return "insufficient_levels";
    case EvalModRejectionReason::ModulusBudget: return "modulus_budget";
    case EvalModRejectionReason::SecurityBudget: return "security_budget";
    case EvalModRejectionReason::ScaleScheduleFailure: return "scale_schedule_failure";
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

CompiledEvalModCircuit compileEvalModPolynomial(const EvalModPolynomial& polynomial,
                                               const EvalModProblem& problem,
                                               std::size_t babyStep) {
    if (polynomial.basis != PolynomialBasis::Monomial || polynomial.decimalCoefficients.empty()
        || babyStep < 2) throw std::invalid_argument("invalid polynomial compiler input");
    CompiledEvalModCircuit compiled;
    compiled.babyStep = babyStep;
    compiled.normalizationGain = ExactScale::rational(
        problem.coeffToSlotScale.numerator,
        problem.coeffToSlotScale.denominator * problem.qSource);
    compiled.denormalizationGain = ExactScale::rational(
        problem.qSource * problem.outputScale.denominator, problem.outputScale.numerator);
    const ExactScale workingScale = problem.coeffToSlotScale;
    auto addNode = [&](EvalModOperation operation, std::vector<std::size_t> inputs,
                       std::size_t chain, ExactScale input, ExactScale output,
                       std::string constant = {}) {
        compiled.nodes.push_back({operation, std::move(inputs), chain, std::move(input),
                                  std::move(output), std::move(constant)});
        return compiled.nodes.size() - 1;
    };
    const std::size_t input = addNode(EvalModOperation::Input, {}, 0, workingScale, workingScale);
    const std::size_t normalization = addNode(EvalModOperation::EncodeConstant, {}, 0,
                                              workingScale, workingScale);
    std::size_t x = addNode(EvalModOperation::MultiplyPlain, {input, normalization}, 0,
                            workingScale, multiplyScale(workingScale, workingScale));
    x = addNode(EvalModOperation::Rescale, {x}, 1, compiled.nodes[x].outputScale, workingScale);

    const std::size_t degree = polynomial.decimalCoefficients.size() - 1;
    std::vector<std::size_t> powers(babyStep + 1, x);
    std::vector<std::size_t> depths(compiled.nodes.size(), 0);
    depths.resize(compiled.nodes.size(), 0);
    for (std::size_t power = 2; power <= babyStep; ++power) {
        auto multiplied = addNode(EvalModOperation::MultiplyCipher, {powers[power - 1], x}, power - 1,
                                  workingScale, multiplyScale(workingScale, workingScale));
        depths.resize(compiled.nodes.size()); depths[multiplied] = depths[powers[power - 1]] + 1;
        auto relin = addNode(EvalModOperation::Relinearize, {multiplied}, power - 1,
                             compiled.nodes[multiplied].outputScale, compiled.nodes[multiplied].outputScale);
        depths.resize(compiled.nodes.size()); depths[relin] = depths[multiplied];
        powers[power] = addNode(EvalModOperation::Rescale, {relin}, power, compiled.nodes[relin].outputScale,
                                workingScale);
        depths.resize(compiled.nodes.size()); depths[powers[power]] = depths[relin];
    }
    std::vector<std::size_t> giantPowers{input, powers[babyStep]};
    const std::size_t blocks = degree / babyStep + 1;
    for (std::size_t block = 2; block < blocks; ++block) {
        auto mul = addNode(EvalModOperation::MultiplyCipher, {giantPowers.back(), powers[babyStep]}, block,
                           workingScale, multiplyScale(workingScale, workingScale));
        depths.resize(compiled.nodes.size()); depths[mul] = depths[giantPowers.back()] + 1;
        auto relin = addNode(EvalModOperation::Relinearize, {mul}, block,
                             compiled.nodes[mul].outputScale, compiled.nodes[mul].outputScale);
        depths.resize(compiled.nodes.size()); depths[relin] = depths[mul];
        auto rescale = addNode(EvalModOperation::Rescale, {relin}, block + 1,
                               compiled.nodes[relin].outputScale, workingScale);
        depths.resize(compiled.nodes.size()); depths[rescale] = depths[relin];
        giantPowers.push_back(rescale);
    }
    std::optional<std::size_t> total;
    for (std::size_t block = 0; block < blocks; ++block) {
        std::optional<std::size_t> blockValue;
        for (std::size_t offset = 0; offset < babyStep; ++offset) {
            const std::size_t coefficient = block * babyStep + offset;
            if (coefficient > degree) break;
            auto constant = addNode(EvalModOperation::EncodeConstant, {}, block, workingScale, workingScale,
                                    polynomial.decimalCoefficients[coefficient]);
            std::size_t term = constant;
            if (offset > 0) {
                term = addNode(EvalModOperation::MultiplyPlain, {powers[offset], constant}, block,
                               workingScale, multiplyScale(workingScale, workingScale));
                depths.resize(compiled.nodes.size()); depths[term] = depths[powers[offset]];
                term = addNode(EvalModOperation::Rescale, {term}, block + 1,
                               compiled.nodes[term].outputScale, workingScale);
                depths.resize(compiled.nodes.size()); depths[term] = depths[powers[offset]];
            }
            if (blockValue) {
                const auto previous = *blockValue;
                const auto current = term;
                term = addNode(EvalModOperation::Add, {previous, current}, block + 1,
                               workingScale, workingScale);
                depths.resize(compiled.nodes.size());
                depths[term] = std::max(depths[previous], depths[current]);
            } else depths.resize(compiled.nodes.size());
            blockValue = term;
        }
        std::size_t term = *blockValue;
        if (block > 0) {
            auto mul = addNode(EvalModOperation::MultiplyCipher, {term, giantPowers[block]}, block + 1,
                               workingScale, multiplyScale(workingScale, workingScale));
            depths.resize(compiled.nodes.size()); depths[mul] = std::max(depths[term], depths[giantPowers[block]]) + 1;
            auto relin = addNode(EvalModOperation::Relinearize, {mul}, block + 1,
                                 compiled.nodes[mul].outputScale, compiled.nodes[mul].outputScale);
            depths.resize(compiled.nodes.size()); depths[relin] = depths[mul];
            term = addNode(EvalModOperation::Rescale, {relin}, block + 2,
                           compiled.nodes[relin].outputScale, workingScale);
            depths.resize(compiled.nodes.size()); depths[term] = depths[relin];
        }
        if (total) {
            const auto previous = *total;
            const auto current = term;
            term = addNode(EvalModOperation::Add, {previous, current}, block + 2,
                           workingScale, workingScale);
            depths.resize(compiled.nodes.size());
            depths[term] = std::max(depths[previous], depths[current]);
        } else depths.resize(compiled.nodes.size());
        total = term;
    }
    const auto denorm = addNode(EvalModOperation::EncodeConstant, {}, 0, problem.outputScale,
                                problem.outputScale);
    auto output = addNode(EvalModOperation::MultiplyPlain, {*total, denorm}, 0, workingScale,
                          multiplyScale(workingScale, problem.outputScale));
    depths.resize(compiled.nodes.size()); depths[output] = depths[*total];
    output = addNode(EvalModOperation::Rescale, {output}, 0, compiled.nodes[output].outputScale,
                     problem.outputScale);
    depths.resize(compiled.nodes.size()); depths[output] = depths[*total];
    compiled.outputNode = output;

    EvalModCircuitCost cost;
    cost.degree = degree;
    cost.multiplicativeDepth = *std::max_element(depths.begin(), depths.end());
    for (const auto& node : compiled.nodes) {
        cost.ciphertextMultiplications += node.operation == EvalModOperation::MultiplyCipher;
        cost.ciphertextPlaintextMultiplications += node.operation == EvalModOperation::MultiplyPlain;
        cost.additions += node.operation == EvalModOperation::Add;
        cost.relinearizations += node.operation == EvalModOperation::Relinearize;
        cost.rescales += node.operation == EvalModOperation::Rescale;
        cost.levelConsumption = std::max(cost.levelConsumption, node.chainIndex);
        cost.preparedPlaintextBytes += node.operation == EvalModOperation::EncodeConstant ? 64 : 0;
    }
    std::vector<std::size_t> remaining(compiled.nodes.size(), 0);
    for (const auto& node : compiled.nodes) for (auto source : node.inputs) ++remaining[source];
    std::vector<bool> live(compiled.nodes.size(), false);
    std::size_t liveCount = 0;
    for (std::size_t index = 0; index < compiled.nodes.size(); ++index) {
        const bool ciphertext = compiled.nodes[index].operation != EvalModOperation::EncodeConstant;
        if (ciphertext) { live[index] = true; ++liveCount; cost.peakLiveCiphertexts = std::max(cost.peakLiveCiphertexts, liveCount); }
        for (auto source : compiled.nodes[index].inputs) {
            if (--remaining[source] == 0 && live[source] && source != compiled.outputNode) {
                live[source] = false; --liveCount;
            }
        }
    }
    cost.backendScratchBytes = 4096;
    compiled.cost = cost;
    return compiled;
}

EvalModSynthesisResult synthesizeEvalMod(const EvalModProblem& problem) {
    if (problem.qSource <= 1 || problem.coeffToSlotScale.numerator <= 0
        || problem.coeffToSlotScale.denominator <= 0 || problem.outputScale.numerator <= 0
        || problem.outputScale.denominator <= 0 || problem.availableLevels == 0
        || problem.targetPrecisionBits < 32
        || !std::isfinite(problem.targetAbsoluteError) || problem.targetAbsoluteError <= 0.0) {
        throw std::invalid_argument("invalid EvalMod synthesis problem");
    }
    EvalModSynthesisResult result{problem, estimateEvalModDomain(problem.ciphertextModel), {}, {}};
    struct Spec { EvalModApproximationFamily family; std::size_t degree; std::size_t babyStep; };
    std::vector<Spec> specs{{EvalModApproximationFamily::PeriodicSineBaseline, 9, 3},
                            {EvalModApproximationFamily::MultiIntervalLeastSquaresPrototype, 15, 4}};
    for (std::size_t degree : std::vector<std::size_t>{7, 9, 11, 13, 15})
        for (std::size_t babyStep : std::vector<std::size_t>{2, 3, 4, 5})
            specs.push_back({EvalModApproximationFamily::MultiIntervalMinimax, degree, babyStep});
    specs.push_back({EvalModApproximationFamily::MinimaxInverseSine, 35, 5});
    for (const auto& spec : specs) {
        EvalModCandidate candidate;
        candidate.family = spec.family;
        candidate.domain = result.domain;
        candidate.polynomial = spec.family == EvalModApproximationFamily::MultiIntervalMinimax
            ? remezOdd(result.domain, spec.degree)
            : polynomial(spec.family == EvalModApproximationFamily::MinimaxInverseSine
                ? inverseSineCorrection(result.domain)
                : fit(result.domain, spec.degree,
                      spec.family == EvalModApproximationFamily::PeriodicSineBaseline));
        candidate.compiledCircuit = compileEvalModPolynomial(candidate.polynomial, problem, spec.babyStep);
        const std::size_t scaleBits = problem.targetPrecisionBits + 8;
        const std::size_t initialModulusBits = (std::max(candidate.compiledCircuit.cost.multiplicativeDepth,
            candidate.compiledCircuit.cost.levelConsumption) + 1) * scaleBits;
        for (const auto& node : candidate.compiledCircuit.nodes) {
            if (node.operation == EvalModOperation::Rescale) {
                const std::size_t inputBits = exactScaleBitsUp(node.inputScale);
                const std::size_t outputBits = exactScaleBitsUp(node.outputScale);
                const std::size_t primeBits = inputBits > outputBits ? inputBits - outputBits : scaleBits;
                candidate.scaleSchedule.push_back({inputBits, inputBits > outputBits
                    ? inputBits - outputBits : 0, primeBits, outputBits,
                    node.chainIndex * scaleBits >= initialModulusBits ? 0
                        : initialModulusBits - node.chainIndex * scaleBits, 8});
            }
        }
        candidate.diagnostic = diagnoseEvalModPolynomialOnGrid(candidate.polynomial, result.domain,
                                                               129, "1e-12", 0.0,
                                                               std::max<std::size_t>(256, problem.targetPrecisionBits * 4));
        candidate.intervalCertificate = certifyEvalModPolynomialIntervals(
            candidate.polynomial, result.domain, "1e-12", 512,
            std::max<std::size_t>(384, problem.targetPrecisionBits * 4));
        const double gain = exactScaleUp(candidate.compiledCircuit.denormalizationGain);
        candidate.propagationBounds = {problem.ciphertextModel.normalizedCoeffToSlotErrorAbsBound,
            candidate.intervalCertificate.approximationErrorUpperBoundDouble, 0.0, 0.0,
            result.domain.errors.periodMismatch,
            candidate.intervalCertificate.derivativeUpperBoundDouble, gain, problem.slotToCoeffOperatorNorm,
            problem.slotToCoeffAdditive, problem.finalAdditive};
        const auto& circuit = candidate.compiledCircuit.cost;
        const double operationCount = static_cast<double>(circuit.ciphertextMultiplications
            + circuit.ciphertextPlaintextMultiplications + circuit.relinearizations
            + circuit.rescales + circuit.additions);
        candidate.polynomialArithmeticError = std::ldexp(operationCount
            * (1.0 + candidate.intervalCertificate.derivativeUpperBoundDouble),
            -static_cast<int>(problem.targetPrecisionBits));
        candidate.propagationBounds.polynomialArithmetic = *candidate.polynomialArithmeticError;
        candidate.predictedBootstrapError = propagatedBootstrapError(candidate.propagationBounds);
        candidate.cost = estimateEvalModCost(candidate.compiledCircuit.cost, problem.backendCost, scaleBits);
        candidate.stage = EvalModCandidateStage::ScaleScheduled;
        candidate.satisfiesLevelBudget = candidate.compiledCircuit.cost.levelConsumption <= problem.availableLevels;
        candidate.intervalCertified = candidate.intervalCertificate.proved;
        candidate.executable = false;
        candidate.satisfiesNumericalTarget = candidate.intervalCertified
            && candidate.polynomialArithmeticError.has_value()
            && candidate.predictedBootstrapError <= problem.targetAbsoluteError;
        const bool scaleScheduleValid = std::all_of(candidate.scaleSchedule.begin(),
            candidate.scaleSchedule.end(), [](const EvalModScaleStage& stage) {
                return stage.availableModulusBits >= stage.outputScaleBits + stage.requiredHeadroomBits;
            });
        candidate.rejectionReason = !scaleScheduleValid
            ? EvalModRejectionReason::ScaleScheduleFailure
            : !candidate.satisfiesLevelBudget
            ? EvalModRejectionReason::InsufficientLevels
            : candidate.cost.dataModulusBitsLowerBound > problem.maxDataModulusBits
                ? EvalModRejectionReason::ModulusBudget
            : !candidate.intervalCertified ? EvalModRejectionReason::Uncertified
            : !candidate.satisfiesNumericalTarget ? EvalModRejectionReason::ApproximationError
            : EvalModRejectionReason::None;
        result.candidates.push_back(std::move(candidate));
    }
    for (std::size_t i = 0; i < result.candidates.size(); ++i) {
        const bool selectableFamily = result.candidates[i].family
            != EvalModApproximationFamily::MultiIntervalLeastSquaresPrototype;
        if (selectableFamily && result.candidates[i].satisfiesLevelBudget
            && result.candidates[i].satisfiesNumericalTarget
            && result.candidates[i].rejectionReason == EvalModRejectionReason::None
            && result.candidates[i].cost.latencyMs <= problem.maxLatencyMs
            && result.candidates[i].cost.peakWorkingSetBytes <= problem.maxWorkingSetBytes
            && (!result.provisionalSelection
            || result.candidates[i].cost.latencyMs
                < result.candidates[*result.provisionalSelection].cost.latencyMs)) result.provisionalSelection = i;
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
    out << "family,K,rho,degree,basis,strategy,baby_step,depth,ct_ct,ct_pt,rescales,level_consumption,"
           "modulus_bits,normalization_gain,denormalization_gain,real_error,complex_error,"
           "complex_derivative,interval_error,arithmetic_error_known,predicted_error,levels_ok,"
           "certified,executable,stage,rejection_reason,provisional\n";
    for (std::size_t i = 0; i < result.candidates.size(); ++i) {
        const auto& c = result.candidates[i];
        out << familyName(c.family) << ',' << c.domain.integerBound << ',' << std::setprecision(17)
            << c.domain.normalizedResidualBound << ',' << c.compiledCircuit.cost.degree
            << ",monomial,paterson_stockmeyer," << c.compiledCircuit.babyStep << ','
            << c.compiledCircuit.cost.multiplicativeDepth << ','
            << c.compiledCircuit.cost.ciphertextMultiplications << ','
            << c.compiledCircuit.cost.ciphertextPlaintextMultiplications << ','
            << c.compiledCircuit.cost.rescales << ',' << c.compiledCircuit.cost.levelConsumption << ','
            << c.cost.dataModulusBitsLowerBound << ',' << exactScaleUp(c.compiledCircuit.normalizationGain)
            << ',' << exactScaleUp(c.compiledCircuit.denormalizationGain) << ','
            << c.diagnostic.approximationMaxError << ','
            << c.diagnostic.complexBoundaryErrorMax << ',' << c.diagnostic.realDerivativeMax << ','
            << c.intervalCertificate.approximationErrorUpperBoundDouble << ','
            << c.polynomialArithmeticError.has_value() << ',' << c.predictedBootstrapError << ','
            << c.satisfiesLevelBudget << ',' << c.intervalCertified << ',' << c.executable << ','
            << stageName(c.stage) << ',' << rejectionName(c.rejectionReason) << ','
            << (result.provisionalSelection && *result.provisionalSelection == i) << '\n';
    }
    return out.str();
}

std::string evalModSynthesisJson(const EvalModSynthesisResult& result) {
    std::ostringstream out;
    out << "{\"K\":" << result.domain.integerBound << ",\"rho\":"
        << std::setprecision(17) << result.domain.normalizedResidualBound << ",\"provisional_selection\":";
    if (result.provisionalSelection) out << *result.provisionalSelection; else out << "null";
    out << ",\"tail_model_provenance\":{\"derivation\":\""
        << jsonEscape(result.problem.ciphertextModel.provenance.derivation)
        << "\",\"included_noise_sources\":\""
        << jsonEscape(result.problem.ciphertextModel.provenance.includedNoiseSources)
        << "\",\"assumptions\":\"" << jsonEscape(result.problem.ciphertextModel.provenance.assumptions)
        << "\"},\"candidates\":[";
    for (std::size_t i = 0; i < result.candidates.size(); ++i) {
        if (i) out << ',';
        const auto& c = result.candidates[i];
        out << "{\"family\":\"" << familyName(c.family) << "\",\"degree\":"
            << c.compiledCircuit.cost.degree << ",\"basis\":\"monomial\",\"evaluation_strategy\":"
            << "\"paterson_stockmeyer\",\"baby_step\":" << c.compiledCircuit.babyStep
            << ",\"depth\":" << c.compiledCircuit.cost.multiplicativeDepth << ",\"level_consumption\":"
            << c.compiledCircuit.cost.levelConsumption << ",\"real_error\":"
            << c.diagnostic.approximationMaxError << ",\"complex_error\":"
            << c.diagnostic.complexBoundaryErrorMax << ",\"complex_derivative\":"
            << c.diagnostic.complexDerivativeMax << ",\"predicted_error\":"
            << c.predictedBootstrapError << ",\"interval_error_bound\":\""
            << c.intervalCertificate.approximationErrorUpperBound << "\",\"certified\":"
            << (c.intervalCertified ? "true" : "false") << ",\"executable\":"
            << (c.executable ? "true" : "false") << ",\"stage\":\"" << stageName(c.stage)
            << "\",\"normalization_gain\":" << exactScaleUp(c.compiledCircuit.normalizationGain)
            << ",\"denormalization_gain\":" << exactScaleUp(c.compiledCircuit.denormalizationGain)
            << ",\"arithmetic_error_known\":"
            << (c.polynomialArithmeticError ? "true" : "false") << ",\"rejection_reason\":\""
            << rejectionName(c.rejectionReason) << "\",\"scale_schedule\":[";
        for (std::size_t stage = 0; stage < c.scaleSchedule.size(); ++stage) {
            if (stage) out << ',';
            out << "{\"input\":" << c.scaleSchedule[stage].inputScaleBits
                << ",\"growth\":" << c.scaleSchedule[stage].multiplicationGrowthBits
                << ",\"prime\":" << c.scaleSchedule[stage].rescalePrimeBits
                << ",\"output\":" << c.scaleSchedule[stage].outputScaleBits
                << ",\"available_modulus\":" << c.scaleSchedule[stage].availableModulusBits
                << ",\"headroom\":" << c.scaleSchedule[stage].requiredHeadroomBits << '}';
        }
        out << "]}";
    }
    out << "]}";
    return out.str();
}

} // namespace m2424::experimental
