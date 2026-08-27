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

struct RemezResult { EvalModPolynomial polynomial; bool converged{}; };

RemezResult remezOdd(const EvalModDomain& domain, std::size_t degree) {
    const std::size_t terms = (degree + 1) / 2;
    const std::size_t unknowns = terms + 1;
    std::vector<MpReal> grid;
    MpReal rho;
    if (mpfr_set_str(rho.get(), domain.normalizedResidualBoundDecimal.c_str(), 10, MPFR_RNDN) != 0)
        throw std::invalid_argument("invalid Remez domain radius");
    for (std::size_t integer = 0; integer <= domain.integerBound; ++integer) {
        for (std::size_t sample = 0; sample < 2048; ++sample) {
            MpReal residual, x, fraction;
            mpfr_set_ui(fraction.get(), 2 * sample, MPFR_RNDN);
            mpfr_div_ui(fraction.get(), fraction.get(), 2047, MPFR_RNDN);
            mpfr_sub_ui(fraction.get(), fraction.get(), 1, MPFR_RNDN);
            mpfr_mul(residual.get(), rho.get(), fraction.get(), MPFR_RNDN);
            mpfr_add_ui(x.get(), residual.get(), integer, MPFR_RNDN);
            if (mpfr_sgn(x.get()) >= 0) grid.push_back(std::move(x));
        }
    }
    std::vector<std::size_t> extrema(unknowns);
    for (std::size_t i = 0; i < unknowns; ++i) extrema[i] = i * (grid.size() - 1) / (unknowns - 1);
    std::vector<MpReal> solution(unknowns);
    bool converged = false;
    MpReal previousMaximum;
    mpfr_set_inf(previousMaximum.get(), 1);
    for (std::size_t iteration = 0; iteration < 24; ++iteration) {
        std::vector<std::vector<MpReal>> a;
        a.reserve(unknowns);
        for (std::size_t row = 0; row < unknowns; ++row) {
            a.emplace_back(unknowns + 1);
            MpReal base, residual, nearest;
            mpfr_set(base.get(), grid[extrema[row]].get(), MPFR_RNDN);
            mpfr_add_d(nearest.get(), base.get(), 0.5, MPFR_RNDN);
            mpfr_floor(nearest.get(), nearest.get());
            mpfr_sub(residual.get(), base.get(), nearest.get(), MPFR_RNDN);
            for (std::size_t col = 0; col < terms; ++col) {
                mpfr_pow_ui(a[row][col].get(), base.get(), 2 * col + 1, MPFR_RNDN);
            }
            mpfr_set_si(a[row][terms].get(), row % 2 ? -1 : 1, MPFR_RNDN);
            mpfr_set(a[row][unknowns].get(), residual.get(), MPFR_RNDN);
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
        const auto previousExtrema = extrema;
        std::vector<MpReal> extremaMagnitudes(unknowns);
        std::vector<int> extremaSigns(unknowns);
        MpReal globalMaximum;
        mpfr_set_zero(globalMaximum.get(), 0);
        int alternationBase = 0;
        for (std::size_t segment = 0; segment < unknowns; ++segment) {
            const std::size_t begin = segment * grid.size() / unknowns;
            const std::size_t end = (segment + 1) * grid.size() / unknowns;
            MpReal maximum;
            mpfr_set_si(maximum.get(), -1, MPFR_RNDN);
            MpReal fallbackMaximum;
            mpfr_set_si(fallbackMaximum.get(), -1, MPFR_RNDN);
            std::size_t fallbackPoint = begin;
            int fallbackSign = 0;
            for (std::size_t point = begin; point < end; ++point) {
                MpReal value, square, target, nearest, error;
                mpfr_mul(square.get(), grid[point].get(), grid[point].get(), MPFR_RNDN);
                mpfr_set_zero(value.get(), 0);
                for (std::size_t col = terms; col-- > 0;) {
                    mpfr_mul(value.get(), value.get(), square.get(), MPFR_RNDN);
                    mpfr_add(value.get(), value.get(), solution[col].get(), MPFR_RNDN);
                }
                mpfr_mul(value.get(), value.get(), grid[point].get(), MPFR_RNDN);
                mpfr_add_d(nearest.get(), grid[point].get(), 0.5, MPFR_RNDN);
                mpfr_floor(nearest.get(), nearest.get());
                mpfr_sub(target.get(), grid[point].get(), nearest.get(), MPFR_RNDN);
                mpfr_sub(error.get(), value.get(), target.get(), MPFR_RNDN);
                const int errorSign = mpfr_sgn(error.get());
                mpfr_abs(error.get(), error.get(), MPFR_RNDN);
                if (mpfr_greater_p(error.get(), fallbackMaximum.get())) {
                    mpfr_set(fallbackMaximum.get(), error.get(), MPFR_RNDN);
                    fallbackPoint = point; fallbackSign = errorSign;
                }
                if (alternationBase == 0 && errorSign != 0) alternationBase = errorSign;
                const int desiredSign = segment % 2 == 0 ? alternationBase : -alternationBase;
                if (errorSign == desiredSign && mpfr_greater_p(error.get(), maximum.get())) {
                    mpfr_set(maximum.get(), error.get(), MPFR_RNDN); extrema[segment] = point;
                    extremaSigns[segment] = errorSign;
                }
            }
            if (mpfr_sgn(maximum.get()) < 0) {
                mpfr_set(maximum.get(), fallbackMaximum.get(), MPFR_RNDN);
                extrema[segment] = fallbackPoint; extremaSigns[segment] = fallbackSign;
            }
            mpfr_set(extremaMagnitudes[segment].get(), maximum.get(), MPFR_RNDN);
            if (mpfr_greater_p(maximum.get(), globalMaximum.get()))
                mpfr_set(globalMaximum.get(), maximum.get(), MPFR_RNDN);
        }
        MpReal difference, tolerance;
        mpfr_sub(difference.get(), globalMaximum.get(), previousMaximum.get(), MPFR_RNDN);
        mpfr_abs(difference.get(), difference.get(), MPFR_RNDN);
        mpfr_mul_d(tolerance.get(), globalMaximum.get(), 5e-2, MPFR_RNDU);
        if (mpfr_cmp_d(tolerance.get(), 1e-12) < 0) mpfr_set_d(tolerance.get(), 1e-12, MPFR_RNDU);
        bool alternating = true;
        for (std::size_t i = 1; i < unknowns; ++i)
            alternating = alternating && extremaSigns[i] != 0
                && extremaSigns[i - 1] != 0 && extremaSigns[i] != extremaSigns[i - 1];
        MpReal minimumMagnitude, magnitudeRatio;
        mpfr_set(minimumMagnitude.get(), extremaMagnitudes.front().get(), MPFR_RNDN);
        for (std::size_t i = 1; i < unknowns; ++i)
            if (mpfr_less_p(extremaMagnitudes[i].get(), minimumMagnitude.get()))
                mpfr_set(minimumMagnitude.get(), extremaMagnitudes[i].get(), MPFR_RNDN);
        const bool balanced = mpfr_sgn(minimumMagnitude.get()) > 0
            && (mpfr_div(magnitudeRatio.get(), globalMaximum.get(), minimumMagnitude.get(), MPFR_RNDU),
                mpfr_cmp_d(magnitudeRatio.get(), 1.25) <= 0);
        const std::size_t exchangeTolerance = std::max<std::size_t>(1, grid.size() / 100);
        bool exchangeStable = true;
        for (std::size_t i = 0; i < unknowns; ++i) {
            const auto movement = extrema[i] > previousExtrema[i]
                ? extrema[i] - previousExtrema[i] : previousExtrema[i] - extrema[i];
            exchangeStable = exchangeStable && movement <= exchangeTolerance;
        }
        if (mpfr_number_p(previousMaximum.get()) && alternating && balanced && exchangeStable
            && mpfr_lessequal_p(difference.get(), tolerance.get())) { converged = true; break; }
        mpfr_set(previousMaximum.get(), globalMaximum.get(), MPFR_RNDN);
    }
    EvalModPolynomial result;
    result.basis = PolynomialBasis::Monomial;
    result.decimalCoefficients.assign(degree + 1, "0");
    for (std::size_t index = 0; index < terms; ++index) {
        char* text = nullptr;
        mpfr_asprintf(&text, "%.120Rg", solution[index].get());
        result.decimalCoefficients[2 * index + 1] = text ? text : "0";
        mpfr_free_str(text);
    }
    return {std::move(result), converged};
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

ExactInteger exactIntegerFromUint64(std::uint64_t value) {
    ExactInteger result;
    mpz_import(result.get_mpz_t(), 1, 1, sizeof(value), 0, 0, &value);
    return result;
}

ExactScale divideScale(const ExactScale& scale, std::uint64_t divisor) {
    if (divisor < 2) throw std::invalid_argument("invalid EvalMod rescale prime");
    return ExactScale::rational(scale.numerator,
                                scale.denominator * exactIntegerFromUint64(divisor));
}

bool equalScale(const ExactScale& left, const ExactScale& right) {
    return left.numerator == right.numerator && left.denominator == right.denominator;
}

mpq_class parseExactDecimal(const std::string& text) {
    if (text.empty()) throw std::invalid_argument("empty exact decimal");
    std::size_t position = 0;
    bool negative = false;
    if (text[position] == '+' || text[position] == '-') {
        negative = text[position] == '-';
        if (++position == text.size()) throw std::invalid_argument("invalid exact decimal");
    }
    const auto exponentPosition = text.find_first_of("eE", position);
    const auto mantissaEnd = exponentPosition == std::string::npos ? text.size() : exponentPosition;
    long exponent = 0;
    if (exponentPosition != std::string::npos) {
        std::size_t consumed = 0;
        try {
            exponent = std::stol(text.substr(exponentPosition + 1), &consumed, 10);
        } catch (const std::exception&) {
            throw std::invalid_argument("invalid exact decimal exponent");
        }
        if (consumed != text.size() - exponentPosition - 1
            || exponent < -100000 || exponent > 100000)
            throw std::invalid_argument("invalid exact decimal exponent");
    }
    std::string digits;
    std::size_t fractionalDigits = 0;
    bool decimalPoint = false;
    for (; position < mantissaEnd; ++position) {
        const char character = text[position];
        if (character == '.') {
            if (decimalPoint) throw std::invalid_argument("invalid exact decimal point");
            decimalPoint = true;
        } else if (character >= '0' && character <= '9') {
            digits.push_back(character);
            if (decimalPoint) ++fractionalDigits;
        } else {
            throw std::invalid_argument("invalid exact decimal digit");
        }
    }
    if (digits.empty()) throw std::invalid_argument("exact decimal has no digits");
    ExactInteger numerator(digits, 10);
    if (negative) numerator = -numerator;
    ExactInteger denominator = 1;
    if (fractionalDigits > 0)
        mpz_ui_pow_ui(denominator.get_mpz_t(), 10, static_cast<unsigned long>(fractionalDigits));
    if (exponent > 0) {
        ExactInteger power;
        mpz_ui_pow_ui(power.get_mpz_t(), 10, static_cast<unsigned long>(exponent));
        numerator *= power;
    } else if (exponent < 0) {
        ExactInteger power;
        mpz_ui_pow_ui(power.get_mpz_t(), 10, static_cast<unsigned long>(-exponent));
        denominator *= power;
    }
    mpq_class result(numerator, denominator);
    result.canonicalize();
    return result;
}

ExactInteger roundRationalAwayFromZero(const ExactInteger& numerator,
                                       const ExactInteger& denominator) {
    if (denominator <= 0) throw std::invalid_argument("invalid exact rounding denominator");
    const bool negative = numerator < 0;
    const ExactInteger magnitude = negative ? -numerator : numerator;
    ExactInteger quotient = magnitude / denominator;
    const ExactInteger remainder = magnitude % denominator;
    if (2 * remainder >= denominator) ++quotient;
    return negative ? -quotient : quotient;
}

std::pair<std::string, double> rationalUpperDecimal(const ExactInteger& numerator,
                                                    const ExactInteger& denominator) {
    mpfr_t value;
    mpfr_init2(value, 384);
    mpq_class rational(numerator, denominator);
    rational.canonicalize();
    mpfr_set_q(value, rational.get_mpq_t(), MPFR_RNDU);
    char* text = nullptr;
    mpfr_asprintf(&text, "%.120RUg", value);
    std::string decimal = text ? text : "0";
    mpfr_free_str(text);
    const double upper = mpfr_get_d(value, MPFR_RNDU);
    mpfr_clear(value);
    return {std::move(decimal), upper};
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

std::string exactScaleDecimal(const ExactScale& scale) {
    mpfr_t value; mpfr_init2(value, 384);
    mpq_class rational(scale.numerator, scale.denominator); rational.canonicalize();
    mpfr_set_q(value, rational.get_mpq_t(), MPFR_RNDN);
    char* text = nullptr; mpfr_asprintf(&text, "%.120Rg", value);
    std::string result = text ? text : "0"; mpfr_free_str(text); mpfr_clear(value);
    return result;
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
    case EvalModCandidateStage::BackendMeasured: return "backend_measured";
    case EvalModCandidateStage::BackendValidated: return "backend_validated";
    }
    return "unknown";
}

const char* rejectionName(EvalModRejectionReason reason) {
    switch (reason) {
    case EvalModRejectionReason::None: return "none";
    case EvalModRejectionReason::ApproximationError: return "approximation_error";
    case EvalModRejectionReason::ApproximationNotConverged: return "approximation_not_converged";
    case EvalModRejectionReason::Uncertified: return "uncertified";
    case EvalModRejectionReason::ArithmeticErrorUnknown: return "arithmetic_error_unknown";
    case EvalModRejectionReason::InsufficientLevels: return "insufficient_levels";
    case EvalModRejectionReason::ModulusBudget: return "modulus_budget";
    case EvalModRejectionReason::SecurityBudget: return "security_budget";
    case EvalModRejectionReason::ScaleScheduleFailure: return "scale_schedule_failure";
    case EvalModRejectionReason::HeadroomViolation: return "headroom_violation";
    }
    return "unknown";
}

std::vector<EvalModNodeErrorState> propagateNodeErrors(
    const CompiledEvalModCircuit& circuit, double rawInputBound, double inputErrorBound,
    std::size_t precisionBits, const EvalModArithmeticErrorModel& model) {
    const double unit = std::ldexp(1.0, -static_cast<int>(precisionBits));
    const double encodingError = model.encodingAbsolute > 0.0 ? model.encodingAbsolute : unit;
    const double additionError = model.additionAbsolute > 0.0 ? model.additionAbsolute : unit;
    const double multiplicationRelative = model.multiplicationRelative > 0.0
        ? model.multiplicationRelative : 2.0 * unit;
    std::vector<EvalModNodeErrorState> states(circuit.nodes.size());
    for (std::size_t index = 0; index < circuit.nodes.size(); ++index) {
        const auto& node = circuit.nodes[index];
        auto& state = states[index];
        if (node.operation == EvalModOperation::Input) {
            state = {rawInputBound, inputErrorBound, 0.0, inputErrorBound};
        } else if (node.operation == EvalModOperation::EncodeConstant) {
            const double value = std::abs(std::stod(node.constantDecimal));
            const double encoding = encodingError * std::max(1.0, value);
            state = {value, encoding, encoding / std::max(1.0, value), encoding};
        } else if (node.operation == EvalModOperation::Add
                   || node.operation == EvalModOperation::AddPlain) {
            const auto& left = states[node.inputs[0]], &right = states[node.inputs[1]];
            state.valueAbsBound = left.valueAbsBound + right.valueAbsBound;
            state.absoluteErrorBound = left.absoluteErrorBound + right.absoluteErrorBound + additionError;
            state.relativeScaleErrorBound = left.relativeScaleErrorBound
                + right.relativeScaleErrorBound;
            state.noiseBound = left.noiseBound + right.noiseBound + additionError;
        } else if (node.operation == EvalModOperation::MultiplyPlain
                   || node.operation == EvalModOperation::MultiplyCipher) {
            const auto& left = states[node.inputs[0]], &right = states[node.inputs[1]];
            state.valueAbsBound = left.valueAbsBound * right.valueAbsBound;
            state.absoluteErrorBound = left.valueAbsBound * right.absoluteErrorBound
                + right.valueAbsBound * left.absoluteErrorBound
                + left.absoluteErrorBound * right.absoluteErrorBound
                + multiplicationRelative * std::max(1.0, state.valueAbsBound);
            state.relativeScaleErrorBound = left.relativeScaleErrorBound
                + right.relativeScaleErrorBound + multiplicationRelative;
            state.noiseBound = left.valueAbsBound * right.noiseBound
                + right.valueAbsBound * left.noiseBound + multiplicationRelative;
        } else {
            state = states[node.inputs[0]];
            double absoluteConstant = unit;
            if (node.operation == EvalModOperation::Relinearize && model.relinearizationAbsolute > 0.0)
                absoluteConstant = model.relinearizationAbsolute;
            else if (node.operation == EvalModOperation::Rescale && model.rescaleAbsolute > 0.0)
                absoluteConstant = model.rescaleAbsolute;
            else if (node.operation == EvalModOperation::ModSwitch)
                absoluteConstant = model.modSwitchAbsolute;
            const double operationError = absoluteConstant * std::max(1.0, state.valueAbsBound);
            state.absoluteErrorBound += operationError;
            state.noiseBound += operationError;
            if (node.operation == EvalModOperation::AlignScale) {
                const double relative = model.metadataScaleRelative > 0.0
                    ? model.metadataScaleRelative
                    : std::exp2(circuit.maxMetadataScaleCorrectionLog2) - 1.0;
                state.relativeScaleErrorBound += relative;
                state.absoluteErrorBound += state.valueAbsBound * relative;
            } else if (node.operation == EvalModOperation::Rescale
                       || node.operation == EvalModOperation::Relinearize) {
                state.relativeScaleErrorBound += unit;
            }
        }
        if (!std::isfinite(state.valueAbsBound) || !std::isfinite(state.absoluteErrorBound)
            || !std::isfinite(state.relativeScaleErrorBound) || !std::isfinite(state.noiseBound))
            throw std::overflow_error("EvalMod node error propagation overflow");
    }
    return states;
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
                       std::string constant = {}, double correctionLog2 = 0.0,
                       double relativeScaleError = 0.0) {
        compiled.nodes.push_back({operation, std::move(inputs), chain, std::move(input),
                                  std::move(output), std::move(constant), correctionLog2,
                                  relativeScaleError});
        return compiled.nodes.size() - 1;
    };
    const std::size_t input = addNode(EvalModOperation::Input, {}, 0, workingScale, workingScale);
    std::size_t x = input;
    if (compiled.normalizationGain.numerator != compiled.normalizationGain.denominator) {
        const std::size_t normalization = addNode(EvalModOperation::EncodeConstant, {}, 0,
                                                  workingScale, workingScale,
                                                  exactScaleDecimal(compiled.normalizationGain));
        x = addNode(EvalModOperation::MultiplyPlain, {input, normalization}, 0,
                    workingScale, multiplyScale(workingScale, workingScale));
        x = addNode(EvalModOperation::Rescale, {x}, 1, compiled.nodes[x].outputScale, workingScale);
    }

    const std::size_t degree = polynomial.decimalCoefficients.size() - 1;
    std::vector<std::size_t> powers(babyStep + 1, x);
    std::vector<std::size_t> depths(compiled.nodes.size(), 0);
    depths.resize(compiled.nodes.size(), 0);
    auto alignPair = [&](std::size_t left, std::size_t right) {
        const std::size_t targetChain = std::max(compiled.nodes[left].chainIndex,
                                                 compiled.nodes[right].chainIndex);
        if (compiled.nodes[left].chainIndex < targetChain) {
            left = addNode(EvalModOperation::ModSwitch, {left, right}, targetChain,
                           compiled.nodes[left].outputScale, compiled.nodes[left].outputScale);
            depths.resize(compiled.nodes.size()); depths[left] = depths[compiled.nodes[left].inputs[0]];
        }
        if (compiled.nodes[right].chainIndex < targetChain) {
            right = addNode(EvalModOperation::ModSwitch, {right, left}, targetChain,
                            compiled.nodes[right].outputScale, compiled.nodes[right].outputScale);
            depths.resize(compiled.nodes.size()); depths[right] = depths[compiled.nodes[right].inputs[0]];
        }
        const double correction = std::log2(
            exactScaleUp(compiled.nodes[left].outputScale)
            / exactScaleUp(compiled.nodes[right].outputScale));
        right = addNode(EvalModOperation::AlignScale, {right, left}, targetChain,
                        compiled.nodes[right].outputScale, compiled.nodes[left].outputScale,
                        {}, correction, std::abs(std::exp2(correction) - 1.0));
        depths.resize(compiled.nodes.size()); depths[right] = depths[compiled.nodes[right].inputs[0]];
        return std::pair<std::size_t, std::size_t>{left, right};
    };
    const std::size_t maximumBabyPower = std::min(babyStep, degree);
    for (std::size_t power = 2; power <= maximumBabyPower; ++power) {
        const auto operands = alignPair(powers[power - 1], x);
        const auto multiplyChain = compiled.nodes[operands.first].chainIndex;
        auto multiplied = addNode(EvalModOperation::MultiplyCipher, {operands.first, operands.second}, multiplyChain,
                                  workingScale, multiplyScale(workingScale, workingScale));
        depths.resize(compiled.nodes.size()); depths[multiplied] = depths[powers[power - 1]] + 1;
        auto relin = addNode(EvalModOperation::Relinearize, {multiplied}, multiplyChain,
                             compiled.nodes[multiplied].outputScale, compiled.nodes[multiplied].outputScale);
        depths.resize(compiled.nodes.size()); depths[relin] = depths[multiplied];
        powers[power] = addNode(EvalModOperation::Rescale, {relin}, multiplyChain + 1, compiled.nodes[relin].outputScale,
                                workingScale);
        depths.resize(compiled.nodes.size()); depths[powers[power]] = depths[relin];
    }
    const std::size_t blocks = degree / babyStep + 1;
    std::vector<std::size_t> giantPowers{input};
    if (blocks > 1) giantPowers.push_back(powers[babyStep]);
    for (std::size_t block = 2; block < blocks; ++block) {
        const auto operands = alignPair(giantPowers.back(), powers[babyStep]);
        const auto multiplyChain = compiled.nodes[operands.first].chainIndex;
        auto mul = addNode(EvalModOperation::MultiplyCipher, {operands.first, operands.second}, multiplyChain,
                           workingScale, multiplyScale(workingScale, workingScale));
        depths.resize(compiled.nodes.size()); depths[mul] = depths[giantPowers.back()] + 1;
        auto relin = addNode(EvalModOperation::Relinearize, {mul}, multiplyChain,
                             compiled.nodes[mul].outputScale, compiled.nodes[mul].outputScale);
        depths.resize(compiled.nodes.size()); depths[relin] = depths[mul];
        auto rescale = addNode(EvalModOperation::Rescale, {relin}, multiplyChain + 1,
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
            const long double coefficientValue = std::stold(polynomial.decimalCoefficients[coefficient]);
            if (coefficientValue == 0.0L) continue;
            auto constant = addNode(EvalModOperation::EncodeConstant, {}, block, workingScale, workingScale,
                                    polynomial.decimalCoefficients[coefficient]);
            std::size_t term = constant;
            if (offset > 0) {
                if (coefficientValue == 1.0L) {
                    term = powers[offset];
                } else {
                    const auto termChain = compiled.nodes[powers[offset]].chainIndex;
                    term = addNode(EvalModOperation::MultiplyPlain, {powers[offset], constant}, termChain,
                                   workingScale, multiplyScale(workingScale, workingScale));
                    depths.resize(compiled.nodes.size()); depths[term] = depths[powers[offset]];
                    term = addNode(EvalModOperation::Rescale, {term}, termChain + 1,
                                   compiled.nodes[term].outputScale, workingScale);
                    depths.resize(compiled.nodes.size()); depths[term] = depths[powers[offset]];
                }
            }
            if (blockValue) {
                const auto previous = *blockValue;
                const auto current = term;
                const bool previousPlain = compiled.nodes[previous].operation == EvalModOperation::EncodeConstant;
                const bool currentPlain = compiled.nodes[current].operation == EvalModOperation::EncodeConstant;
                if (previousPlain != currentPlain) {
                    const auto cipher = previousPlain ? current : previous;
                    const auto plain = previousPlain ? previous : current;
                    term = addNode(EvalModOperation::AddPlain, {cipher, plain},
                                   compiled.nodes[cipher].chainIndex, workingScale, workingScale);
                } else {
                    const auto operands = alignPair(previous, current);
                    term = addNode(EvalModOperation::Add, {operands.first, operands.second},
                                   std::max(compiled.nodes[operands.first].chainIndex,
                                            compiled.nodes[operands.second].chainIndex),
                                   workingScale, workingScale);
                }
                depths.resize(compiled.nodes.size());
                depths[term] = std::max(depths[previous], depths[current]);
            } else depths.resize(compiled.nodes.size());
            blockValue = term;
        }
        if (!blockValue) continue;
        std::size_t term = *blockValue;
        if (block > 0) {
            if (compiled.nodes[term].operation == EvalModOperation::EncodeConstant) {
                const auto multiplyChain = compiled.nodes[giantPowers[block]].chainIndex;
                auto mul = addNode(EvalModOperation::MultiplyPlain,
                                   {giantPowers[block], term}, multiplyChain, workingScale,
                                   multiplyScale(workingScale, workingScale));
                depths.resize(compiled.nodes.size()); depths[mul] = depths[giantPowers[block]];
                term = addNode(EvalModOperation::Rescale, {mul}, multiplyChain + 1,
                               compiled.nodes[mul].outputScale, workingScale);
                depths.resize(compiled.nodes.size()); depths[term] = depths[mul];
            } else {
                const auto operands = alignPair(term, giantPowers[block]);
                const auto multiplyChain = compiled.nodes[operands.first].chainIndex;
                auto mul = addNode(EvalModOperation::MultiplyCipher,
                                   {operands.first, operands.second}, multiplyChain, workingScale,
                                   multiplyScale(workingScale, workingScale));
                depths.resize(compiled.nodes.size());
                depths[mul] = std::max(depths[term], depths[giantPowers[block]]) + 1;
                auto relin = addNode(EvalModOperation::Relinearize, {mul}, multiplyChain,
                                     compiled.nodes[mul].outputScale,
                                     compiled.nodes[mul].outputScale);
                depths.resize(compiled.nodes.size()); depths[relin] = depths[mul];
                term = addNode(EvalModOperation::Rescale, {relin}, multiplyChain + 1,
                               compiled.nodes[relin].outputScale, workingScale);
                depths.resize(compiled.nodes.size()); depths[term] = depths[relin];
            }
        }
        if (total) {
            const auto previous = *total;
            const auto current = term;
            const auto operands = alignPair(previous, current);
            term = addNode(EvalModOperation::Add, {operands.first, operands.second},
                           std::max(compiled.nodes[operands.first].chainIndex,
                                    compiled.nodes[operands.second].chainIndex),
                           workingScale, workingScale);
            depths.resize(compiled.nodes.size());
            depths[term] = std::max(depths[previous], depths[current]);
        } else depths.resize(compiled.nodes.size());
        total = term;
    }
    if (!total) throw std::invalid_argument("zero polynomial has no ciphertext circuit");
    auto output = *total;
    if (compiled.denormalizationGain.numerator != compiled.denormalizationGain.denominator) {
        const std::size_t outputChain = compiled.nodes[*total].chainIndex;
        const bool exactIntegerGain = compiled.denormalizationGain.denominator == 1;
        const ExactScale plaintextScale = exactIntegerGain
            ? ExactScale::rational(1, 1) : problem.outputScale;
        const auto denorm = addNode(EvalModOperation::EncodeConstant, {}, outputChain,
                                    plaintextScale, plaintextScale,
                                    exactScaleDecimal(compiled.denormalizationGain));
        output = addNode(EvalModOperation::MultiplyPlain, {*total, denorm}, outputChain,
                         workingScale, multiplyScale(workingScale, plaintextScale));
        depths.resize(compiled.nodes.size()); depths[output] = depths[*total];
        if (!exactIntegerGain) {
            output = addNode(EvalModOperation::Rescale, {output}, outputChain + 1,
                             compiled.nodes[output].outputScale, problem.outputScale);
            depths.resize(compiled.nodes.size()); depths[output] = depths[*total];
        }
    }
    compiled.outputNode = output;

    EvalModCircuitCost cost;
    cost.degree = degree;
    cost.multiplicativeDepth = *std::max_element(depths.begin(), depths.end());
    for (const auto& node : compiled.nodes) {
        cost.ciphertextMultiplications += node.operation == EvalModOperation::MultiplyCipher;
        cost.ciphertextPlaintextMultiplications += node.operation == EvalModOperation::MultiplyPlain;
        cost.additions += node.operation == EvalModOperation::Add
            || node.operation == EvalModOperation::AddPlain;
        cost.plaintextAdditions += node.operation == EvalModOperation::AddPlain;
        cost.modulusSwitches += node.operation == EvalModOperation::ModSwitch;
        cost.scaleAlignments += node.operation == EvalModOperation::AlignScale;
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

EvalModExactScaleSchedule buildExactEvalModScaleSchedule(
    const CompiledEvalModCircuit& circuit,
    const EvalModExactModulusContext& context,
    bool rigorous) {
    EvalModExactScaleSchedule result;
    result.available = true;
    result.nodeScales.resize(circuit.nodes.size());
    result.nodeModuli.resize(circuit.nodes.size());
    result.rescalePrimes.assign(circuit.nodes.size(), 0);

    const auto reject = [&](std::size_t node, std::string failure) {
        result.valid = false;
        result.rigorous = false;
        result.failingNode = node;
        result.failure = std::move(failure);
        return result;
    };
    if (context.dataPrimes.empty() || context.specialPrime < 2)
        return reject(0, "invalid_exact_modulus_context");
    for (std::size_t index = 0; index < context.dataPrimes.size(); ++index) {
        const auto prime = context.dataPrimes[index];
        if (prime < 2 || prime == context.specialPrime)
            return reject(0, "invalid_or_repeated_modulus_prime");
        for (std::size_t previous = 0; previous < index; ++previous)
            if (prime == context.dataPrimes[previous])
                return reject(0, "repeated_data_modulus_prime");
    }

    for (std::size_t index = 0; index < circuit.nodes.size(); ++index) {
        const auto& node = circuit.nodes[index];
        if (node.chainIndex >= context.dataPrimes.size())
            return reject(index, "node_chain_exceeds_data_modulus");
        for (const auto input : node.inputs)
            if (input >= index) return reject(index, "non_topological_evalmod_dag");

        ExactInteger modulus = 1;
        const std::size_t activePrimes = context.dataPrimes.size() - node.chainIndex;
        for (std::size_t prime = 0; prime < activePrimes; ++prime)
            modulus *= exactIntegerFromUint64(context.dataPrimes[prime]);
        result.nodeModuli[index] = std::move(modulus);

        if (node.operation == EvalModOperation::Input
            || node.operation == EvalModOperation::EncodeConstant) {
            result.nodeScales[index] = node.outputScale;
            continue;
        }
        if (node.inputs.empty()) return reject(index, "operation_has_no_input");
        const auto& left = result.nodeScales[node.inputs[0]];
        if (node.operation == EvalModOperation::Relinearize) {
            result.nodeScales[index] = left;
        } else if (node.operation == EvalModOperation::ModSwitch) {
            if (node.inputs.size() != 2
                || circuit.nodes[node.inputs[0]].chainIndex >= node.chainIndex
                || circuit.nodes[node.inputs[1]].chainIndex != node.chainIndex)
                return reject(index, "invalid_modswitch_target_level");
            result.nodeScales[index] = left;
        } else if (node.operation == EvalModOperation::Rescale) {
            if (node.chainIndex == 0)
                return reject(index, "rescale_does_not_advance_chain");
            const auto source = node.inputs[0];
            if (circuit.nodes[source].chainIndex + 1 != node.chainIndex)
                return reject(index, "rescale_chain_step_is_not_one");
            const std::size_t droppedIndex = context.dataPrimes.size() - node.chainIndex;
            const auto droppedPrime = context.dataPrimes[droppedIndex];
            result.rescalePrimes[index] = droppedPrime;
            result.nodeScales[index] = divideScale(left, droppedPrime);
        } else if (node.operation == EvalModOperation::AlignScale) {
            if (rigorous) return reject(index, "metadata_scale_alignment_prohibited");
            if (node.inputs.size() != 2)
                return reject(index, "metadata_scale_alignment_target_missing");
            result.nodeScales[index] = result.nodeScales[node.inputs[1]];
        } else {
            if (node.inputs.size() != 2)
                return reject(index, "binary_operation_input_count_mismatch");
            const auto& right = result.nodeScales[node.inputs[1]];
            if (node.operation == EvalModOperation::Add
                || node.operation == EvalModOperation::AddPlain) {
                if (!equalScale(left, right) && rigorous)
                    return reject(index, "addition_scale_mismatch");
                result.nodeScales[index] = left;
            } else if (node.operation == EvalModOperation::MultiplyPlain
                       || node.operation == EvalModOperation::MultiplyCipher) {
                result.nodeScales[index] = multiplyScale(left, right);
            } else {
                return reject(index, "unsupported_evalmod_operation");
            }
        }
    }
    result.valid = true;
    result.rigorous = rigorous;
    return result;
}

EvalModHeadroomCertificate certifyEvalModModSwitchHeadroom(
    const CompiledEvalModCircuit& circuit,
    const EvalModExactScaleSchedule& schedule,
    const std::vector<EvalModNodeErrorState>& nodeErrors,
    bool errorBoundsRigorous) {
    EvalModHeadroomCertificate result;
    result.available = schedule.available;
    if (!schedule.valid || schedule.nodeScales.size() != circuit.nodes.size()
        || schedule.nodeModuli.size() != circuit.nodes.size()
        || nodeErrors.size() != circuit.nodes.size()) {
        result.failure = "exact_scale_schedule_or_node_bounds_unavailable";
        return result;
    }
    for (std::size_t index = 0; index < circuit.nodes.size(); ++index) {
        if (circuit.nodes[index].operation != EvalModOperation::ModSwitch) continue;
        const double semanticMagnitude = nodeErrors[index].valueAbsBound
            + nodeErrors[index].absoluteErrorBound;
        if (!std::isfinite(semanticMagnitude) || semanticMagnitude < 0.0) {
            result.failingNode = index;
            result.failure = "invalid_modswitch_semantic_bound";
            return result;
        }
        ExactInteger encodedMagnitude = 0;
        if (semanticMagnitude > 0.0) {
            const auto semantic = ExactScale::fromBinaryDouble(semanticMagnitude);
            const auto& scale = schedule.nodeScales[index];
            const ExactInteger numerator = semantic.numerator * scale.numerator;
            const ExactInteger denominator = semantic.denominator * scale.denominator;
            encodedMagnitude = (numerator + denominator - 1) / denominator;
        }
        const auto& targetModulus = schedule.nodeModuli[index];
        const bool proved = 2 * encodedMagnitude < targetModulus;
        double headroomBits = std::numeric_limits<double>::infinity();
        if (encodedMagnitude > 0) {
            headroomBits = std::log2(mpz_get_d(targetModulus.get_mpz_t())) - 1.0
                - std::log2(mpz_get_d(encodedMagnitude.get_mpz_t()));
        }
        result.modSwitchGates.push_back({index, targetModulus, encodedMagnitude,
                                         headroomBits, proved,
                                         proved && schedule.rigorous && errorBoundsRigorous});
        if (!proved) {
            result.failingNode = index;
            result.failure = "modswitch_headroom_violation";
            return result;
        }
    }
    result.valid = true;
    result.rigorous = schedule.rigorous && errorBoundsRigorous;
    if (!result.rigorous) result.failure = "headroom_bounds_not_rigorous";
    return result;
}

PreparedEvalModConstants prepareEvalModConstants(
    SealAdapter& adapter,
    const CompiledEvalModCircuit& circuit,
    const Cipher& evalModInput) {
    if (!isCompiledEvalModCircuitValid(circuit))
        throw std::invalid_argument("invalid compiled EvalMod circuit");
    PreparedEvalModConstants result;
    result.inputDataPrimes = adapter.coeffModulusValues(evalModInput);
    if (result.inputDataPrimes.empty())
        throw std::invalid_argument("EvalMod input has no active data primes");
    const EvalModExactModulusContext exactContext{
        result.inputDataPrimes, adapter.specialKeyModulusValue()};
    const auto runtimeSchedule = buildExactEvalModScaleSchedule(circuit, exactContext, false);
    if (!runtimeSchedule.valid)
        throw std::invalid_argument("cannot prepare constants for invalid runtime scale schedule: "
                                    + runtimeSchedule.failure);
    std::vector<std::optional<std::size_t>> constantLevels(circuit.nodes.size());
    std::vector<std::optional<ExactScale>> constantScales(circuit.nodes.size());
    for (const auto& consumer : circuit.nodes) {
        for (const auto input : consumer.inputs) {
            if (circuit.nodes[input].operation != EvalModOperation::EncodeConstant) continue;
            if (constantLevels[input] && *constantLevels[input] != consumer.chainIndex)
                throw std::invalid_argument("EvalMod constant is consumed at multiple levels");
            constantLevels[input] = consumer.chainIndex;
            const ExactScale encodingScale = consumer.operation == EvalModOperation::AddPlain
                ? runtimeSchedule.nodeScales[consumer.inputs[0]]
                : circuit.nodes[input].outputScale;
            if (constantScales[input] && !equalScale(*constantScales[input], encodingScale))
                throw std::invalid_argument("EvalMod constant requires multiple encoding scales");
            constantScales[input] = encodingScale;
        }
    }
    for (std::size_t index = 0; index < circuit.nodes.size(); ++index) {
        const auto& node = circuit.nodes[index];
        if (node.operation != EvalModOperation::EncodeConstant) continue;
        if (!constantLevels[index])
            throw std::invalid_argument("unused EvalMod constant cannot be prepared");
        if (!constantScales[index])
            throw std::invalid_argument("EvalMod constant encoding scale is unavailable");
        const std::size_t targetLevel = *constantLevels[index];
        if (targetLevel >= result.inputDataPrimes.size())
            throw std::invalid_argument("EvalMod constant level exceeds active modulus chain");
        const auto coefficient = parseExactDecimal(node.constantDecimal);
        const auto encodingScale = ExactScale::fromBinaryDouble(
            exactScaleUp(*constantScales[index]));
        const ExactInteger scaledNumerator = coefficient.get_num() * encodingScale.numerator;
        const ExactInteger scaledDenominator = coefficient.get_den() * encodingScale.denominator;
        const ExactInteger rounded = roundRationalAwayFromZero(scaledNumerator, scaledDenominator);

        const std::size_t activePrimes = result.inputDataPrimes.size() - targetLevel;
        std::vector<std::uint64_t> residues;
        residues.reserve(activePrimes);
        for (std::size_t limb = 0; limb < activePrimes; ++limb) {
            const auto modulus = exactIntegerFromUint64(result.inputDataPrimes[limb]);
            ExactInteger residue = rounded % modulus;
            if (residue < 0) residue += modulus;
            residues.push_back(static_cast<std::uint64_t>(mpz_get_ui(residue.get_mpz_t())));
        }

        const ExactInteger errorNumerator = abs(
            rounded * encodingScale.denominator * coefficient.get_den()
            - coefficient.get_num() * encodingScale.numerator);
        const ExactInteger errorDenominator = encodingScale.numerator * coefficient.get_den();
        const auto error = rationalUpperDecimal(errorNumerator, errorDenominator);
        PreparedEvalModConstant prepared;
        prepared.node = index;
        prepared.decimal = node.constantDecimal;
        prepared.encodingScale = encodingScale;
        prepared.roundedScaledInteger = rounded;
        prepared.encodingErrorUpperBoundDecimal = error.first;
        prepared.encodingErrorUpperBound = error.second;
        prepared.plaintext = adapter.encodeScalarRnsAtScaleFor(
            residues, exactScaleUp(encodingScale), evalModInput, targetLevel);
        prepared.rigorous = true;
        result.constants.push_back(std::move(prepared));
    }
    result.rigorous = true;
    return result;
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
        if (spec.family == EvalModApproximationFamily::MultiIntervalMinimax) {
            try {
                auto remez = remezOdd(result.domain, spec.degree);
                candidate.polynomial = std::move(remez.polynomial);
                candidate.approximationConverged = remez.converged;
            } catch (const std::exception&) {
                candidate.polynomial.basis = PolynomialBasis::Monomial;
                candidate.polynomial.decimalCoefficients.assign(spec.degree + 1, "0");
                candidate.polynomial.decimalCoefficients[1] = "1";
                candidate.approximationConverged = false;
            }
        } else {
            candidate.polynomial = polynomial(spec.family == EvalModApproximationFamily::MinimaxInverseSine
                ? inverseSineCorrection(result.domain)
                : fit(result.domain, spec.degree,
                      spec.family == EvalModApproximationFamily::PeriodicSineBaseline));
            candidate.approximationConverged = true;
        }
        candidate.compiledCircuit = compileEvalModPolynomial(candidate.polynomial, problem, spec.babyStep);
        const std::size_t scaleBits = problem.targetPrecisionBits + 8;
        if (problem.exactModulusContext) {
            candidate.exactScaleSchedule = buildExactEvalModScaleSchedule(
                candidate.compiledCircuit, *problem.exactModulusContext,
                problem.requireCertifiedScaleSchedule);
            if (candidate.exactScaleSchedule.valid) {
                for (std::size_t index = 0; index < candidate.compiledCircuit.nodes.size(); ++index) {
                    if (candidate.compiledCircuit.nodes[index].operation != EvalModOperation::Rescale)
                        continue;
                    const auto source = candidate.compiledCircuit.nodes[index].inputs[0];
                    const auto prime = candidate.exactScaleSchedule.rescalePrimes[index];
                    const auto inputBits = exactScaleBitsUp(candidate.exactScaleSchedule.nodeScales[source]);
                    const auto outputBits = exactScaleBitsUp(candidate.exactScaleSchedule.nodeScales[index]);
                    const auto primeBits = static_cast<std::size_t>(
                        std::ceil(std::log2(static_cast<double>(prime))));
                    const auto modulusBits = static_cast<std::size_t>(mpz_sizeinbase(
                        candidate.exactScaleSchedule.nodeModuli[index].get_mpz_t(), 2));
                    candidate.scaleSchedule.push_back({inputBits,
                                                       inputBits > outputBits ? inputBits - outputBits : 0,
                                                       primeBits,
                                                       outputBits, modulusBits, 8});
                }
            }
        } else {
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
        const double normalization = exactScaleUp(candidate.compiledCircuit.normalizationGain);
        const double normalizedBound = static_cast<double>(result.domain.integerBound)
            + result.domain.normalizedResidualBound;
        candidate.nodeErrorStates = propagateNodeErrors(
            candidate.compiledCircuit, normalizedBound / normalization,
            problem.ciphertextModel.normalizedCoeffToSlotErrorAbsBound / normalization,
            problem.targetPrecisionBits, problem.arithmeticErrorModel);
        if (candidate.exactScaleSchedule.valid) {
            candidate.headroomCertificate = certifyEvalModModSwitchHeadroom(
                candidate.compiledCircuit, candidate.exactScaleSchedule,
                candidate.nodeErrorStates, problem.arithmeticErrorModel.rigorous);
        }
        candidate.polynomialArithmeticError =
            candidate.nodeErrorStates[candidate.compiledCircuit.outputNode].absoluteErrorBound;
        candidate.arithmeticErrorRigorous = problem.arithmeticErrorModel.rigorous;
        auto estimatedBounds = candidate.propagationBounds;
        estimatedBounds.approximation = candidate.diagnostic.approximationMaxError;
        estimatedBounds.polynomialArithmetic = *candidate.polynomialArithmeticError;
        candidate.estimatedBootstrapError = propagatedBootstrapError(estimatedBounds);
        candidate.propagationBounds.polynomialArithmetic = *candidate.polynomialArithmeticError;
        candidate.predictedBootstrapError = propagatedBootstrapError(candidate.propagationBounds);
        candidate.cost = estimateEvalModCost(candidate.compiledCircuit.cost, problem.backendCost, scaleBits);
        candidate.stage = EvalModCandidateStage::ScaleScheduled;
        candidate.satisfiesLevelBudget = candidate.compiledCircuit.cost.levelConsumption <= problem.availableLevels;
        candidate.intervalCertified = candidate.intervalCertificate.proved;
        candidate.circuitValid = isCompiledEvalModCircuitValid(candidate.compiledCircuit);
        candidate.backendRunnable = candidate.circuitValid && candidate.satisfiesLevelBudget;
        candidate.executable = candidate.backendRunnable;
        candidate.approximationCertified = candidate.intervalCertified;
        candidate.arithmeticErrorCertified = candidate.arithmeticErrorRigorous;
        candidate.analyticalArithmeticBound = candidate.polynomialArithmeticError;
        candidate.analyticalArithmeticBoundRigorous = candidate.arithmeticErrorRigorous;
        if (candidate.intervalCertified && candidate.arithmeticErrorRigorous)
            candidate.certifiedBootstrapError = candidate.predictedBootstrapError;
        candidate.satisfiesNumericalTarget = candidate.intervalCertified
            && candidate.arithmeticErrorRigorous && candidate.certifiedBootstrapError.has_value()
            && *candidate.certifiedBootstrapError <= problem.targetAbsoluteError;
        const bool scaleScheduleValid = std::all_of(candidate.scaleSchedule.begin(),
            candidate.scaleSchedule.end(), [](const EvalModScaleStage& stage) {
                return stage.availableModulusBits >= stage.outputScaleBits + stage.requiredHeadroomBits;
            });
        const bool exactScaleScheduleValid = !problem.requireCertifiedScaleSchedule
            || (candidate.exactScaleSchedule.available
                && candidate.exactScaleSchedule.valid
                && candidate.exactScaleSchedule.rigorous);
        const bool headroomValid = !problem.requireCertifiedScaleSchedule
            || (candidate.headroomCertificate.available
                && candidate.headroomCertificate.valid
                && candidate.headroomCertificate.rigorous);
        candidate.rejectionReason = !exactScaleScheduleValid || !scaleScheduleValid
            ? EvalModRejectionReason::ScaleScheduleFailure
            : !headroomValid ? EvalModRejectionReason::HeadroomViolation
            : !candidate.satisfiesLevelBudget
            ? EvalModRejectionReason::InsufficientLevels
            : candidate.cost.dataModulusBitsLowerBound > problem.maxDataModulusBits
                ? EvalModRejectionReason::ModulusBudget
            : !candidate.approximationConverged ? EvalModRejectionReason::ApproximationNotConverged
            : !candidate.intervalCertified ? EvalModRejectionReason::Uncertified
            : !candidate.arithmeticErrorRigorous ? EvalModRejectionReason::ArithmeticErrorUnknown
            : !candidate.satisfiesNumericalTarget ? EvalModRejectionReason::ApproximationError
            : EvalModRejectionReason::None;
        result.candidates.push_back(std::move(candidate));
    }
    for (std::size_t i = 0; i < result.candidates.size(); ++i) {
        const bool selectableFamily = result.candidates[i].family
            != EvalModApproximationFamily::MultiIntervalLeastSquaresPrototype;
        if (selectableFamily && result.candidates[i].approximationConverged
            && result.candidates[i].satisfiesLevelBudget
            && (!problem.requireCertifiedScaleSchedule
                || (result.candidates[i].exactScaleSchedule.valid
                    && result.candidates[i].exactScaleSchedule.rigorous
                    && result.candidates[i].headroomCertificate.valid
                    && result.candidates[i].headroomCertificate.rigorous))
            && result.candidates[i].estimatedBootstrapError <= problem.targetAbsoluteError
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

std::vector<double> evaluateEvalModReferenceMpfr(const EvalModPolynomial& polynomial,
                                                 const ExactScale& normalizationGain,
                                                 const ExactScale& denormalizationGain,
                                                 const std::vector<double>& rawInput,
                                                 std::size_t precisionBits) {
    if (polynomial.basis != PolynomialBasis::Monomial || precisionBits < 128)
        throw std::invalid_argument("invalid MPFR EvalMod reference input");
    const auto precision = static_cast<mpfr_prec_t>(precisionBits);
    mpq_class normalization(normalizationGain.numerator, normalizationGain.denominator);
    mpq_class denormalization(denormalizationGain.numerator, denormalizationGain.denominator);
    normalization.canonicalize(); denormalization.canonicalize();
    MpReal norm(precision), denorm(precision), x(precision), value(precision), coefficient(precision);
    mpfr_set_q(norm.get(), normalization.get_mpq_t(), MPFR_RNDN);
    mpfr_set_q(denorm.get(), denormalization.get_mpq_t(), MPFR_RNDN);
    std::vector<double> result;
    result.reserve(rawInput.size());
    for (double input : rawInput) {
        mpfr_set_d(x.get(), input, MPFR_RNDN);
        mpfr_mul(x.get(), x.get(), norm.get(), MPFR_RNDN);
        mpfr_set_zero(value.get(), 0);
        for (auto it = polynomial.decimalCoefficients.rbegin();
             it != polynomial.decimalCoefficients.rend(); ++it) {
            if (mpfr_set_str(coefficient.get(), it->c_str(), 10, MPFR_RNDN) != 0)
                throw std::invalid_argument("invalid MPFR polynomial coefficient");
            mpfr_mul(value.get(), value.get(), x.get(), MPFR_RNDN);
            mpfr_add(value.get(), value.get(), coefficient.get(), MPFR_RNDN);
        }
        mpfr_mul(value.get(), value.get(), denorm.get(), MPFR_RNDN);
        result.push_back(mpfr_get_d(value.get(), MPFR_RNDN));
    }
    return result;
}

std::vector<double> evaluateExactEvalModTargetMpfr(
    const ExactScale& normalizationGain, const ExactScale& denormalizationGain,
    const std::vector<double>& rawInput, std::size_t precisionBits) {
    if (precisionBits < 128) throw std::invalid_argument("invalid MPFR EvalMod target precision");
    const auto precision = static_cast<mpfr_prec_t>(precisionBits);
    mpq_class normalization(normalizationGain.numerator, normalizationGain.denominator);
    mpq_class denormalization(denormalizationGain.numerator, denormalizationGain.denominator);
    normalization.canonicalize(); denormalization.canonicalize();
    MpReal norm(precision), denorm(precision), x(precision), nearest(precision);
    mpfr_set_q(norm.get(), normalization.get_mpq_t(), MPFR_RNDN);
    mpfr_set_q(denorm.get(), denormalization.get_mpq_t(), MPFR_RNDN);
    std::vector<double> result;
    result.reserve(rawInput.size());
    for (double input : rawInput) {
        if (!std::isfinite(input)) throw std::invalid_argument("non-finite EvalMod target input");
        mpfr_set_d(x.get(), input, MPFR_RNDN);
        mpfr_mul(x.get(), x.get(), norm.get(), MPFR_RNDN);
        mpfr_add_d(nearest.get(), x.get(), 0.5, MPFR_RNDN);
        mpfr_floor(nearest.get(), nearest.get());
        mpfr_sub(x.get(), x.get(), nearest.get(), MPFR_RNDN);
        mpfr_mul(x.get(), x.get(), denorm.get(), MPFR_RNDN);
        result.push_back(mpfr_get_d(x.get(), MPFR_RNDN));
    }
    return result;
}

bool isCompiledEvalModCircuitValid(const CompiledEvalModCircuit& circuit) {
    if (circuit.nodes.empty() || circuit.outputNode >= circuit.nodes.size()
        || !std::isfinite(circuit.maxMetadataScaleCorrectionLog2)
        || circuit.maxMetadataScaleCorrectionLog2 < 0.0
        || !std::isfinite(circuit.maxPlannedScaleDriftLog2)
        || circuit.maxPlannedScaleDriftLog2 < circuit.maxMetadataScaleCorrectionLog2) return false;
    std::size_t inputs = 0;
    enum class ValueKind { Invalid, Scalar, Ciphertext };
    std::vector<ValueKind> kinds(circuit.nodes.size(), ValueKind::Invalid);
    for (std::size_t index = 0; index < circuit.nodes.size(); ++index) {
        const auto& node = circuit.nodes[index];
        for (auto source : node.inputs) if (source >= index) return false;
        const bool unary = node.operation == EvalModOperation::Relinearize
            || node.operation == EvalModOperation::Rescale;
        const bool binary = node.operation == EvalModOperation::Add
            || node.operation == EvalModOperation::AddPlain
            || node.operation == EvalModOperation::MultiplyPlain
            || node.operation == EvalModOperation::MultiplyCipher
            || node.operation == EvalModOperation::ModSwitch
            || node.operation == EvalModOperation::AlignScale;
        if ((node.operation == EvalModOperation::Input && (!node.inputs.empty() || ++inputs != 1))
            || (node.operation == EvalModOperation::EncodeConstant
                && (!node.inputs.empty() || node.constantDecimal.empty()))
            || (unary && node.inputs.size() != 1) || (binary && node.inputs.size() != 2)) return false;
        if (node.operation == EvalModOperation::Input) kinds[index] = ValueKind::Ciphertext;
        else if (node.operation == EvalModOperation::EncodeConstant) kinds[index] = ValueKind::Scalar;
        else if (unary) {
            if (kinds[node.inputs[0]] != ValueKind::Ciphertext) return false;
            kinds[index] = ValueKind::Ciphertext;
        } else if (node.operation == EvalModOperation::MultiplyPlain) {
            const auto left = kinds[node.inputs[0]], right = kinds[node.inputs[1]];
            if (!((left == ValueKind::Ciphertext && right == ValueKind::Scalar)
                  || (left == ValueKind::Scalar && right == ValueKind::Ciphertext))) return false;
            kinds[index] = ValueKind::Ciphertext;
        } else if (node.operation == EvalModOperation::AddPlain) {
            if (kinds[node.inputs[0]] != ValueKind::Ciphertext
                || kinds[node.inputs[1]] != ValueKind::Scalar) return false;
            kinds[index] = ValueKind::Ciphertext;
        } else {
            if (kinds[node.inputs[0]] != ValueKind::Ciphertext
                || kinds[node.inputs[1]] != ValueKind::Ciphertext) return false;
            kinds[index] = ValueKind::Ciphertext;
        }
    }
    return inputs == 1 && kinds[circuit.outputNode] == ValueKind::Ciphertext;
}

struct MpDagValue {
    explicit MpDagValue(mpfr_prec_t precision) : scalarValue(precision) {}
    bool scalar{};
    MpReal scalarValue;
    std::vector<MpReal> slots;
};

std::vector<MpDagValue> executeEvalModDagMpfrNodes(
    const CompiledEvalModCircuit& circuit, const std::vector<double>& input,
    mpfr_prec_t precision = 512) {
    std::vector<MpDagValue> values;
    values.reserve(circuit.nodes.size());
    for (std::size_t index = 0; index < circuit.nodes.size(); ++index) {
        values.emplace_back(precision);
        const auto& node = circuit.nodes[index];
        auto& output = values.back();
        if (node.operation == EvalModOperation::Input) {
            output.slots.reserve(input.size());
            for (double value : input) {
                output.slots.emplace_back(precision);
                mpfr_set_d(output.slots.back().get(), value, MPFR_RNDN);
            }
        } else if (node.operation == EvalModOperation::EncodeConstant) {
            output.scalar = true;
            if (mpfr_set_str(output.scalarValue.get(), node.constantDecimal.c_str(), 10, MPFR_RNDN) != 0)
                throw std::invalid_argument("invalid MPFR DAG constant");
        } else if (node.operation == EvalModOperation::Relinearize
                   || node.operation == EvalModOperation::Rescale
                   || node.operation == EvalModOperation::ModSwitch
                   || node.operation == EvalModOperation::AlignScale) {
            output.slots = values[node.inputs[0]].slots;
        } else {
            const auto& left = values[node.inputs[0]];
            const auto& right = values[node.inputs[1]];
            const auto& source = left.scalar ? right.slots : left.slots;
            output.slots.reserve(source.size());
            for (std::size_t slot = 0; slot < source.size(); ++slot) {
                output.slots.emplace_back(precision);
                mpfr_srcptr lhs = left.scalar ? left.scalarValue.get() : left.slots[slot].get();
                mpfr_srcptr rhs = right.scalar ? right.scalarValue.get() : right.slots[slot].get();
                if (node.operation == EvalModOperation::Add
                    || node.operation == EvalModOperation::AddPlain)
                    mpfr_add(output.slots.back().get(), lhs, rhs, MPFR_RNDN);
                else
                    mpfr_mul(output.slots.back().get(), lhs, rhs, MPFR_RNDN);
            }
        }
    }
    return values;
}

Cipher executeEvalModCircuit(SealAdapter& adapter, const CompiledEvalModCircuit& circuit,
                             const Cipher& coeffToSlotOutput, EvalModExecutionTrace* trace,
                             const std::vector<double>* semanticInput,
                             std::vector<EvalModNodeDifferential>* differentialTrace,
                             const PreparedEvalModConstants* preparedConstants) {
    if (!isCompiledEvalModCircuitValid(circuit))
        throw std::invalid_argument("invalid compiled EvalMod circuit");
    std::vector<const PreparedEvalModConstant*> preparedByNode(circuit.nodes.size(), nullptr);
    if (preparedConstants) {
        if (!preparedConstants->rigorous
            || preparedConstants->inputDataPrimes != adapter.coeffModulusValues(coeffToSlotOutput))
            throw std::invalid_argument("prepared EvalMod constants do not match input context and level");
        for (const auto& prepared : preparedConstants->constants) {
            if (!prepared.rigorous || prepared.node >= circuit.nodes.size()
                || circuit.nodes[prepared.node].operation != EvalModOperation::EncodeConstant
                || prepared.decimal != circuit.nodes[prepared.node].constantDecimal
                || preparedByNode[prepared.node])
                throw std::invalid_argument("invalid prepared EvalMod constant");
            preparedByNode[prepared.node] = &prepared;
        }
        for (std::size_t index = 0; index < circuit.nodes.size(); ++index)
            if (circuit.nodes[index].operation == EvalModOperation::EncodeConstant
                && !preparedByNode[index])
                throw std::invalid_argument("prepared EvalMod constant is missing");
    }
    struct Value { std::optional<Cipher> cipher; std::optional<double> scalar; };
    std::vector<Value> values(circuit.nodes.size());
    const std::size_t inputChainIndex = adapter.chainIndex(coeffToSlotOutput);
    if (trace) { trace->nodeStates.assign(circuit.nodes.size(), {}); trace->levelsConsumed = 0; }
    std::vector<MpDagValue> references;
    if (semanticInput) references = executeEvalModDagMpfrNodes(circuit, *semanticInput);
    if (differentialTrace) differentialTrace->clear();
    double previousCipherError = 0.0;
    for (std::size_t index = 0; index < circuit.nodes.size(); ++index) {
        const auto& node = circuit.nodes[index];
        if (node.operation == EvalModOperation::Input) values[index].cipher = coeffToSlotOutput;
        else if (node.operation == EvalModOperation::EncodeConstant) {
            if (!preparedConstants) values[index].scalar = std::stod(node.constantDecimal);
        }
        else if (node.operation == EvalModOperation::Relinearize)
            values[index].cipher = adapter.relinearize(*values[node.inputs[0]].cipher);
        else if (node.operation == EvalModOperation::Rescale)
            values[index].cipher = adapter.rescaleToNext(*values[node.inputs[0]].cipher);
        else if (node.operation == EvalModOperation::ModSwitch)
            values[index].cipher = adapter.modSwitchTo(*values[node.inputs[0]].cipher,
                                                       *values[node.inputs[1]].cipher);
        else if (node.operation == EvalModOperation::AlignScale) {
            const auto actual = adapter.scale(*values[node.inputs[0]].cipher);
            const auto target = adapter.scale(*values[node.inputs[1]].cipher);
            const double correction = std::abs(std::log2(actual / target));
            if (correction > circuit.maxMetadataScaleCorrectionLog2)
                throw std::runtime_error("branch scales require explicit arithmetic correction");
            values[index].cipher = adapter.normalizeScale(*values[node.inputs[0]].cipher, target);
        } else if (node.operation == EvalModOperation::MultiplyPlain) {
            const auto leftNode = node.inputs[0], rightNode = node.inputs[1];
            const bool leftConstant = circuit.nodes[leftNode].operation
                == EvalModOperation::EncodeConstant;
            const auto scalarNode = leftConstant ? leftNode : rightNode;
            const auto cipherNode = leftConstant ? rightNode : leftNode;
            const Cipher& cipher = *values[cipherNode].cipher;
            if (preparedConstants) {
                values[index].cipher = adapter.multiplyPlain(
                    cipher, preparedByNode[scalarNode]->plaintext);
            } else {
                values[index].cipher = adapter.multiplyPlain(
                    cipher, adapter.encodeScalarAtScaleFor(
                        *values[scalarNode].scalar,
                        exactScaleUp(circuit.nodes[scalarNode].outputScale), cipher));
            }
        } else if (node.operation == EvalModOperation::MultiplyCipher)
            values[index].cipher = adapter.multiply(*values[node.inputs[0]].cipher,
                                                    *values[node.inputs[1]].cipher);
        else if (node.operation == EvalModOperation::Add)
            values[index].cipher = adapter.add(*values[node.inputs[0]].cipher,
                                               *values[node.inputs[1]].cipher);
        else if (node.operation == EvalModOperation::AddPlain) {
            const Cipher& cipher = *values[node.inputs[0]].cipher;
            if (preparedConstants) {
                values[index].cipher = adapter.addPlain(
                    cipher, preparedByNode[node.inputs[1]]->plaintext);
            } else {
                values[index].cipher = adapter.addPlain(
                    cipher, adapter.encodeScalarFor(*values[node.inputs[1]].scalar, cipher));
            }
        }
        if (values[index].cipher) {
            const auto info = adapter.info(*values[index].cipher);
            if (inputChainIndex < info.chainIndex
                || inputChainIndex - info.chainIndex != node.chainIndex)
                throw std::runtime_error("planned chain index mismatch at node " + std::to_string(index));
            const double plannedScale = exactScaleUp(node.outputScale);
            if (std::abs(std::log2(info.scale / plannedScale))
                > circuit.maxPlannedScaleDriftLog2)
                throw std::runtime_error("planned scale mismatch at node " + std::to_string(index));
            if (trace) trace->nodeStates[index] = info;
            if (semanticInput && differentialTrace) {
                const auto decoded = adapter.decode(adapter.decrypt(*values[index].cipher));
                double maximumError = 0.0, referenceAtMaximum = 0.0, actualAtMaximum = 0.0;
                for (std::size_t slot = 0; slot < semanticInput->size(); ++slot) {
                    const double reference = mpfr_get_d(references[index].slots[slot].get(), MPFR_RNDN);
                    const double error = std::abs(decoded[slot] - reference);
                    if (error >= maximumError) {
                        maximumError = error; referenceAtMaximum = reference;
                        actualAtMaximum = decoded[slot];
                    }
                }
                differentialTrace->push_back({index, node.operation, info.chainIndex, info.scale,
                    referenceAtMaximum, actualAtMaximum, maximumError,
                    std::max(0.0, maximumError - previousCipherError)});
                previousCipherError = maximumError;
            }
        }
    }
    Cipher output = *values[circuit.outputNode].cipher;
    if (trace) {
        trace->outputScale = adapter.scale(output);
        trace->levelsConsumed = inputChainIndex - adapter.chainIndex(output);
    }
    return output;
}

EvalModCoeffToSlotResult executeEvalModAfterCoeffToSlot(
    SealAdapter& adapter, const CompiledEvalModCircuit& circuit,
    CoeffToSlotResult coeffToSlotOutput) {
    EvalModCoeffToSlotResult result;
    result.slotCipherFirst = executeEvalModCircuit(
        adapter, circuit, coeffToSlotOutput.slotCipherFirst, &result.firstTrace);
    result.slotCipherSecond = executeEvalModCircuit(
        adapter, circuit, coeffToSlotOutput.slotCipherSecond, &result.secondTrace);
    return result;
}

std::vector<double> executeEvalModDagPlaintext(const CompiledEvalModCircuit& circuit,
                                               const std::vector<double>& input) {
    struct Value { std::vector<double> slots; std::optional<double> scalar; };
    std::vector<Value> values(circuit.nodes.size());
    for (std::size_t index = 0; index < circuit.nodes.size(); ++index) {
        const auto& node = circuit.nodes[index];
        if (node.operation == EvalModOperation::Input) values[index].slots = input;
        else if (node.operation == EvalModOperation::EncodeConstant)
            values[index].scalar = std::stod(node.constantDecimal);
        else if (node.operation == EvalModOperation::Relinearize
                 || node.operation == EvalModOperation::Rescale
                 || node.operation == EvalModOperation::ModSwitch
                 || node.operation == EvalModOperation::AlignScale) values[index] = values[node.inputs[0]];
        else {
            const auto apply = [&](double left, double right) {
                return node.operation == EvalModOperation::Add
                    || node.operation == EvalModOperation::AddPlain ? left + right : left * right;
            };
            const auto& left = values[node.inputs[0]];
            const auto& right = values[node.inputs[1]];
            if (left.scalar && right.scalar) values[index].scalar = apply(*left.scalar, *right.scalar);
            else {
                const auto& source = left.scalar ? right.slots : left.slots;
                values[index].slots.resize(source.size());
                for (std::size_t slot = 0; slot < source.size(); ++slot) {
                    const double a = left.scalar ? *left.scalar : left.slots[slot];
                    const double b = right.scalar ? *right.scalar : right.slots[slot];
                    values[index].slots[slot] = apply(a, b);
                }
            }
        }
    }
    return values[circuit.outputNode].slots;
}

EvalModBackendValidation validateEvalModCandidateBackend(EvalModCandidate& candidate,
                                                        const EvalModProblem& problem,
                                                        const CkksProfile& profile,
                                                        const std::vector<double>& rawInput) {
    EvalModBackendValidation result;
    candidate.executable = false;
    candidate.circuitValid = isCompiledEvalModCircuitValid(candidate.compiledCircuit);
    candidate.backendRunnable = false;
    candidate.backendMeasured = false;
    if (rawInput.empty() || !candidate.circuitValid
        || candidate.family == EvalModApproximationFamily::MultiIntervalLeastSquaresPrototype
        || !candidate.satisfiesLevelBudget || !candidate.approximationConverged) {
        result.failure = "candidate_not_backend_eligible";
        return result;
    }
    if (profile.coeffModulusBits.size() <= candidate.compiledCircuit.cost.levelConsumption) {
        result.failure = "profile_chain_too_short";
        return result;
    }
    candidate.backendRunnable = true;
    try {
        SealAdapter adapter = SealAdapter::create(profile);
        adapter.generateKeys(std::vector<int>{}, true);
        const Cipher encrypted = adapter.encrypt(adapter.encode(rawInput));
        const auto preparedConstants = prepareEvalModConstants(
            adapter, candidate.compiledCircuit, encrypted);
        result.preparedConstantsUsed = preparedConstants.rigorous;
        for (const auto& prepared : preparedConstants.constants)
            result.maxPreparedConstantEncodingError = std::max(
                result.maxPreparedConstantEncodingError,
                prepared.encodingErrorUpperBound);
        EvalModExecutionTrace trace;
        const Cipher output = executeEvalModCircuit(adapter, candidate.compiledCircuit, encrypted,
                                                    &trace, &rawInput, &result.differentialTrace,
                                                    &preparedConstants);
        const auto decoded = adapter.decode(adapter.decrypt(output));
        const auto polynomialReference = evaluateEvalModReferenceMpfr(
            candidate.polynomial, candidate.compiledCircuit.normalizationGain,
            candidate.compiledCircuit.denormalizationGain, rawInput, 512);
        const auto evalModTarget = evaluateExactEvalModTargetMpfr(
            candidate.compiledCircuit.normalizationGain,
            candidate.compiledCircuit.denormalizationGain, rawInput, 512);
        for (std::size_t index = 0; index < rawInput.size(); ++index) {
            result.implementationError = std::max(result.implementationError,
                std::abs(decoded[index] - polynomialReference[index]));
            result.approximationError = std::max(result.approximationError,
                std::abs(polynomialReference[index] - evalModTarget[index]));
            result.totalMeasuredError = std::max(result.totalMeasuredError,
                std::abs(decoded[index] - evalModTarget[index]));
        }
        result.maxAbsoluteError = result.implementationError;
        result.executedNodes = candidate.compiledCircuit.nodes.size();
        result.outputChainIndex = adapter.chainIndex(output);
        result.outputScale = adapter.scale(output);
        result.executionSucceeded = true;
        const double implementationBudget = problem.precisionBudget.implementation > 0.0
            ? problem.precisionBudget.implementation : problem.targetAbsoluteError;
        const double approximationBudget = problem.precisionBudget.approximation > 0.0
            ? problem.precisionBudget.approximation : problem.targetAbsoluteError;
        result.matchesPolynomialReference = result.implementationError <= implementationBudget;
        const bool approximationWithinBudget = result.approximationError <= approximationBudget;
        result.matchesEvalModTarget = result.totalMeasuredError <= problem.targetAbsoluteError;
        result.predictionCoveredMeasurement = result.totalMeasuredError <= candidate.predictedBootstrapError;
        for (const auto& node : result.differentialTrace) {
            if (node.absoluteError > implementationBudget) {
                result.firstImplementationBudgetExceedingNode = node.node;
                break;
            }
        }
        result.passed = result.executionSucceeded && result.matchesPolynomialReference
            && approximationWithinBudget && result.matchesEvalModTarget;
        if (!result.matchesPolynomialReference) result.failure = "implementation_budget_exceeded";
        else if (!approximationWithinBudget) result.failure = "approximation_budget_exceeded";
        else if (!result.matchesEvalModTarget) result.failure = "total_budget_exceeded";
        candidate.measuredBackendError = result.implementationError;
        candidate.backendMeasured = result.executionSucceeded;
        candidate.approximationCertified = candidate.intervalCertified;
        candidate.arithmeticErrorCertified = candidate.arithmeticErrorRigorous;
        candidate.rigorouslyValidated = candidate.backendMeasured && candidate.approximationCertified
            && candidate.arithmeticErrorCertified && candidate.satisfiesNumericalTarget;
        candidate.executable = candidate.backendRunnable;
        if (candidate.backendMeasured) candidate.stage = EvalModCandidateStage::BackendMeasured;
        if (candidate.rigorouslyValidated) {
            candidate.stage = EvalModCandidateStage::BackendValidated;
            candidate.rejectionReason = EvalModRejectionReason::None;
        }
    } catch (const std::exception& error) {
        result.failure = error.what();
        candidate.backendRunnable = false;
        candidate.executable = false;
    }
    return result;
}

EvalModArithmeticErrorModel calibrateEvalModArithmeticModel(
    const std::vector<EvalModBackendValidation>& measurements, double safetyFactor) {
    if (measurements.empty() || !std::isfinite(safetyFactor) || safetyFactor < 1.0)
        throw std::invalid_argument("invalid EvalMod arithmetic calibration input");
    double maximum = 0.0;
    for (const auto& measurement : measurements) {
        if (!measurement.executionSucceeded || !std::isfinite(measurement.implementationError)
            || measurement.implementationError < 0.0)
            throw std::invalid_argument("calibration requires successful finite measurements");
        maximum = std::max(maximum, measurement.implementationError);
    }
    const double conservative = std::max(std::numeric_limits<double>::epsilon(),
                                         maximum * safetyFactor);
    return {conservative, conservative, conservative, conservative, conservative,
            conservative, conservative, true, false,
            "maximum measured implementation error multiplied by safety factor "
                + std::to_string(safetyFactor)};
}

std::string evalModSynthesisCsv(const EvalModSynthesisResult& result) {
    std::ostringstream out;
    out << "family,K,rho,degree,basis,strategy,baby_step,depth,ct_ct,ct_pt,rescales,level_consumption,"
           "mod_switches,scale_alignments,plaintext_additions,modulus_bits,normalization_gain,denormalization_gain,real_error,complex_error,"
           "complex_derivative,interval_error,approximation_converged,arithmetic_error_known,arithmetic_error_rigorous,estimated_error,certified_error,measured_backend_error,levels_ok,"
           "certified,circuit_valid,backend_runnable,backend_measured,rigorously_validated,executable,"
           "exact_scale_available,exact_scale_valid,exact_scale_rigorous,headroom_valid,headroom_rigorous,"
           "stage,rejection_reason,provisional\n";
    for (std::size_t i = 0; i < result.candidates.size(); ++i) {
        const auto& c = result.candidates[i];
        out << familyName(c.family) << ',' << c.domain.integerBound << ',' << std::setprecision(17)
            << c.domain.normalizedResidualBound << ',' << c.compiledCircuit.cost.degree
            << ",monomial,paterson_stockmeyer," << c.compiledCircuit.babyStep << ','
            << c.compiledCircuit.cost.multiplicativeDepth << ','
            << c.compiledCircuit.cost.ciphertextMultiplications << ','
            << c.compiledCircuit.cost.ciphertextPlaintextMultiplications << ','
            << c.compiledCircuit.cost.rescales << ',' << c.compiledCircuit.cost.levelConsumption << ','
            << c.compiledCircuit.cost.modulusSwitches << ','
            << c.compiledCircuit.cost.scaleAlignments << ','
            << c.compiledCircuit.cost.plaintextAdditions << ','
            << c.cost.dataModulusBitsLowerBound << ',' << exactScaleUp(c.compiledCircuit.normalizationGain)
            << ',' << exactScaleUp(c.compiledCircuit.denormalizationGain) << ','
            << c.diagnostic.approximationMaxError << ','
            << c.diagnostic.complexBoundaryErrorMax << ',' << c.diagnostic.realDerivativeMax << ','
            << c.intervalCertificate.approximationErrorUpperBoundDouble << ','
            << c.approximationConverged << ',' << c.polynomialArithmeticError.has_value() << ','
            << c.arithmeticErrorRigorous << ',' << c.estimatedBootstrapError << ',';
        if (c.certifiedBootstrapError) out << *c.certifiedBootstrapError;
        out << ',';
        if (c.measuredBackendError) out << *c.measuredBackendError;
        out << ','
            << c.satisfiesLevelBudget << ',' << c.intervalCertified << ',' << c.circuitValid << ','
            << c.backendRunnable << ',' << c.backendMeasured << ',' << c.rigorouslyValidated << ','
            << c.executable << ',' << c.exactScaleSchedule.available << ','
            << c.exactScaleSchedule.valid << ',' << c.exactScaleSchedule.rigorous << ','
            << c.headroomCertificate.valid << ',' << c.headroomCertificate.rigorous << ','
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
            << c.compiledCircuit.cost.levelConsumption << ",\"mod_switches\":"
            << c.compiledCircuit.cost.modulusSwitches << ",\"scale_alignments\":"
            << c.compiledCircuit.cost.scaleAlignments << ",\"plaintext_additions\":"
            << c.compiledCircuit.cost.plaintextAdditions << ",\"real_error\":"
            << c.diagnostic.approximationMaxError << ",\"complex_error\":"
            << c.diagnostic.complexBoundaryErrorMax << ",\"complex_derivative\":"
            << c.diagnostic.complexDerivativeMax << ",\"estimated_error\":"
            << c.estimatedBootstrapError << ",\"certified_error\":";
        if (c.certifiedBootstrapError) out << *c.certifiedBootstrapError; else out << "null";
        out << ",\"interval_error_bound\":\""
            << c.intervalCertificate.approximationErrorUpperBound << "\",\"certified\":"
            << (c.intervalCertified ? "true" : "false") << ",\"executable\":"
            << (c.executable ? "true" : "false") << ",\"stage\":\"" << stageName(c.stage)
            << "\",\"circuit_valid\":" << (c.circuitValid ? "true" : "false")
            << ",\"backend_runnable\":" << (c.backendRunnable ? "true" : "false")
            << ",\"backend_measured\":" << (c.backendMeasured ? "true" : "false")
            << ",\"rigorously_validated\":" << (c.rigorouslyValidated ? "true" : "false")
            << ",\"approximation_converged\":"
            << (c.approximationConverged ? "true" : "false")
            << ",\"arithmetic_error_rigorous\":"
            << (c.arithmeticErrorRigorous ? "true" : "false")
            << ",\"exact_scale_schedule\":{\"available\":"
            << (c.exactScaleSchedule.available ? "true" : "false")
            << ",\"valid\":" << (c.exactScaleSchedule.valid ? "true" : "false")
            << ",\"rigorous\":" << (c.exactScaleSchedule.rigorous ? "true" : "false")
            << ",\"failure\":\"" << jsonEscape(c.exactScaleSchedule.failure) << "\"}"
            << ",\"headroom_certificate\":{\"available\":"
            << (c.headroomCertificate.available ? "true" : "false")
            << ",\"valid\":" << (c.headroomCertificate.valid ? "true" : "false")
            << ",\"rigorous\":" << (c.headroomCertificate.rigorous ? "true" : "false")
            << ",\"mod_switch_gates\":" << c.headroomCertificate.modSwitchGates.size()
            << ",\"failure\":\"" << jsonEscape(c.headroomCertificate.failure) << "\"}"
            << ",\"measured_backend_error\":";
        if (c.measuredBackendError) out << *c.measuredBackendError; else out << "null";
        out
            << ",\"normalization_gain\":" << exactScaleUp(c.compiledCircuit.normalizationGain)
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
