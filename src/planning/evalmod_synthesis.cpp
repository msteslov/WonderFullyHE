#include "m2424/experimental/evalmod_analysis/synthesis.hpp"

#include <mpfr.h>
#include <seal/util/config.h>

#include <algorithm>
#include <cmath>
#include <cstring>
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

void chebyshevValue(mpfr_ptr output, mpfr_srcptr input, std::size_t degree) {
    if (degree == 0) { mpfr_set_ui(output, 1, MPFR_RNDN); return; }
    if (degree == 1) { mpfr_set(output, input, MPFR_RNDN); return; }
    MpReal previous(mpfr_get_prec(input)), current(mpfr_get_prec(input));
    MpReal next(mpfr_get_prec(input)), product(mpfr_get_prec(input));
    mpfr_set_ui(previous.get(), 1, MPFR_RNDN);
    mpfr_set(current.get(), input, MPFR_RNDN);
    for (std::size_t index = 2; index <= degree; ++index) {
        mpfr_mul(product.get(), input, current.get(), MPFR_RNDN);
        mpfr_mul_ui(product.get(), product.get(), 2, MPFR_RNDN);
        mpfr_sub(next.get(), product.get(), previous.get(), MPFR_RNDN);
        mpfr_set(previous.get(), current.get(), MPFR_RNDN);
        mpfr_set(current.get(), next.get(), MPFR_RNDN);
    }
    mpfr_set(output, current.get(), MPFR_RNDN);
}

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

RemezResult remezOdd(const EvalModDomain& domain, std::size_t degree,
                     bool chebyshevBasis = false) {
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
                if (chebyshevBasis)
                    chebyshevValue(a[row][col].get(), base.get(), 2 * col + 1);
                else
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
                mpfr_set_zero(value.get(), 0);
                if (chebyshevBasis) {
                    for (std::size_t col = 0; col < terms; ++col) {
                        MpReal basisValue, contribution;
                        chebyshevValue(basisValue.get(), grid[point].get(), 2 * col + 1);
                        mpfr_mul(contribution.get(), solution[col].get(),
                                 basisValue.get(), MPFR_RNDN);
                        mpfr_add(value.get(), value.get(), contribution.get(), MPFR_RNDN);
                    }
                } else {
                    mpfr_mul(square.get(), grid[point].get(), grid[point].get(), MPFR_RNDN);
                    for (std::size_t col = terms; col-- > 0;) {
                        mpfr_mul(value.get(), value.get(), square.get(), MPFR_RNDN);
                        mpfr_add(value.get(), value.get(), solution[col].get(), MPFR_RNDN);
                    }
                    mpfr_mul(value.get(), value.get(), grid[point].get(), MPFR_RNDN);
                }
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
    result.basis = chebyshevBasis ? PolynomialBasis::Chebyshev
                                  : PolynomialBasis::Monomial;
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

double exactScaleUp(const ExactScale& scale);

std::uint64_t binary64Bits(double value) {
    std::uint64_t bits{};
    static_assert(sizeof(bits) == sizeof(value), "unexpected binary64 size");
    std::memcpy(&bits, &value, sizeof(bits));
    return bits;
}

double binary64Value(std::uint64_t bits) {
    double value{};
    std::memcpy(&value, &bits, sizeof(value));
    return value;
}

EvalModScaleValue scaleValue(const ExactScale& ideal, double runtime) {
    if (!std::isfinite(runtime) || runtime <= 0.0)
        throw std::overflow_error("non-finite EvalMod runtime scale");
    return {ideal, binary64Bits(runtime), ExactScale::fromBinaryDouble(runtime)};
}

EvalModScaleValue initialScaleValue(const ExactScale& scale) {
    return scaleValue(scale, exactScaleUp(scale));
}

EvalModScaleValue multiplyScaleValue(const EvalModScaleValue& left,
                                     const EvalModScaleValue& right) {
    volatile double runtime = binary64Value(left.runtimeBinary64Bits)
        * binary64Value(right.runtimeBinary64Bits);
    return scaleValue(multiplyScale(left.ideal, right.ideal), runtime);
}

EvalModScaleValue divideScaleValue(const EvalModScaleValue& input, std::uint64_t prime) {
    volatile double runtime = binary64Value(input.runtimeBinary64Bits)
        / static_cast<double>(prime);
    return scaleValue(divideScale(input.ideal, prime), runtime);
}

std::optional<double> identityCorrectionScale(double branchScale, double targetScale,
                                              std::uint64_t rescalePrime) {
    if (!std::isfinite(branchScale) || branchScale <= 0.0
        || !std::isfinite(targetScale) || targetScale <= 0.0 || rescalePrime < 2)
        return std::nullopt;
    volatile double estimate = targetScale * static_cast<double>(rescalePrime) / branchScale;
    if (!std::isfinite(estimate) || estimate <= 0.0) return std::nullopt;
    const auto matches = [&](double correction) {
        volatile double multiplied = branchScale * correction;
        volatile double rescaled = multiplied / static_cast<double>(rescalePrime);
        return binary64Bits(rescaled) == binary64Bits(targetScale);
    };
    if (matches(estimate)) return estimate;
    double lower = estimate, upper = estimate;
    for (std::size_t step = 0; step < 4096; ++step) {
        lower = std::nextafter(lower, 0.0);
        if (matches(lower)) return lower;
        upper = std::nextafter(upper, std::numeric_limits<double>::infinity());
        if (matches(upper)) return upper;
    }
    return std::nullopt;
}

std::optional<double> identityMultiplyScale(double branchScale, double targetScale) {
    if (!std::isfinite(branchScale) || branchScale <= 0.0
        || !std::isfinite(targetScale) || targetScale <= 0.0)
        return std::nullopt;
    volatile double estimate = targetScale / branchScale;
    const auto matches = [&](double correction) {
        volatile double multiplied = branchScale * correction;
        return binary64Bits(multiplied) == binary64Bits(targetScale);
    };
    if (std::isfinite(estimate) && estimate > 0.0 && matches(estimate)) return estimate;
    double lower = estimate, upper = estimate;
    for (std::size_t step = 0; step < 4096; ++step) {
        lower = std::nextafter(lower, 0.0);
        if (matches(lower)) return lower;
        upper = std::nextafter(upper, std::numeric_limits<double>::infinity());
        if (matches(upper)) return upper;
    }
    return std::nullopt;
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

std::string exactTerminatingDecimal(const mpq_class& input) {
    mpq_class value = input;
    value.canonicalize();
    ExactInteger denominator = value.get_den();
    std::size_t twos = 0, fives = 0;
    while (mpz_divisible_ui_p(denominator.get_mpz_t(), 2)) {
        denominator /= 2; ++twos;
    }
    while (mpz_divisible_ui_p(denominator.get_mpz_t(), 5)) {
        denominator /= 5; ++fives;
    }
    if (denominator != 1)
        throw std::invalid_argument("Chebyshev decimal conversion is not terminating");
    const std::size_t digits = std::max(twos, fives);
    ExactInteger scaled = value.get_num();
    if (digits > twos) {
        ExactInteger power;
        mpz_ui_pow_ui(power.get_mpz_t(), 2, static_cast<unsigned long>(digits - twos));
        scaled *= power;
    }
    if (digits > fives) {
        ExactInteger power;
        mpz_ui_pow_ui(power.get_mpz_t(), 5, static_cast<unsigned long>(digits - fives));
        scaled *= power;
    }
    const bool negative = scaled < 0;
    std::string text = (negative ? -scaled : scaled).get_str();
    if (digits == 0) return negative ? "-" + text : text;
    if (text.size() <= digits)
        text.insert(0, digits + 1 - text.size(), '0');
    text.insert(text.size() - digits, 1, '.');
    while (text.size() > 1 && text.back() == '0') text.pop_back();
    if (!text.empty() && text.back() == '.') text.pop_back();
    return negative ? "-" + text : text;
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
    case EvalModApproximationFamily::MultiIntervalChebyshev: return "multi_interval_chebyshev";
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

double upperAdd(double left, double right) {
    if (!std::isfinite(left) || !std::isfinite(right))
        return std::numeric_limits<double>::infinity();
    const double value = left + right;
    return std::isfinite(value)
        ? std::nextafter(value, std::numeric_limits<double>::infinity())
        : std::numeric_limits<double>::infinity();
}

double upperMultiply(double left, double right) {
    if ((!std::isfinite(left) && right != 0.0) || (!std::isfinite(right) && left != 0.0))
        return std::numeric_limits<double>::infinity();
    const double value = left * right;
    if (value == 0.0) return 0.0;
    return std::isfinite(value)
        ? std::nextafter(value, std::numeric_limits<double>::infinity())
        : std::numeric_limits<double>::infinity();
}

double upperDivide(double numerator, double denominator) {
    if (!std::isfinite(numerator) || !std::isfinite(denominator)
        || numerator < 0.0 || denominator <= 0.0)
        return std::numeric_limits<double>::infinity();
    if (numerator == 0.0) return 0.0;
    const double value = numerator / denominator;
    return std::isfinite(value)
        ? std::nextafter(value, std::numeric_limits<double>::infinity())
        : std::numeric_limits<double>::infinity();
}

double deterministicDivideRoundSlotBound(std::size_t polyModulusDegree,
                                         std::size_t secretCoefficientAbsSupport,
                                         std::size_t ciphertextComponents,
                                         double outputScale) {
    if (!polyModulusDegree || !secretCoefficientAbsSupport || !ciphertextComponents
        || !std::isfinite(outputScale) || outputScale <= 0.0)
        return std::numeric_limits<double>::infinity();
    // For each component, coefficientwise divide-and-round leaves a residual
    // polynomial with |r_i| <= 1/2. At every canonical root,
    // |r(zeta)| <= N/2 and |s(zeta)| <= N * support(s_i).
    const double residualEmbedding = 0.5 * static_cast<double>(polyModulusDegree);
    const double secretEmbedding = upperMultiply(
        static_cast<double>(polyModulusDegree),
        static_cast<double>(secretCoefficientAbsSupport));
    double power = 1.0;
    double total = 0.0;
    for (std::size_t component = 0; component < ciphertextComponents; ++component) {
        total = upperAdd(total, upperMultiply(residualEmbedding, power));
        power = upperMultiply(power, secretEmbedding);
    }
    return upperDivide(total, outputScale);
}

double deterministicKeySwitchSlotBound(
    const PreparedEvalModConstants& preparedConstants,
    std::size_t chainIndex,
    std::size_t relinearizationSteps,
    double scale) {
    if (!relinearizationSteps) return 0.0;
    if (!preparedConstants.deterministicEvaluationKeyNoiseSupport
        || !preparedConstants.polyModulusDegree || !preparedConstants.specialPrime
        || chainIndex >= preparedConstants.inputDataPrimes.size()
        || !std::isfinite(scale) || scale <= 0.0)
        return std::numeric_limits<double>::infinity();
    const std::size_t activePrimes = preparedConstants.inputDataPrimes.size() - chainIndex;
    double primeSum = 0.0;
    for (std::size_t index = 0; index < activePrimes; ++index)
        primeSum = upperAdd(primeSum,
                            static_cast<double>(preparedConstants.inputDataPrimes[index] - 1));
    const double degree = static_cast<double>(preparedConstants.polyModulusDegree);
    // Each decomposition limb t_i has coefficient magnitude <= q_i-1. Each
    // evaluation-key error polynomial has coefficient magnitude <= B_e.
    // The canonical embedding of t_i*e_i is at most
    // N(q_i-1) * N B_e. Key switching divides this sum by special prime p.
    double keyNoiseEmbedding = upperMultiply(
        upperMultiply(upperMultiply(degree, degree),
                      static_cast<double>(
                          preparedConstants.evaluationKeyNoiseCoefficientAbsSupport)),
        primeSum);
    keyNoiseEmbedding = upperDivide(
        keyNoiseEmbedding, static_cast<double>(preparedConstants.specialPrime));
    // SEAL's final special-prime ModDown divide-and-round affects both output
    // components; apply the same coefficient-support/canonical-embedding proof.
    const double modDownSlotError = deterministicDivideRoundSlotBound(
        preparedConstants.polyModulusDegree,
        preparedConstants.secretCoefficientAbsSupport, 2, scale);
    const double keyNoiseSlotError = upperDivide(keyNoiseEmbedding, scale);
    return upperMultiply(
        upperAdd(keyNoiseSlotError, modDownSlotError),
        static_cast<double>(relinearizationSteps));
}

std::string upperDecimal(double value) {
    if (!std::isfinite(value)) return "inf";
    std::ostringstream output;
    output << std::setprecision(std::numeric_limits<double>::max_digits10)
           << std::scientific << value;
    return output.str();
}

CertifiedBound deterministicBound(double value, std::string provenance) {
    if (!std::isfinite(value) || value < 0.0)
        throw std::invalid_argument("invalid deterministic EvalMod bound");
    return {value, upperDecimal(value), BoundKind::Deterministic,
            -std::numeric_limits<double>::infinity(), true, std::move(provenance)};
}

CertifiedBound unknownBound(std::string provenance) {
    return {std::numeric_limits<double>::infinity(), "inf", BoundKind::Unknown,
            0.0, false, std::move(provenance)};
}

double exactRationalAbsUpper(const mpq_class& value) {
    MpReal converted;
    mpfr_set_q(converted.get(), value.get_mpq_t(), MPFR_RNDU);
    mpfr_abs(converted.get(), converted.get(), MPFR_RNDU);
    return mpfr_get_d(converted.get(), MPFR_RNDU);
}

double modulusLog2Lower(const ExactInteger& modulus) {
    MpReal converted;
    mpfr_set_z(converted.get(), modulus.get_mpz_t(), MPFR_RNDD);
    mpfr_log2(converted.get(), converted.get(), MPFR_RNDD);
    return mpfr_get_d(converted.get(), MPFR_RNDD);
}

double certifiedHeadroomBits(const EvalModNodeCertificate& certificate) {
    if (!certificate.semanticError.rigorous) return -std::numeric_limits<double>::infinity();
    const double semanticMagnitude = upperAdd(certificate.valueAbs.upperBound,
                                              certificate.semanticError.upperBound);
    if (semanticMagnitude == 0.0) return std::numeric_limits<double>::infinity();
    const double scale = binary64Value(certificate.scale.runtimeBinary64Bits);
    const double encodedMagnitude = upperMultiply(semanticMagnitude, scale);
    if (!std::isfinite(encodedMagnitude) || encodedMagnitude <= 0.0)
        return -std::numeric_limits<double>::infinity();
    return modulusLog2Lower(certificate.modulus) - 1.0 - std::log2(encodedMagnitude);
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

EvalModPolynomial convertChebyshevToMonomial(const EvalModPolynomial& polynomial) {
    if (polynomial.basis != PolynomialBasis::Chebyshev
        || polynomial.decimalCoefficients.empty())
        throw std::invalid_argument("expected non-empty Chebyshev polynomial");
    using RationalPolynomial = std::vector<mpq_class>;
    RationalPolynomial result(polynomial.decimalCoefficients.size(), mpq_class(0));
    RationalPolynomial previous{mpq_class(1)};
    RationalPolynomial current{mpq_class(0), mpq_class(1)};
    for (std::size_t degree = 0; degree < polynomial.decimalCoefficients.size(); ++degree) {
        const auto coefficient = parseExactDecimal(polynomial.decimalCoefficients[degree]);
        const auto& basis = degree == 0 ? previous : current;
        for (std::size_t power = 0; power < basis.size(); ++power)
            result[power] += coefficient * basis[power];
        if (degree + 1 == polynomial.decimalCoefficients.size()) break;
        if (degree == 0) continue;
        RationalPolynomial next(current.size() + 1, mpq_class(0));
        for (std::size_t power = 0; power < current.size(); ++power)
            next[power + 1] += 2 * current[power];
        for (std::size_t power = 0; power < previous.size(); ++power)
            next[power] -= previous[power];
        previous = std::move(current);
        current = std::move(next);
    }
    EvalModPolynomial converted;
    converted.basis = PolynomialBasis::Monomial;
    converted.decimalCoefficients.reserve(result.size());
    for (auto& coefficient : result) {
        coefficient.canonicalize();
        converted.decimalCoefficients.push_back(exactTerminatingDecimal(coefficient));
    }
    while (converted.decimalCoefficients.size() > 1
           && converted.decimalCoefficients.back() == "0")
        converted.decimalCoefficients.pop_back();
    return converted;
}

CompiledEvalModCircuit compileEvalModPolynomial(const EvalModPolynomial& polynomial,
                                               const EvalModProblem& problem,
                                               std::size_t babyStep) {
    if (polynomial.basis == PolynomialBasis::Chebyshev)
        return compileEvalModPolynomial(
            convertChebyshevToMonomial(polynomial), problem, babyStep);
    if (problem.requireCertifiedScaleSchedule && !problem.exactModulusContext)
        throw std::invalid_argument("MissingExactModulusContext");
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
        if (problem.requireCertifiedScaleSchedule) {
            if (!problem.exactModulusContext)
                throw std::invalid_argument("MissingExactModulusContext");
            const auto schedule = buildExactEvalModScaleSchedule(
                compiled, *problem.exactModulusContext, false);
            if (!schedule.valid)
                throw std::runtime_error("certified scale prefix is invalid: " + schedule.failure);
            if (schedule.scaleValues[left].runtimeBinary64Bits
                == schedule.scaleValues[right].runtimeBinary64Bits)
                return std::pair<std::size_t, std::size_t>{left, right};
            const std::size_t correctedChain = targetChain + 1;
            auto correctWithRescale = [&](std::size_t branch, std::size_t target)
                -> std::optional<std::pair<std::size_t, std::size_t>> {
                if (correctedChain >= problem.exactModulusContext->dataPrimes.size())
                    return std::nullopt;
                const auto prime = problem.exactModulusContext->dataPrimes[
                    problem.exactModulusContext->dataPrimes.size() - correctedChain];
                const double branchScale = binary64Value(
                    schedule.scaleValues[branch].runtimeBinary64Bits);
                const double targetScale = binary64Value(
                    schedule.scaleValues[target].runtimeBinary64Bits);
                const auto correction = identityCorrectionScale(branchScale, targetScale, prime);
                if (!correction) return std::nullopt;
                const auto correctionScale = ExactScale::fromBinaryDouble(*correction);
                const auto constant = addNode(EvalModOperation::EncodeConstant, {}, targetChain,
                                              correctionScale, correctionScale, "1");
                depths.resize(compiled.nodes.size()); depths[constant] = depths[branch];
                const auto multipliedScale = multiplyScale(
                    schedule.scaleValues[branch].ideal, correctionScale);
                auto multiplied = addNode(EvalModOperation::MultiplyPlain, {branch, constant},
                                          targetChain, schedule.scaleValues[branch].ideal,
                                          multipliedScale);
                depths.resize(compiled.nodes.size()); depths[multiplied] = depths[branch];
                auto corrected = addNode(EvalModOperation::Rescale, {multiplied}, correctedChain,
                                         multipliedScale, divideScale(multipliedScale, prime));
                depths.resize(compiled.nodes.size()); depths[corrected] = depths[branch];
                auto loweredTarget = addNode(EvalModOperation::ModSwitch, {target, corrected},
                                             correctedChain, schedule.scaleValues[target].ideal,
                                             schedule.scaleValues[target].ideal);
                depths.resize(compiled.nodes.size()); depths[loweredTarget] = depths[target];
                return std::pair<std::size_t, std::size_t>{loweredTarget, corrected};
            };
            auto correctWithoutRescale = [&](std::size_t branch, std::size_t target)
                -> std::optional<std::pair<std::size_t, std::size_t>> {
                const double branchScale = binary64Value(
                    schedule.scaleValues[branch].runtimeBinary64Bits);
                const double targetScale = binary64Value(
                    schedule.scaleValues[target].runtimeBinary64Bits);
                const auto correction = identityMultiplyScale(branchScale, targetScale);
                if (!correction) return std::nullopt;
                const auto correctionScale = ExactScale::fromBinaryDouble(*correction);
                const auto constant = addNode(EvalModOperation::EncodeConstant, {}, targetChain,
                                              correctionScale, correctionScale, "1");
                depths.resize(compiled.nodes.size()); depths[constant] = depths[branch];
                const auto multipliedScale = multiplyScale(
                    schedule.scaleValues[branch].ideal, correctionScale);
                auto corrected = addNode(EvalModOperation::MultiplyPlain, {branch, constant},
                                         targetChain, schedule.scaleValues[branch].ideal,
                                         multipliedScale);
                depths.resize(compiled.nodes.size()); depths[corrected] = depths[branch];
                return std::pair<std::size_t, std::size_t>{target, corrected};
            };
            if (auto corrected = correctWithoutRescale(right, left)) return *corrected;
            if (auto corrected = correctWithoutRescale(left, right))
                return std::pair<std::size_t, std::size_t>{corrected->second, corrected->first};
            if (auto corrected = correctWithRescale(right, left)) return *corrected;
            if (auto corrected = correctWithRescale(left, right))
                return std::pair<std::size_t, std::size_t>{corrected->second, corrected->first};
            throw std::runtime_error("ScaleScheduleInfeasible: exact binary64 identity correction not found");
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
            std::size_t term{};
            if (offset > 0 && coefficientValue == 1.0L) {
                term = powers[offset];
            } else {
                const auto constant = addNode(
                    EvalModOperation::EncodeConstant, {}, block, workingScale, workingScale,
                    polynomial.decimalCoefficients[coefficient]);
                term = constant;
                if (offset > 0) {
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
                    if (problem.requireCertifiedScaleSchedule) {
                        if (!problem.exactModulusContext)
                            throw std::invalid_argument("MissingExactModulusContext");
                        const auto schedule = buildExactEvalModScaleSchedule(
                            compiled, *problem.exactModulusContext, false);
                        if (!schedule.valid)
                            throw std::runtime_error(
                                "certified AddPlain scale prefix is invalid: " + schedule.failure);
                        compiled.nodes[plain].chainIndex = compiled.nodes[cipher].chainIndex;
                        compiled.nodes[plain].inputScale =
                            schedule.scaleValues[cipher].runtimeExact;
                        compiled.nodes[plain].outputScale =
                            schedule.scaleValues[cipher].runtimeExact;
                    }
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
    result.scaleValues.resize(circuit.nodes.size());
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
            result.scaleValues[index] = initialScaleValue(node.outputScale);
            result.nodeScales[index] = result.scaleValues[index].ideal;
            continue;
        }
        if (node.inputs.empty()) return reject(index, "operation_has_no_input");
        const auto& left = result.scaleValues[node.inputs[0]];
        if (node.operation == EvalModOperation::Relinearize) {
            result.scaleValues[index] = left;
        } else if (node.operation == EvalModOperation::ModSwitch) {
            if (node.inputs.size() != 2
                || circuit.nodes[node.inputs[0]].chainIndex >= node.chainIndex
                || circuit.nodes[node.inputs[1]].chainIndex != node.chainIndex)
                return reject(index, "invalid_modswitch_target_level");
            result.scaleValues[index] = left;
        } else if (node.operation == EvalModOperation::Rescale) {
            if (node.chainIndex == 0)
                return reject(index, "rescale_does_not_advance_chain");
            const auto source = node.inputs[0];
            if (circuit.nodes[source].chainIndex + 1 != node.chainIndex)
                return reject(index, "rescale_chain_step_is_not_one");
            const std::size_t droppedIndex = context.dataPrimes.size() - node.chainIndex;
            const auto droppedPrime = context.dataPrimes[droppedIndex];
            result.rescalePrimes[index] = droppedPrime;
            result.scaleValues[index] = divideScaleValue(left, droppedPrime);
        } else if (node.operation == EvalModOperation::AlignScale) {
            if (rigorous) return reject(index, "metadata_scale_alignment_prohibited");
            if (node.inputs.size() != 2)
                return reject(index, "metadata_scale_alignment_target_missing");
            result.scaleValues[index] = result.scaleValues[node.inputs[1]];
        } else {
            if (node.inputs.size() != 2)
                return reject(index, "binary_operation_input_count_mismatch");
            const auto& right = result.scaleValues[node.inputs[1]];
            if (node.operation == EvalModOperation::Add
                || node.operation == EvalModOperation::AddPlain) {
                if (rigorous && left.runtimeBinary64Bits != right.runtimeBinary64Bits)
                    return reject(index, "addition_scale_mismatch");
                result.scaleValues[index] = left;
            } else if (node.operation == EvalModOperation::MultiplyPlain
                       || node.operation == EvalModOperation::MultiplyCipher) {
                result.scaleValues[index] = multiplyScaleValue(left, right);
            } else {
                return reject(index, "unsupported_evalmod_operation");
            }
        }
        result.nodeScales[index] = result.scaleValues[index].ideal;
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
    if (!schedule.valid || schedule.scaleValues.size() != circuit.nodes.size()
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
            const auto& scale = schedule.scaleValues[index].runtimeExact;
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
    result.specialPrime = adapter.specialKeyModulusValue();
    result.polyModulusDegree = 2 * adapter.slotCount();
    result.secretCoefficientAbsSupport = 1;
#ifdef SEAL_USE_GAUSSIAN_NOISE
    result.evaluationKeyNoiseCoefficientAbsSupport = 0;
    result.deterministicEvaluationKeyNoiseSupport = false;
#else
    // Vendored SEAL's centered-binomial sampler is the difference of two
    // 21-bit Hamming weights, hence its exact coefficient support is [-21, 21].
    result.evaluationKeyNoiseCoefficientAbsSupport = 21;
    result.deterministicEvaluationKeyNoiseSupport = true;
#endif
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
                ? runtimeSchedule.scaleValues[consumer.inputs[0]].runtimeExact
                : circuit.nodes[input].outputScale;
            if (constantScales[input] && !equalScale(*constantScales[input], encodingScale))
                throw std::invalid_argument("EvalMod constant requires multiple encoding scales");
            constantScales[input] = encodingScale;
        }
    }
    for (std::size_t index = 0; index < circuit.nodes.size(); ++index) {
        const auto& node = circuit.nodes[index];
        if (node.operation != EvalModOperation::EncodeConstant) continue;
        // The compiler may retain a coefficient node that was optimized out of
        // the ciphertext path (for example, coefficient 1 for the x term).
        // Such a node has no runtime plaintext and therefore needs no prepared
        // representation.
        if (!constantLevels[index]) continue;
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

EvalModArithmeticCertificate certifyEvalModDagArithmetic(
    const CompiledEvalModCircuit& circuit,
    const EvalModExactScaleSchedule& schedule,
    const PreparedEvalModConstants& preparedConstants,
    double inputValueAbsUpperBound,
    double inputSemanticErrorUpperBound) {
    EvalModArithmeticCertificate result;
    result.outputError = unknownBound("arithmetic certificate not constructed");
    if (!isCompiledEvalModCircuitValid(circuit)
        || !schedule.available || !schedule.valid || !schedule.rigorous
        || schedule.scaleValues.size() != circuit.nodes.size()
        || schedule.nodeModuli.size() != circuit.nodes.size()
        || !preparedConstants.rigorous
        || !preparedConstants.polyModulusDegree
        || !preparedConstants.secretCoefficientAbsSupport
        || preparedConstants.specialPrime < 2
        || !std::isfinite(inputValueAbsUpperBound) || inputValueAbsUpperBound < 0.0
        || !std::isfinite(inputSemanticErrorUpperBound) || inputSemanticErrorUpperBound < 0.0) {
        result.status = !schedule.available || !schedule.valid || !schedule.rigorous
            ? EvalModCertificationStatus::InvalidScaleSchedule
            : EvalModCertificationStatus::InvalidInput;
        result.detail = "invalid rigorous arithmetic certificate input";
        return result;
    }

    ExactInteger expectedInputModulus = 1;
    for (const auto prime : preparedConstants.inputDataPrimes)
        expectedInputModulus *= exactIntegerFromUint64(prime);
    if (preparedConstants.inputDataPrimes.empty()
        || expectedInputModulus != schedule.nodeModuli.front()) {
        result.status = EvalModCertificationStatus::PreparedConstantMismatch;
        result.detail = "prepared constants and exact modulus schedule use different contexts";
        return result;
    }
    result.keyNoiseMetadata = {
#ifdef SEAL_USE_GAUSSIAN_NOISE
        EvaluationKeySamplerKind::ClippedRoundedGaussian,
#else
        EvaluationKeySamplerKind::CenteredBinomial,
#endif
        3.2,
        preparedConstants.evaluationKeyNoiseCoefficientAbsSupport,
        preparedConstants.secretCoefficientAbsSupport,
        preparedConstants.polyModulusDegree,
        preparedConstants.specialPrime,
        preparedConstants.inputDataPrimes,
        0,
        preparedConstants.deterministicEvaluationKeyNoiseSupport,
        preparedConstants.deterministicEvaluationKeyNoiseSupport
            ? "vendored SEAL centered-binomial support [-21,21] and ternary secret support [-1,1]"
            : "finite deterministic support for configured SEAL sampler unavailable"
    };

    std::vector<const PreparedEvalModConstant*> preparedByNode(circuit.nodes.size(), nullptr);
    for (const auto& prepared : preparedConstants.constants) {
        if (!prepared.rigorous || prepared.node >= circuit.nodes.size()
            || circuit.nodes[prepared.node].operation != EvalModOperation::EncodeConstant
            || prepared.decimal != circuit.nodes[prepared.node].constantDecimal
            || !std::isfinite(prepared.encodingErrorUpperBound)
            || prepared.encodingErrorUpperBound < 0.0 || preparedByNode[prepared.node]) {
            result.status = EvalModCertificationStatus::PreparedConstantMismatch;
            result.detail = "invalid prepared constant certificate";
            return result;
        }
        preparedByNode[prepared.node] = &prepared;
    }

    std::vector<bool> reachable(circuit.nodes.size(), false);
    reachable[circuit.outputNode] = true;
    for (std::size_t index = circuit.nodes.size(); index-- > 0;)
        if (reachable[index])
            for (const auto input : circuit.nodes[index].inputs) reachable[input] = true;
    for (std::size_t index = 0; index < circuit.nodes.size(); ++index)
        if (reachable[index] && circuit.nodes[index].operation == EvalModOperation::EncodeConstant
            && !preparedByNode[index]) {
            result.status = EvalModCertificationStatus::PreparedConstantMismatch;
            result.failingNode = index;
            result.detail = "reachable EvalMod constant has no prepared encoding certificate";
            return result;
        }

    result.nodes.resize(circuit.nodes.size());
    std::vector<std::size_t> ciphertextComponents(circuit.nodes.size(), 0);
    bool operationBoundUnavailable = false;
    for (std::size_t index = 0; index < circuit.nodes.size(); ++index) {
        const auto& node = circuit.nodes[index];
        auto& certificate = result.nodes[index];
        certificate.scale = schedule.scaleValues[index];
        certificate.modulus = schedule.nodeModuli[index];
        if (!reachable[index]) {
            certificate.valueAbs = deterministicBound(0.0, "unreachable DAG node");
            certificate.semanticError = deterministicBound(0.0, "unreachable DAG node");
            certificate.localAddedError = deterministicBound(0.0, "unreachable DAG node");
            certificate.headroomBits = std::numeric_limits<double>::infinity();
            continue;
        }

        if (node.operation == EvalModOperation::Input) {
            ciphertextComponents[index] = 2;
            certificate.valueAbs = deterministicBound(
                inputValueAbsUpperBound, "certified EvalMod input value bound");
            certificate.semanticError = deterministicBound(
                inputSemanticErrorUpperBound, "upstream CoeffToSlot input error certificate");
            certificate.localAddedError = deterministicBound(0.0, "input adds no local error");
        } else if (node.operation == EvalModOperation::EncodeConstant) {
            const auto intended = parseExactDecimal(node.constantDecimal);
            certificate.valueAbs = deterministicBound(
                exactRationalAbsUpper(intended), "exact decimal coefficient magnitude");
            certificate.semanticError = deterministicBound(
                preparedByNode[index]->encodingErrorUpperBound,
                "prepared RNS plaintext coefficient quantization certificate");
            certificate.localAddedError = certificate.semanticError;
        } else if (node.operation == EvalModOperation::Add
                   || node.operation == EvalModOperation::AddPlain) {
            const auto& left = result.nodes[node.inputs[0]];
            const auto& right = result.nodes[node.inputs[1]];
            ciphertextComponents[index] = std::max(
                ciphertextComponents[node.inputs[0]],
                ciphertextComponents[node.inputs[1]]);
            certificate.valueAbs = deterministicBound(
                upperAdd(left.valueAbs.upperBound, right.valueAbs.upperBound),
                "M_out = M_a + M_b");
            certificate.localAddedError = deterministicBound(
                0.0, "RNS addition adds no numerical error");
            if (left.semanticError.rigorous && right.semanticError.rigorous)
                certificate.semanticError = deterministicBound(
                    upperAdd(left.semanticError.upperBound, right.semanticError.upperBound),
                    "E_out = E_a + E_b");
            else
                certificate.semanticError = unknownBound(
                    "addition depends on an unknown operand error");
        } else if (node.operation == EvalModOperation::MultiplyPlain
                   || node.operation == EvalModOperation::MultiplyCipher) {
            const auto& left = result.nodes[node.inputs[0]];
            const auto& right = result.nodes[node.inputs[1]];
            ciphertextComponents[index] = node.operation == EvalModOperation::MultiplyCipher
                ? ciphertextComponents[node.inputs[0]]
                    + ciphertextComponents[node.inputs[1]] - 1
                : std::max(ciphertextComponents[node.inputs[0]],
                           ciphertextComponents[node.inputs[1]]);
            certificate.valueAbs = deterministicBound(
                upperMultiply(left.valueAbs.upperBound, right.valueAbs.upperBound),
                "M_out = M_a M_b");
            certificate.localAddedError = deterministicBound(
                0.0, "ring multiplication adds no separate semantic error");
            if (left.semanticError.rigorous && right.semanticError.rigorous) {
                const double cross = upperAdd(
                    upperMultiply(left.valueAbs.upperBound, right.semanticError.upperBound),
                    upperMultiply(right.valueAbs.upperBound, left.semanticError.upperBound));
                certificate.semanticError = deterministicBound(
                    upperAdd(cross, upperMultiply(left.semanticError.upperBound,
                                                  right.semanticError.upperBound)),
                    node.operation == EvalModOperation::MultiplyPlain
                        ? "|c|E_x + M_x|delta_c| + E_x|delta_c|"
                        : "M_aE_b + M_bE_a + E_aE_b");
            } else {
                certificate.semanticError = unknownBound(
                    "multiplication depends on an unknown operand error");
            }
        } else {
            const auto& input = result.nodes[node.inputs[0]];
            certificate.valueAbs = input.valueAbs;
            if (node.operation == EvalModOperation::Relinearize) {
                const auto inputComponents = ciphertextComponents[node.inputs[0]];
                ciphertextComponents[index] = 2;
                const auto steps = inputComponents > 2 ? inputComponents - 2 : 0;
                result.keyNoiseMetadata.relevantEvaluationKeyComponents += steps;
                const double local = deterministicKeySwitchSlotBound(
                    preparedConstants, node.chainIndex, steps,
                    binary64Value(certificate.scale.runtimeBinary64Bits));
                if (std::isfinite(local) && input.semanticError.rigorous) {
                    certificate.localAddedError = deterministicBound(
                        local,
                        "SEAL CBD/ternary-support key-switch and special-prime ModDown bound");
                    certificate.semanticError = deterministicBound(
                        upperAdd(input.semanticError.upperBound, local),
                        "E_out = E_in + deterministic key-switch bound");
                } else {
                    certificate.localAddedError = unknownBound(
                        "finite-support SEAL key-switch bound unavailable");
                    certificate.semanticError = unknownBound(
                        "relinearization depends on an unknown key-switch bound");
                    if (!operationBoundUnavailable) {
                        result.status = EvalModCertificationStatus::RigorousKeySwitchBoundUnavailable;
                        result.failingNode = index;
                        result.detail = "deterministic/probabilistic SEAL key-switch certificate unavailable";
                    }
                    operationBoundUnavailable = true;
                }
            } else if (node.operation == EvalModOperation::Rescale) {
                ciphertextComponents[index] = ciphertextComponents[node.inputs[0]];
                const double local = deterministicDivideRoundSlotBound(
                    preparedConstants.polyModulusDegree,
                    preparedConstants.secretCoefficientAbsSupport,
                    ciphertextComponents[index],
                    binary64Value(certificate.scale.runtimeBinary64Bits));
                if (std::isfinite(local) && input.semanticError.rigorous) {
                    certificate.localAddedError = deterministicBound(
                        local,
                        "SEAL divide_and_round_q_last_ntt coefficient support and canonical-embedding bound");
                    certificate.semanticError = deterministicBound(
                        upperAdd(input.semanticError.upperBound, local),
                        "E_out = E_in + deterministic rescale residual bound");
                } else {
                    certificate.localAddedError = unknownBound(
                        "finite SEAL divide-and-round residual bound unavailable");
                    certificate.semanticError = unknownBound("rescale semantic error is unknown");
                    if (!operationBoundUnavailable) {
                        result.status = EvalModCertificationStatus::RigorousRescaleBoundUnavailable;
                        result.failingNode = index;
                        result.detail = "SEAL CKKS rescale certificate unavailable";
                    }
                    operationBoundUnavailable = true;
                }
            } else if (node.operation == EvalModOperation::AlignScale) {
                ciphertextComponents[index] = ciphertextComponents[node.inputs[0]];
                certificate.localAddedError = unknownBound(
                    "metadata-only scale alignment is prohibited");
                certificate.semanticError = unknownBound(
                    "metadata-only scale alignment invalidates semantics");
                if (!operationBoundUnavailable) {
                    result.status = EvalModCertificationStatus::ScaleScheduleInfeasible;
                    result.failingNode = index;
                    result.detail = "certified DAG contains AlignScale";
                }
                operationBoundUnavailable = true;
            } else {
                ciphertextComponents[index] = ciphertextComponents[node.inputs[0]];
                certificate.localAddedError = deterministicBound(
                    0.0, "ModSwitch preserves represented integer when no-wrap is proved");
                certificate.semanticError = input.semanticError;
            }
        }
        certificate.headroomBits = certifiedHeadroomBits(certificate);
        if (node.operation == EvalModOperation::ModSwitch
            && certificate.semanticError.rigorous && certificate.headroomBits <= 0.0
            && !operationBoundUnavailable) {
            result.status = EvalModCertificationStatus::ModSwitchHeadroomViolation;
            result.failingNode = index;
            result.detail = "ModSwitch represented value is not proved inside centered target modulus";
            operationBoundUnavailable = true;
        }
    }

    result.outputError = result.nodes[circuit.outputNode].semanticError;
    result.rigorous = !operationBoundUnavailable && result.outputError.rigorous;
    result.log2FailureProbability = result.rigorous
        ? -std::numeric_limits<double>::infinity() : 0.0;
    if (result.rigorous) {
        result.status = EvalModCertificationStatus::Certified;
        result.detail = "all DAG arithmetic bounds are deterministic";
    } else if (!operationBoundUnavailable) {
        result.status = EvalModCertificationStatus::UnknownOperationBound;
        result.detail = "output depends on an unknown arithmetic bound";
    }
    return result;
}

EvalModSynthesisResult synthesizeEvalMod(const EvalModProblem& problem) {
    if (problem.qSource <= 1 || problem.coeffToSlotScale.numerator <= 0
        || problem.coeffToSlotScale.denominator <= 0 || problem.outputScale.numerator <= 0
        || problem.outputScale.denominator <= 0 || problem.availableLevels == 0
        || problem.targetPrecisionBits < 32
        || !std::isfinite(problem.targetAbsoluteError) || problem.targetAbsoluteError <= 0.0
        || problem.minimumSearchDegree < 1
        || problem.minimumSearchDegree > problem.maximumSearchDegree
        || problem.minimumBabyStep < 2
        || problem.minimumBabyStep > problem.maximumBabyStep
        || problem.maximumBabyStep > problem.maximumSearchDegree) {
        throw std::invalid_argument("invalid EvalMod synthesis problem");
    }
    EvalModSynthesisResult result{problem, estimateEvalModDomain(problem.ciphertextModel), {}, {}};
    struct Spec { EvalModApproximationFamily family; std::size_t degree; std::size_t babyStep; };
    std::vector<Spec> specs{{EvalModApproximationFamily::PeriodicSineBaseline, 9, 3},
                            {EvalModApproximationFamily::MultiIntervalLeastSquaresPrototype, 15, 4}};
    std::size_t firstDegree = problem.minimumSearchDegree;
    if ((firstDegree & 1U) == 0) ++firstDegree;
    for (const auto family : {EvalModApproximationFamily::MultiIntervalMinimax,
                              EvalModApproximationFamily::MultiIntervalChebyshev}) {
        for (std::size_t degree = firstDegree; degree <= problem.maximumSearchDegree;) {
            for (std::size_t babyStep = problem.minimumBabyStep;
                 babyStep <= std::min(problem.maximumBabyStep, degree); ++babyStep)
                specs.push_back({family, degree, babyStep});
            if (degree > problem.maximumSearchDegree - std::min<std::size_t>(2, degree)) break;
            degree += 2;
        }
    }
    specs.push_back({EvalModApproximationFamily::MinimaxInverseSine, 35, 5});
    const double searchOutputGain = exactScaleUp(
        ExactScale::rational(problem.qSource * problem.outputScale.denominator,
                             problem.outputScale.numerator));
    double normalizedApproximationBudget = problem.targetAbsoluteError
        - problem.slotToCoeffAdditive - problem.finalAdditive;
    if (normalizedApproximationBudget > 0.0
        && problem.slotToCoeffOperatorNorm > 0.0 && searchOutputGain > 0.0) {
        normalizedApproximationBudget /= problem.slotToCoeffOperatorNorm * searchOutputGain;
        normalizedApproximationBudget -=
            problem.ciphertextModel.normalizedCoeffToSlotErrorAbsBound;
        normalizedApproximationBudget -= result.domain.errors.periodMismatch;
    } else {
        normalizedApproximationBudget = 0.0;
    }
    std::optional<std::size_t> minimaxTargetDegree;
    std::optional<std::size_t> chebyshevTargetDegree;
    for (const auto& spec : specs) {
        const auto reachedDegree = spec.family == EvalModApproximationFamily::MultiIntervalMinimax
            ? minimaxTargetDegree
            : spec.family == EvalModApproximationFamily::MultiIntervalChebyshev
                ? chebyshevTargetDegree : std::optional<std::size_t>{};
        if (reachedDegree && spec.degree > *reachedDegree) continue;
        EvalModCandidate candidate;
        candidate.family = spec.family;
        candidate.domain = result.domain;
        if (spec.family == EvalModApproximationFamily::MultiIntervalMinimax
            || spec.family == EvalModApproximationFamily::MultiIntervalChebyshev) {
            try {
                const bool chebyshev = spec.family
                    == EvalModApproximationFamily::MultiIntervalChebyshev;
                auto remez = remezOdd(result.domain, spec.degree, chebyshev);
                candidate.polynomial = chebyshev
                    ? convertChebyshevToMonomial(remez.polynomial)
                    : std::move(remez.polynomial);
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
                    const auto inputBits = exactScaleBitsUp(
                        candidate.exactScaleSchedule.scaleValues[source].runtimeExact);
                    const auto outputBits = exactScaleBitsUp(
                        candidate.exactScaleSchedule.scaleValues[index].runtimeExact);
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
            1.0, gain, problem.slotToCoeffOperatorNorm,
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
                candidate.nodeErrorStates, false);
        }
        candidate.polynomialArithmeticError =
            candidate.nodeErrorStates[candidate.compiledCircuit.outputNode].absoluteErrorBound;
        // This model is empirical/diagnostic even if legacy callers set its
        // historical `rigorous` bit. Certified status can only come from
        // `certifyEvalModDagArithmetic` and its per-operation certificates.
        candidate.arithmeticErrorRigorous = false;
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
        candidate.arithmeticErrorCertified = false;
        candidate.analyticalArithmeticBound = candidate.polynomialArithmeticError;
        candidate.analyticalArithmeticBoundRigorous = false;
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
        if (candidate.approximationConverged && candidate.intervalCertified
            && candidate.satisfiesLevelBudget && normalizedApproximationBudget > 0.0
            && candidate.intervalCertificate.approximationErrorUpperBoundDouble
                <= normalizedApproximationBudget) {
            if (candidate.family == EvalModApproximationFamily::MultiIntervalMinimax)
                minimaxTargetDegree = spec.degree;
            else if (candidate.family == EvalModApproximationFamily::MultiIntervalChebyshev)
                chebyshevTargetDegree = spec.degree;
        }
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

static Cipher executeEvalModCircuitImpl(
    SealAdapter& adapter, const CompiledEvalModCircuit& circuit,
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
        std::vector<bool> requiredByCipherOperation(circuit.nodes.size(), false);
        for (const auto& consumer : circuit.nodes)
            for (const auto input : consumer.inputs)
                if (circuit.nodes[input].operation == EvalModOperation::EncodeConstant)
                    requiredByCipherOperation[input] = true;
        for (std::size_t index = 0; index < circuit.nodes.size(); ++index)
            if (requiredByCipherOperation[index] && !preparedByNode[index])
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

Cipher executeEvalModCircuitDiagnostic(
    SealAdapter& adapter, const CompiledEvalModCircuit& circuit,
    const Cipher& coeffToSlotOutput, EvalModExecutionTrace* trace,
    const std::vector<double>* semanticInput,
    std::vector<EvalModNodeDifferential>* differentialTrace) {
    return executeEvalModCircuitImpl(adapter, circuit, coeffToSlotOutput, trace,
                                     semanticInput, differentialTrace, nullptr);
}

Cipher executePreparedEvalMod(
    SealAdapter& adapter, const CompiledEvalModCircuit& circuit,
    const PreparedEvalModConstants& preparedConstants,
    const Cipher& coeffToSlotOutput, EvalModExecutionTrace* trace,
    const std::vector<double>* semanticInput,
    std::vector<EvalModNodeDifferential>* differentialTrace) {
    return executeEvalModCircuitImpl(adapter, circuit, coeffToSlotOutput, trace,
                                     semanticInput, differentialTrace, &preparedConstants);
}

EvalModCoeffToSlotResult executeEvalModAfterCoeffToSlot(
    SealAdapter& adapter, const CompiledEvalModCircuit& circuit,
    const PreparedEvalModConstants& preparedConstants,
    CoeffToSlotResult coeffToSlotOutput) {
    EvalModCoeffToSlotResult result;
    result.slotCipherFirst = executePreparedEvalMod(
        adapter, circuit, preparedConstants, coeffToSlotOutput.slotCipherFirst,
        &result.firstTrace);
    result.slotCipherSecond = executePreparedEvalMod(
        adapter, circuit, preparedConstants, coeffToSlotOutput.slotCipherSecond,
        &result.secondTrace);
    return result;
}

PreparedEvalModPlan prepareEvalMod(
    SealAdapter& adapter,
    const EvalModCandidate& candidate,
    const EvalModProblem& problem,
    const Cipher& evalModInput) {
    PreparedEvalModPlan plan;
    plan.contextFingerprint = adapter.contextFingerprint();
    plan.inputChainIndex = adapter.chainIndex(evalModInput);
    plan.inputScaleBinary64Bits = binary64Bits(adapter.scale(evalModInput));
    plan.domain = candidate.domain;
    plan.normalizationGain = candidate.compiledCircuit.normalizationGain;
    plan.denormalizationGain = candidate.compiledCircuit.denormalizationGain;
    plan.approximation = candidate.polynomial;
    plan.circuit = candidate.compiledCircuit;
    plan.dataPrimes = adapter.coeffModulusValues(evalModInput);
    plan.specialPrime = adapter.specialKeyModulusValue();
    plan.levelsConsumed = candidate.compiledCircuit.cost.levelConsumption;
    plan.securityBits = adapter.securityLevelBits();
    if (!isEvalModDomainValid(candidate.domain)) {
        plan.status = EvalModCertificationStatus::DiscontinuityMarginViolation;
        plan.detail = "expanded EvalMod domain reaches or crosses a rounding discontinuity";
        return plan;
    }
    plan.inputError = deterministicBound(
        problem.ciphertextModel.normalizedCoeffToSlotErrorAbsBound,
        "upstream normalized CoeffToSlot error certificate");
    plan.discontinuityMargin = deterministicBound(
        candidate.domain.discontinuityMargin(),
        "0.5 minus the MPFR outward-rounded expanded residual domain");

    const EvalModExactModulusContext context{plan.dataPrimes, plan.specialPrime};
    plan.scaleSchedule = buildExactEvalModScaleSchedule(
        plan.circuit, context, true);
    if (!plan.scaleSchedule.valid || !plan.scaleSchedule.rigorous) {
        plan.status = EvalModCertificationStatus::ScaleScheduleInfeasible;
        plan.detail = plan.scaleSchedule.failure;
        return plan;
    }
    if (plan.levelsConsumed >= plan.dataPrimes.size()
        || plan.levelsConsumed > problem.availableLevels) {
        plan.status = EvalModCertificationStatus::InsufficientLevels;
        plan.detail = "compiled EvalMod DAG exceeds available data-prime levels";
        return plan;
    }
    try {
        plan.constants = prepareEvalModConstants(adapter, plan.circuit, evalModInput);
    } catch (const std::exception& error) {
        plan.status = EvalModCertificationStatus::PreparedConstantMismatch;
        plan.detail = error.what();
        return plan;
    }
    const double normalization = exactScaleUp(plan.normalizationGain);
    const double normalizedInputBound = static_cast<double>(candidate.domain.integerBound)
        + candidate.domain.normalizedResidualBound;
    plan.arithmeticCertificate = certifyEvalModDagArithmetic(
        plan.circuit, plan.scaleSchedule, plan.constants,
        normalizedInputBound / normalization,
        problem.ciphertextModel.normalizedCoeffToSlotErrorAbsBound / normalization);
    plan.arithmeticError = plan.arithmeticCertificate.outputError;
    if (!plan.arithmeticCertificate.rigorous) {
        plan.status = plan.arithmeticCertificate.status;
        plan.detail = plan.arithmeticCertificate.detail;
        return plan;
    }
    if (!candidate.approximationConverged || !candidate.intervalCertificate.proved) {
        plan.status = EvalModCertificationStatus::ApproximationInsufficient;
        plan.detail = "expanded-domain executable polynomial approximation is not certified";
        return plan;
    }
    plan.approximationError = deterministicBound(
        candidate.intervalCertificate.approximationErrorUpperBoundDouble,
        "MPFR interval certificate on the expanded multi-interval domain");
    plan.normalizedEvalModError = deterministicBound(
        upperAdd(plan.approximationError.upperBound, plan.arithmeticError.upperBound),
        "certified approximation plus DAG arithmetic error");
    const double outputGain = exactScaleUp(plan.denormalizationGain);
    plan.denormalizedEvalModError = deterministicBound(
        upperMultiply(plan.normalizedEvalModError.upperBound, outputGain),
        "normalized EvalMod error multiplied by exact output gain");

    const auto subtractDown = [](double left, double right) {
        return std::nextafter(left - right, -std::numeric_limits<double>::infinity());
    };
    const auto divideDown = [](double numerator, double denominator) {
        return std::nextafter(numerator / denominator, 0.0);
    };
    double remaining = subtractDown(problem.targetAbsoluteError,
                                    problem.slotToCoeffAdditive);
    remaining = subtractDown(remaining, problem.finalAdditive);
    if (!(remaining > 0.0) || !(problem.slotToCoeffOperatorNorm > 0.0)
        || !(outputGain > 0.0)) {
        plan.status = EvalModCertificationStatus::NoEvalModErrorBudgetRemaining;
        plan.detail = "SlotToCoeff/final additive errors consume the bootstrap target";
        return plan;
    }
    remaining = divideDown(
        remaining, upperMultiply(problem.slotToCoeffOperatorNorm, outputGain));
    remaining = subtractDown(
        remaining, problem.ciphertextModel.normalizedCoeffToSlotErrorAbsBound);
    remaining = subtractDown(remaining, candidate.domain.errors.periodMismatch);
    plan.normalizedEvalModBudget = remaining;
    if (!(remaining > 0.0)) {
        plan.status = EvalModCertificationStatus::NoEvalModErrorBudgetRemaining;
        plan.detail = "no normalized EvalMod budget remains after input and period errors";
        return plan;
    }

    const double normalizedContribution = upperAdd(
        upperAdd(problem.ciphertextModel.normalizedCoeffToSlotErrorAbsBound,
                 candidate.domain.errors.periodMismatch),
        plan.normalizedEvalModError.upperBound);
    plan.bootstrapContribution = deterministicBound(
        upperAdd(
            upperMultiply(
                upperMultiply(problem.slotToCoeffOperatorNorm, outputGain),
                normalizedContribution),
            upperAdd(problem.slotToCoeffAdditive, problem.finalAdditive)),
        "full real-domain bootstrap propagation bound");

    plan.log2FailureProbability = problem.ciphertextModel.tailModel == TailModel::Deterministic
        ? -std::numeric_limits<double>::infinity()
        : candidate.domain.failureProbabilityLog2;
    if (plan.log2FailureProbability > -128.0) {
        plan.status = EvalModCertificationStatus::FailureProbabilityExceeded;
        plan.detail = "combined domain/operation failure probability exceeds 2^-128";
        return plan;
    }
    if (plan.securityBits < 128) {
        plan.status = EvalModCertificationStatus::SecurityBudgetExceeded;
        plan.detail = "SEAL context was not validated at 128-bit security";
        return plan;
    }
    if (plan.approximationError.upperBound > remaining) {
        plan.status = EvalModCertificationStatus::ApproximationInsufficient;
        plan.detail = "expanded-domain approximation alone exceeds normalized EvalMod budget";
        return plan;
    }
    if (plan.normalizedEvalModError.upperBound > remaining) {
        plan.status = EvalModCertificationStatus::ArithmeticNoiseTooLarge;
        plan.detail = "deterministic arithmetic certificate exceeds remaining normalized budget";
        return plan;
    }
    plan.status = EvalModCertificationStatus::Certified;
    plan.rigorous = true;
    plan.detail = "all prepared EvalMod certificate gates passed";
    return plan;
}

EvalModPlanSearchResult prepareEvalModSearch(
    SealAdapter& adapter,
    const EvalModSynthesisResult& synthesis,
    const Cipher& evalModInput) {
    EvalModPlanSearchResult result;
    result.minimumDegree = std::numeric_limits<std::size_t>::max();
    auto certifiedProblem = synthesis.problem;
    certifiedProblem.exactModulusContext = EvalModExactModulusContext{
        adapter.coeffModulusValues(evalModInput), adapter.specialKeyModulusValue()};
    certifiedProblem.requireCertifiedScaleSchedule = true;
    bool sawCandidate = false;
    bool firstFailureRecorded = false;
    bool everyFailureWasResource = true;
    bool everyFailureWasUnavailableBound = true;
    double selectedLatency = std::numeric_limits<double>::infinity();
    for (const auto& sourceCandidate : synthesis.candidates) {
        if (sourceCandidate.family
                == EvalModApproximationFamily::MultiIntervalLeastSquaresPrototype
            || !sourceCandidate.approximationConverged)
            continue;
        sawCandidate = true;
        if (std::find(result.familiesSearched.begin(), result.familiesSearched.end(),
                      sourceCandidate.family) == result.familiesSearched.end())
            result.familiesSearched.push_back(sourceCandidate.family);
        result.minimumDegree = std::min(
            result.minimumDegree, sourceCandidate.compiledCircuit.cost.degree);
        result.maximumDegree = std::max(
            result.maximumDegree, sourceCandidate.compiledCircuit.cost.degree);
        result.maximumDepth = std::max(
            result.maximumDepth, sourceCandidate.compiledCircuit.cost.multiplicativeDepth);
        if (sourceCandidate.intervalCertificate.proved)
            result.bestApproximationBound = std::min(
                result.bestApproximationBound,
                sourceCandidate.intervalCertificate.approximationErrorUpperBoundDouble);

        auto candidate = sourceCandidate;
        PreparedEvalModPlan plan;
        try {
            candidate.compiledCircuit = compileEvalModPolynomial(
                candidate.polynomial, certifiedProblem,
                sourceCandidate.compiledCircuit.babyStep);
            candidate.satisfiesLevelBudget =
                candidate.compiledCircuit.cost.levelConsumption
                <= certifiedProblem.availableLevels;
            plan = prepareEvalMod(
                adapter, candidate, certifiedProblem, evalModInput);
        } catch (const std::exception& error) {
            plan.status = EvalModCertificationStatus::ScaleScheduleInfeasible;
            plan.detail = error.what();
        }
        if (plan.arithmeticError.rigorous)
            result.bestArithmeticBound = std::min(
                result.bestArithmeticBound, plan.arithmeticError.upperBound);
        if (!firstFailureRecorded && !plan.rigorous) {
            result.firstFailingGate = plan.status;
            firstFailureRecorded = true;
        }
        everyFailureWasResource = everyFailureWasResource
            && (plan.status == EvalModCertificationStatus::InsufficientLevels
                || plan.status == EvalModCertificationStatus::SecurityBudgetExceeded);
        everyFailureWasUnavailableBound = everyFailureWasUnavailableBound
            && (plan.status == EvalModCertificationStatus::RigorousRescaleBoundUnavailable
                || plan.status
                    == EvalModCertificationStatus::RigorousKeySwitchBoundUnavailable
                || plan.status == EvalModCertificationStatus::UnknownOperationBound);
        if (plan.rigorous && sourceCandidate.cost.latencyMs < selectedLatency) {
            selectedLatency = sourceCandidate.cost.latencyMs;
            result.plan = std::move(plan);
        }
    }
    if (!sawCandidate) result.minimumDegree = 0;
    if (result.plan) {
        result.status = EvalModPlanSearchStatus::Certified;
        result.firstFailingGate = EvalModCertificationStatus::Certified;
        result.detail = "selected the lowest predicted-latency certified plan";
        result.proofScope = "certificate applies to the selected prepared plan and exact SEAL context";
        return result;
    }
    result.globalImpossibilityProved = false;
    std::ostringstream scope;
    scope << "implemented families=" << result.familiesSearched.size()
          << ", degree_range=" << result.minimumDegree << ':' << result.maximumDegree
          << ", maximum_depth=" << result.maximumDepth;
    result.proofScope = scope.str();
    if (sawCandidate && everyFailureWasResource) {
        result.status = EvalModPlanSearchStatus::ProfileResourceInfeasible;
        result.detail = "all searched candidates exceed the exact level/security profile";
        result.recommendedRigorousNextPath =
            "reduce multiplicative depth or select a larger security-valid polynomial modulus degree";
    } else if (sawCandidate && everyFailureWasUnavailableBound) {
        result.status = EvalModPlanSearchStatus::RigorousBoundUnavailable;
        result.detail = "all searched candidates require an operation bound not implemented rigorously";
        result.recommendedRigorousNextPath =
            "derive a sampler-specific deterministic or proven probabilistic operation bound";
    } else {
        result.status = EvalModPlanSearchStatus::NoCertifiedPlanInSearchSpace;
        result.detail = "no certified candidate was found in the explicitly reported search space";
        result.recommendedRigorousNextPath =
            "expand target-driven Chebyshev/composite structures and tighten proven operation norms";
    }
    return result;
}

EvalModPreflightResult preflightEvalMod(
    const SealAdapter& adapter,
    const PreparedEvalModPlan& plan,
    const Cipher& evalModInput) {
    if (!plan.rigorous || plan.status != EvalModCertificationStatus::Certified)
        return {plan.status, false, "EvalMod plan is not certified: " + plan.detail};
    if (adapter.contextFingerprint() != plan.contextFingerprint)
        return {EvalModCertificationStatus::ContextMismatch, false,
                "SEAL context fingerprint differs from prepared plan"};
    if (adapter.coeffModulusValues(evalModInput) != plan.dataPrimes)
        return {EvalModCertificationStatus::ContextMismatch, false,
                "active coefficient modulus primes differ from prepared plan"};
    const auto info = adapter.info(evalModInput);
    if (info.chainIndex != plan.inputChainIndex)
        return {EvalModCertificationStatus::InputLevelMismatch, false,
                "EvalMod input chain level differs from prepared plan"};
    if (binary64Bits(info.scale) != plan.inputScaleBinary64Bits)
        return {EvalModCertificationStatus::InputScaleMismatch, false,
                "EvalMod input binary64 scale differs from prepared plan"};
    if (info.ciphertextSize != 2)
        return {EvalModCertificationStatus::CiphertextRepresentationMismatch, false,
                "prepared EvalMod plan requires a two-component CKKS ciphertext"};
    if (plan.circuit.cost.relinearizations && !adapter.hasRelinKeys())
        return {EvalModCertificationStatus::MissingEvaluationKeys, false,
                "prepared EvalMod plan requires relinearization keys"};
    return {EvalModCertificationStatus::Certified, true, "runtime input matches prepared plan"};
}

Cipher applyEvalMod(
    SealAdapter& adapter,
    const PreparedEvalModPlan& plan,
    const Cipher& evalModInput,
    EvalModExecutionTrace* trace) {
    const auto preflight = preflightEvalMod(adapter, plan, evalModInput);
    if (!preflight.compatible) throw std::invalid_argument(preflight.detail);
    return executePreparedEvalMod(
        adapter, plan.circuit, plan.constants, evalModInput, trace);
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
        const EvalModExactModulusContext runtimeContext{
            adapter.coeffModulusValues(encrypted), adapter.specialKeyModulusValue()};
        const auto preparedConstants = prepareEvalModConstants(
            adapter, candidate.compiledCircuit, encrypted);
        const auto certifiedScaleSchedule = buildExactEvalModScaleSchedule(
            candidate.compiledCircuit, runtimeContext, true);
        if (certifiedScaleSchedule.valid) {
            const double normalization = exactScaleUp(
                candidate.compiledCircuit.normalizationGain);
            const double normalizedInputBound =
                static_cast<double>(candidate.domain.integerBound)
                + candidate.domain.normalizedResidualBound;
            candidate.arithmeticCertificate = certifyEvalModDagArithmetic(
                candidate.compiledCircuit, certifiedScaleSchedule, preparedConstants,
                normalizedInputBound / normalization,
                problem.ciphertextModel.normalizedCoeffToSlotErrorAbsBound / normalization);
            candidate.exactScaleSchedule = certifiedScaleSchedule;
            candidate.arithmeticErrorRigorous = candidate.arithmeticCertificate.rigorous;
            candidate.arithmeticErrorCertified = candidate.arithmeticCertificate.rigorous;
            candidate.analyticalArithmeticBound =
                candidate.arithmeticCertificate.outputError.upperBound;
            candidate.analyticalArithmeticBoundRigorous =
                candidate.arithmeticCertificate.outputError.rigorous;
            if (candidate.arithmeticCertificate.rigorous) {
                candidate.polynomialArithmeticError =
                    candidate.arithmeticCertificate.outputError.upperBound;
                candidate.propagationBounds.polynomialArithmetic =
                    candidate.arithmeticCertificate.outputError.upperBound;
                candidate.predictedBootstrapError =
                    propagatedBootstrapError(candidate.propagationBounds);
                candidate.certifiedBootstrapError = candidate.predictedBootstrapError;
                candidate.satisfiesNumericalTarget = candidate.intervalCertified
                    && candidate.predictedBootstrapError <= problem.targetAbsoluteError;
            }
        }
        result.preparedConstantsUsed = preparedConstants.rigorous;
        for (const auto& prepared : preparedConstants.constants)
            result.maxPreparedConstantEncodingError = std::max(
                result.maxPreparedConstantEncodingError,
                prepared.encodingErrorUpperBound);
        EvalModExecutionTrace trace;
        const Cipher output = executePreparedEvalMod(
            adapter, candidate.compiledCircuit, preparedConstants, encrypted,
            &trace, &rawInput, &result.differentialTrace);
        const auto runtimeScaleSchedule = buildExactEvalModScaleSchedule(
            candidate.compiledCircuit, runtimeContext, false);
        result.runtimeScaleBitsMatch = runtimeScaleSchedule.valid;
        if (runtimeScaleSchedule.valid) {
            for (std::size_t index = 0; index < candidate.compiledCircuit.nodes.size(); ++index) {
                if (candidate.compiledCircuit.nodes[index].operation
                    == EvalModOperation::EncodeConstant) continue;
                if (binary64Bits(trace.nodeStates[index].scale)
                    != runtimeScaleSchedule.scaleValues[index].runtimeBinary64Bits) {
                    result.runtimeScaleBitsMatch = false;
                    result.firstRuntimeScaleMismatchNode = index;
                    break;
                }
            }
        }
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
        result.passed = result.executionSucceeded && result.runtimeScaleBitsMatch
            && result.matchesPolynomialReference
            && approximationWithinBudget && result.matchesEvalModTarget;
        if (!result.runtimeScaleBitsMatch) result.failure = "runtime_scale_schedule_mismatch";
        else if (!result.matchesPolynomialReference) result.failure = "implementation_budget_exceeded";
        else if (!approximationWithinBudget) result.failure = "approximation_budget_exceeded";
        else if (!result.matchesEvalModTarget) result.failure = "total_budget_exceeded";
        candidate.measuredBackendError = result.implementationError;
        candidate.backendMeasured = result.executionSucceeded;
        candidate.approximationCertified = candidate.intervalCertified;
        candidate.rigorouslyValidated = candidate.approximationCertified
            && candidate.arithmeticErrorCertified
            && candidate.exactScaleSchedule.valid
            && candidate.exactScaleSchedule.rigorous
            && candidate.satisfiesNumericalTarget;
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
