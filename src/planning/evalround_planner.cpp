#include "m2424/evalround.hpp"
#include <algorithm>
#include <cfenv>
#include <cmath>
#include <set>
#include <stdexcept>

namespace m2424 {
namespace {
unsigned base(EvalRoundRadix radix) {
    const auto b = static_cast<unsigned>(radix);
    if (b != 2 && b != 3) throw std::invalid_argument("EvalRound radix must be 2 or 3");
    return b;
}
bool nonnegative(double x) { return std::isfinite(x) && x >= 0; }
double up(double x) { return std::nextafter(x, std::numeric_limits<double>::infinity()); }
double add(double a, double b) {
    if (a == 0) return b;
    if (b == 0) return a;
    const double sum = a + b, v = sum - a;
    return (a - (sum - v)) + (b - v) > 0 ? up(sum) : sum;
}
double mul(double a, double b) {
    if (a == 0 || b == 0) return 0;
    if (a == 1) return b;
    if (b == 1) return a;
    const double product = a * b;
    int exponent;
    // Scaling a normal finite result by a power of two is exact in binary64.
    if (std::isnormal(product) && (std::frexp(a, &exponent) == 0.5 || std::frexp(b, &exponent) == 0.5))
        return product;
    return up(product);
}
bool roundingSupported() {
#ifdef __FAST_MATH__
    return false;
#else
    return std::fegetround() == FE_TONEAREST;
#endif
}
EvalRoundRejectionReason problemReason(const EvalRoundProblem& p) {
    if (!nonnegative(p.rho) || p.rho >= 0.5) return EvalRoundRejectionReason::DomainViolation;
    if (!roundingSupported() || !std::isfinite(p.requiredIntegerError) || p.requiredIntegerError <= 0
        || p.maxCleaningRoundsPerDigit > 64 || std::isnan(p.maxFailureProbabilityLog2)
        || p.maxFailureProbabilityLog2 > 0) return EvalRoundRejectionReason::InvalidProblem;
    return EvalRoundRejectionReason::None;
}
BootstrapBound deterministic(double value, const std::string& provenance) {
    return {value, BootstrapBoundKind::Deterministic, provenance, {}};
}
bool known(const BootstrapBound& bound, const std::vector<BootstrapFailureEvent>& events) {
    if (!nonnegative(bound.upperBound) || bound.provenance.empty()) return false;
    if (bound.kind == BootstrapBoundKind::Deterministic) return bound.failureEventIds.empty();
    if (bound.kind != BootstrapBoundKind::Probabilistic || bound.failureEventIds.empty()) return false;
    return std::all_of(bound.failureEventIds.begin(), bound.failureEventIds.end(), [&](const auto& id) {
        return std::any_of(events.begin(), events.end(), [&](const auto& e) { return e.id == id; });
    });
}
// Upward binary64 enclosure of 2/sqrt(3), from 1.15470053837925152901829756...
constexpr double ternaryGain = 0x1.279a74590331dp+0;
double contribution(EvalRoundRadix radix, std::size_t digit, double error, double local) {
    double weight = 1;
    for (std::size_t j = 0; j < digit; ++j) weight *= base(radix); // uint32 K: weights < 2^53.
    return mul(weight, add(radix == EvalRoundRadix::Binary ? error : mul(ternaryGain, error), local));
}
}
std::size_t evalRoundDigitCount(std::uint32_t K, EvalRoundRadix radix) {
    const auto b = base(radix);
    std::uint64_t capacity = 1;
    std::size_t count = 0;
    while (capacity < 2ULL * K + 1) { capacity *= b; ++count; }
    return count;
}
std::vector<int> evalRoundIntegerDigits(std::int64_t I, std::uint32_t K, EvalRoundRadix radix) {
    const auto count = evalRoundDigitCount(K, radix);
    if (I < -static_cast<std::int64_t>(K) || I > K) throw std::invalid_argument("integer outside EvalRound domain");
    std::int64_t value = radix == EvalRoundRadix::Binary ? I + K : I;
    std::vector<int> digits;
    for (std::size_t j = 0; j < count; ++j) {
        int digit = static_cast<int>(value % base(radix));
        if (radix == EvalRoundRadix::BalancedTernary) {
            if (digit == 2) digit = -1;
            if (digit == -2) digit = 1;
        }
        digits.push_back(digit);
        value = (value - digit) / base(radix);
    }
    return digits;
}
std::int64_t reconstructEvalRoundInteger(const std::vector<int>& digits, std::uint32_t K, EvalRoundRadix radix) {
    if (digits.size() != evalRoundDigitCount(K, radix)) throw std::invalid_argument("digit count mismatch");
    std::int64_t result = radix == EvalRoundRadix::Binary ? -static_cast<std::int64_t>(K) : 0;
    std::int64_t weight = 1;
    for (int digit : digits) {
        if (digit > 1 || digit < (radix == EvalRoundRadix::Binary ? 0 : -1))
            throw std::invalid_argument("invalid digit");
        result += weight * digit; weight *= base(radix);
    }
    if (result < -static_cast<std::int64_t>(K) || result > K) throw std::invalid_argument("reconstruction outside domain");
    return result;
}
double evalRoundCleaningErrorUpper(EvalRoundRadix radix, double a, double localError) {
    base(radix);
    if (!roundingSupported() || !nonnegative(a) || a > 1 || !nonnegative(localError))
        throw std::invalid_argument("cleaning requires 0<=a<=1 and a known nonnegative local bound");
    return add(mul(radix == EvalRoundRadix::Binary ? 5 : 3, mul(a, a)), localError);
}
double evalRoundReconstructionErrorUpper(EvalRoundRadix radix,
    const std::vector<double>& errors, const std::vector<double>& localErrors) {
    base(radix);
    if (!roundingSupported() || errors.size() != localErrors.size() || errors.size() > 33)
        throw std::invalid_argument("invalid reconstruction bounds");
    double result = 0;
    for (std::size_t j = 0; j < errors.size(); ++j) {
        if (!nonnegative(errors[j]) || !nonnegative(localErrors[j])) throw std::invalid_argument("invalid reconstruction error");
        result = add(result, contribution(radix, j, errors[j], localErrors[j]));
    }
    return result;
}
EvalRoundCandidate makeEvalRoundReferenceCandidate(const EvalRoundProblem& p, EvalRoundRadix radix,
    EvalRoundExtractionMethod method, const EvalRoundCost& cost) {
    if (problemReason(p) != EvalRoundRejectionReason::None) throw std::invalid_argument("invalid EvalRound problem");
    EvalRoundCandidate c;
    c.radix = radix; c.cost = cost;
    c.id = "reference_" + std::to_string(base(radix)) + "_" + std::to_string(static_cast<int>(method));
    c.extraction = {method, "exact mathematical reference extraction", {}, {}, true, p.K, p.rho, ""};
    c.digits.resize(evalRoundDigitCount(p.K, radix));
    std::vector<double> errors(c.digits.size());
    if (method == EvalRoundExtractionMethod::PiecewiseReference) {
        c.extraction.provenance = "v9 4.1/4.2: disjoint closed intervals rho<1/2 have a unique integer center; exact digit target";
    } else if (method == EvalRoundExtractionMethod::BinaryQuadraticK1 && p.K == 1 && radix == EvalRoundRadix::Binary) {
        // J=I+1: b0=1-I^2, b1=(I^2+I)/2 for I=-1,0,1.
        errors = {add(mul(2, p.rho), mul(p.rho, p.rho)),
                  mul(0.5, add(mul(3, p.rho), mul(p.rho, p.rho)))};
        c.extraction.polynomialDegrees = {2, 2};
        c.extraction.provenance = "All K=1 intervals: b0(x)=1-x^2, b1(x)=(x^2+x)/2; expand x=I+xi, |xi|<=rho";
    } else if (method == EvalRoundExtractionMethod::TernaryPhaseReferenceK1 && p.K == 1 && radix == EvalRoundRadix::BalancedTernary) {
        errors = {mul(up(2 * 0x1.921fb54442d19p+1 / 3), p.rho)};
        c.extraction.provenance = "v9 4.2 phase reference exp(2*pi*i*x/3): |exp(i*t)-1|<=|t| on every K=1 interval";
    } else throw std::invalid_argument("unsupported reference extraction/radix/domain");
    for (std::size_t j = 0; j < c.digits.size(); ++j) {
        c.digits[j].extractionError = deterministic(errors[j], c.extraction.provenance);
        c.digits[j].cleaningLocalErrors.assign(p.maxCleaningRoundsPerDigit,
            deterministic(0, "Exact-arithmetic reference cleaner; not a backend arithmetic certificate"));
        c.digits[j].reconstructionLocalError = deterministic(0, "Exact-arithmetic reference reconstruction; binary -K is exact");
    }
    return c;
}
EvalRoundPlan planEvalRoundCandidate(const EvalRoundProblem& p, const EvalRoundCandidate& c) {
    EvalRoundPlan result; result.problem = p; result.candidateId = c.id; result.radix = c.radix;
    result.extraction = c.extraction; result.cost = c.cost; result.failureEvents = c.failureEvents;
    const auto reject = [&](EvalRoundRejectionReason reason, const std::string& detail) {
        result.rejection = reason; result.provenance = detail; return result;
    };
    if (problemReason(p) != EvalRoundRejectionReason::None) return reject(problemReason(p), "Invalid EvalRound input domain, budget or search limit");
    if ((c.radix != EvalRoundRadix::Binary && c.radix != EvalRoundRadix::BalancedTernary)
        || c.id.empty() || c.digits.size() != evalRoundDigitCount(p.K, c.radix)
        || !nonnegative(c.cost.extraction) || !nonnegative(c.cost.cleaningPerDigitIteration) || !nonnegative(c.cost.reconstruction))
        return reject(EvalRoundRejectionReason::InvalidCandidate, "Invalid radix, digit count, ID or configured cost");
    if ((c.extraction.method == EvalRoundExtractionMethod::BinaryQuadraticK1
            && (p.K != 1 || c.radix != EvalRoundRadix::Binary))
        || (c.extraction.method == EvalRoundExtractionMethod::TernaryPhaseReferenceK1
            && (p.K != 1 || c.radix != EvalRoundRadix::BalancedTernary))
        || (!c.extraction.polynomialDegrees.empty() && c.extraction.polynomialDegrees.size() != c.digits.size()))
        return reject(EvalRoundRejectionReason::InvalidCandidate, "Extraction method does not match radix, K or polynomial count");
    switch (c.extraction.method) {
    case EvalRoundExtractionMethod::PiecewiseReference:
    case EvalRoundExtractionMethod::BinaryQuadraticK1:
    case EvalRoundExtractionMethod::TernaryPhaseReferenceK1:
    case EvalRoundExtractionMethod::ExternalPolynomial:
    case EvalRoundExtractionMethod::DigitExtract: break;
    default: return reject(EvalRoundRejectionReason::InvalidCandidate, "Unknown extraction method");
    }
    if (c.extraction.method == EvalRoundExtractionMethod::DigitExtract) {
        const auto epsilon = c.extraction.digitExtractEpsilon;
        if (!(p.rho > 0 && p.rho <= 0.125 && epsilon && std::isfinite(*epsilon) && 2*p.rho < *epsilon && *epsilon <= 0.25))
            return reject(EvalRoundRejectionReason::DigitExtractDomainViolation, "DigitExtract requires 0<rho<=1/8 and 2*rho<epsilon_de<=1/4");
    }
    if (!c.extraction.verified || c.extraction.provenance.empty()
        || c.extraction.certifiedK != p.K || c.extraction.certifiedRho != p.rho) {
        result.status = EvalRoundPlanStatus::Diagnostic;
        return reject(EvalRoundRejectionReason::ExtractionNotCertified, "Extraction must have a whole-domain certificate bound to K and rho");
    }
    const auto failure = bootstrapFailureLog2UpperBound(c.failureEvents);
    if (!failure || *failure > p.maxFailureProbabilityLog2)
        return reject(EvalRoundRejectionReason::FailureProbabilityExceeded, "Primary extraction/arithmetic assumptions violate failure budget");
    result.log2FailureProbabilityUpper = *failure;
    std::vector<std::vector<double>> errorSchedules;
    bool missingBound = false, cleaningDomain = false;
    for (const auto& d : c.digits) {
        if (!known(d.extractionError, c.failureEvents) || !known(d.reconstructionLocalError, c.failureEvents))
            return reject(EvalRoundRejectionReason::BoundUnavailable, "Unknown extraction or reconstruction bound");
        std::vector<double> errors{d.extractionError.upperBound};
        for (std::size_t r = 0; r < p.maxCleaningRoundsPerDigit; ++r) {
            if (errors.back() > 1) { cleaningDomain = true; break; }
            if (r >= d.cleaningLocalErrors.size() || !known(d.cleaningLocalErrors[r], c.failureEvents)) { missingBound = true; break; }
            const double next = evalRoundCleaningErrorUpper(c.radix, errors.back(), d.cleaningLocalErrors[r].upperBound);
            if (!std::isfinite(next)) { cleaningDomain = true; break; }
            errors.push_back(next);
        }
        errorSchedules.push_back(std::move(errors));
    }
    // For each total iteration count retain the least certified weighted error.
    // All future contributions are nonnegative and independent of past choices.
    struct State { double error; std::vector<std::size_t> counts; };
    std::vector<std::optional<State>> states(1, State{0, {}});
    for (std::size_t j = 0; j < c.digits.size(); ++j) {
        std::vector<std::optional<State>> next(states.size() + errorSchedules[j].size() - 1);
        for (std::size_t used = 0; used < states.size(); ++used) if (states[used]) {
            for (std::size_t r = 0; r < errorSchedules[j].size(); ++r) {
                const double error = add(states[used]->error, contribution(c.radix, j, errorSchedules[j][r], c.digits[j].reconstructionLocalError.upperBound));
                if (!next[used+r] || error < next[used+r]->error) {
                    auto counts = states[used]->counts; counts.push_back(r);
                    next[used+r] = State{error, std::move(counts)};
                }
            }
        }
        states = std::move(next);
    }
    for (std::size_t total = 0; total < states.size(); ++total) if (states[total] && states[total]->error <= p.requiredIntegerError) {
        result.totalCleaningIterations = total; result.integerErrorUpper = states[total]->error;
        result.configuredCost = c.cost.extraction + c.cost.reconstruction + c.cost.cleaningPerDigitIteration * total;
        if (!std::isfinite(result.configuredCost)) return reject(EvalRoundRejectionReason::InvalidCandidate, "Configured cost overflow");
        for (std::size_t j = 0; j < c.digits.size(); ++j) {
            const auto r = states[total]->counts[j];
            result.digits.push_back({r, {errorSchedules[j].begin(), errorSchedules[j].begin()+r+1},
                contribution(c.radix, j, errorSchedules[j][r], c.digits[j].reconstructionLocalError.upperBound),
                c.digits[j].extractionError,
                {c.digits[j].cleaningLocalErrors.begin(), c.digits[j].cleaningLocalErrors.begin()+r},
                c.digits[j].reconstructionLocalError});
        }
        result.status = EvalRoundPlanStatus::Certified; result.rejection = EvalRoundRejectionReason::None;
        result.provenance = "v9 4.3/4.5: minimum total digit-cleaning iterations in configured search; conditional reference error certificate only";
        return result;
    }
    return reject(missingBound ? EvalRoundRejectionReason::BoundUnavailable : cleaningDomain ? EvalRoundRejectionReason::CleaningDomainViolation
        : EvalRoundRejectionReason::ErrorBudgetExceeded, "No certified reconstruction within requiredIntegerError and the configured cleaning search");
}
EvalRoundPlanningResult planEvalRound(const EvalRoundProblem& p, const std::vector<EvalRoundCandidate>& candidates) {
    EvalRoundPlanningResult result;
    std::set<std::string> ids;
    for (const auto& c : candidates) {
        if (!ids.insert(c.id).second) { result.selected.reset(); result.rejection = EvalRoundRejectionReason::InvalidCandidate; result.provenance = "Duplicate candidate ID"; return result; }
        result.candidates.push_back(planEvalRoundCandidate(p, c));
        const auto index = result.candidates.size()-1;
        if (result.candidates.back().status == EvalRoundPlanStatus::Certified
            && (!result.selected || result.candidates[index].configuredCost < result.candidates[*result.selected].configuredCost)) result.selected = index;
    }
    result.rejection = result.selected ? EvalRoundRejectionReason::None : candidates.empty() ? EvalRoundRejectionReason::InvalidCandidate : result.candidates.front().rejection;
    result.provenance = result.selected ? "Minimum configured cost among certified reference candidates" : "No certified candidate in supplied search space";
    return result;
}
} // namespace m2424
