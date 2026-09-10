#include "m2424/bootstrap_contract.hpp"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <stdexcept>
#include <map>
#include <cfenv>

namespace m2424 {
namespace {
using Status = BootstrapCertificationStatus;
constexpr double negativeInfinity = -std::numeric_limits<double>::infinity();
BootstrapContractResult fail(Status status, std::string gate, std::string detail) {
    return {status, std::move(gate), std::move(detail)};
}
bool known(const BootstrapBound& b, const std::vector<BootstrapFailureEvent>& events) {
    if (!std::isfinite(b.upperBound) || b.upperBound < 0 || b.provenance.empty()) return false;
    if (b.kind == BootstrapBoundKind::Deterministic) return b.failureEventIds.empty();
    if (b.kind != BootstrapBoundKind::Probabilistic || b.failureEventIds.empty()) return false;
    return std::all_of(b.failureEventIds.begin(), b.failureEventIds.end(), [&](const auto& id) {
        return std::any_of(events.begin(), events.end(), [&](const auto& e) { return e.id == id; });
    });
}
// IEEE binary64 round-to-nearest basic operations; no reliance on exp/log libm
// error bounds. TwoSum preserves exact additions (notably boundary budgets).
double up(double x) { return std::nextafter(x, std::numeric_limits<double>::infinity()); }
double down(double x) { return std::nextafter(x, negativeInfinity); }
double addUp(double a, double b) {
    const double sum = a + b;
    const double v = sum - a;
    return (a - (sum - v)) + (b - v) > 0 ? up(sum) : sum;
}
constexpr double ln2Lower = 0x1.62e42fefa39efp-1;
constexpr double ln2Upper = 0x1.62e42fefa39f0p-1;
// x <= 0. Range reduction leaves t in [0, ln(2)]. Taylor through
// degree 24 plus twice term 25 bounds the positive exponential tail.
double exp2Upper(double x) {
    if (x < -1074) return std::numeric_limits<double>::denorm_min();
    const int exponent = static_cast<int>(std::floor(x));
    const double fraction = x - exponent;
    if (fraction == 0) return std::scalbn(1.0, exponent);
    const double t = up(fraction * ln2Upper);
    double term = 1, sum = 1;
    for (int k = 1; k <= 24; ++k) {
        term = up(up(term * t) / k);
        sum = addUp(sum, term);
    }
    term = up(up(term * t) / 25);
    return up(std::scalbn(addUp(sum, up(2 * term)), exponent));
}
// log(m) = 2 sum z^(2k+1)/(2k+1), z=(m-1)/(m+1), 1<=m<2.
// After 32 terms, the tail is <= 2*z^65/(65*(1-z^2)).
double log2Upper(double x) {
    int exponent;
    const double m = 2 * std::frexp(x, &exponent);
    --exponent;
    if (m == 1) return exponent;
    const double z = up((m - 1) / down(m + 1));
    const double square = up(z * z);
    double power = z, sum = 0;
    for (int k = 0; k < 32; ++k) {
        sum = addUp(sum, up(power / (2 * k + 1)));
        power = up(power * square);
    }
    const double tail = up(up(2 * power) / down(65 * down(1 - square)));
    const double logarithm = addUp(up(2 * sum), tail);
    return addUp(exponent, up(logarithm / ln2Lower));
}
bool validContext(const BootstrapInputContext& c) {
    double scale;
    static_assert(sizeof(scale) == sizeof(c.scaleBinary64Bits));
    static_assert(std::numeric_limits<double>::is_iec559);
    std::memcpy(&scale, &c.scaleBinary64Bits, sizeof(scale));
    return !c.sourcePrimes.empty() && c.chainIndex == c.sourcePrimes.size() - 1
        && c.raisedPrimes.size() >= c.sourcePrimes.size()
        && std::equal(c.sourcePrimes.begin(), c.sourcePrimes.end(), c.raisedPrimes.begin())
        && std::all_of(c.raisedPrimes.begin(), c.raisedPrimes.end(), [](auto x) { return x > 1; })
        && c.specialPrime > 1 && std::isfinite(scale) && scale > 0
        && std::any_of(c.contextFingerprint.begin(), c.contextFingerprint.end(), [](auto x) { return x != 0; })
        && std::all_of(c.sourcePrimes.begin(), c.sourcePrimes.end(), [](auto x) { return x > 1; });
}
constexpr const char* gateNames[] = {"input", "sparse_secret", "key_switch", "domain",
    "extraction", "arithmetic", "coeff_to_slot_hp", "coeff_to_slot_lp", "reconstruction",
    "combination", "slot_to_coeff", "scale_schedule", "headroom", "evaluation_keys", "security"};
constexpr Status gateStatuses[] = {Status::InvalidInput, Status::SparseSecretCertificateUnavailable,
    Status::KeySwitchBoundUnavailable, Status::DomainViolation, Status::ExtractionNotCertified,
    Status::RequiredBoundUnavailable, Status::RequiredBoundUnavailable, Status::RequiredBoundUnavailable,
    Status::CleaningBoundUnavailable, Status::RequiredBoundUnavailable, Status::SlotToCoeffGainUnavailable,
    Status::ScaleScheduleInfeasible, Status::HeadroomViolation, Status::MissingEvaluationKeys,
    Status::SecurityBudgetExceeded};
}
BootstrapContractResult validateBootstrapTarget(const BootstrapTarget& t) {
    if (!std::isfinite(t.targetAbsoluteError) || t.targetAbsoluteError <= 0
        || t.targetSecurityBits <= 0 || std::isnan(t.maxFailureProbabilityLog2)
        || t.maxFailureProbabilityLog2 > 0 || !std::isfinite(t.maxLatencyMs) || t.maxLatencyMs < 0)
        return fail(Status::InvalidInput, "target", "Expected positive error/security targets and valid failure/resource budgets");
    return {Status::Certified, {}, {}};
}
BootstrapInputResolution resolveBootstrapInput(const SealAdapter& adapter, const Cipher& input,
                                               std::optional<std::uint64_t> expected) {
    BootstrapInputContext c;
    try {
        c.sourcePrimes = adapter.coeffModulusValues(input);
        c.raisedPrimes = adapter.dataModulusValues();
        c.specialPrime = adapter.specialKeyModulusValue();
        c.contextFingerprint = adapter.contextFingerprint();
        c.chainIndex = adapter.chainIndex(input);
        const auto scale = adapter.scale(input);
        std::memcpy(&c.scaleBinary64Bits, &scale, sizeof(scale));
    } catch (const std::exception& e) {
        return {fail(Status::MissingExactModulusContext, "input.context", e.what()), std::nullopt};
    }
    if (!validContext(c))
        return {fail(Status::MissingExactModulusContext, "input.context", "Invalid exact source context"), std::nullopt};
    if (expected && *expected != c.scaleBinary64Bits)
        return {fail(Status::InputScaleMismatch, "input.scale", "Exact binary64 scale differs from contract"), std::nullopt};
    return {{Status::Certified, {}, {}}, c};
}
bool bootstrapGateRequiresNumericalBound(BootstrapGate gate) {
    switch (gate) {
    case BootstrapGate::ScaleSchedule:
    case BootstrapGate::EvaluationKeys:
    case BootstrapGate::Security: // Dedicated minimumSecurityBits field below.
        return false;
    default:
        return true;
    }
}
std::optional<double> bootstrapFailureLog2UpperBound(const std::vector<BootstrapFailureEvent>& events) {
#ifdef __FAST_MATH__
    return std::nullopt; // Reassociation invalidates the outward arithmetic proof.
#endif
    if (std::fegetround() != FE_TONEAREST) return std::nullopt;
    std::map<std::string, BootstrapFailureEvent> unique;
    double maximum = negativeInfinity;
    for (const auto& event : events) {
        if (event.id.empty() || event.provenance.empty()
            || !std::isfinite(event.log2FailureProbability) || event.log2FailureProbability > 0)
            return std::nullopt;
        const auto inserted = unique.emplace(event.id, event);
        if (!inserted.second && (inserted.first->second.log2FailureProbability != event.log2FailureProbability
            || inserted.first->second.provenance != event.provenance)) return std::nullopt;
        maximum = std::max(maximum, event.log2FailureProbability);
    }
    if (unique.empty()) return negativeInfinity;
    double sum = 0;
    for (const auto& entry : unique)
        sum = addUp(sum, exp2Upper(addUp(entry.second.log2FailureProbability, -maximum)));
    return std::min(0.0, addUp(maximum, log2Upper(sum)));
}
BootstrapContractResult validateBootstrapCertificate(const BootstrapTarget& t,
    const BootstrapPlanMetadata& p, const BootstrapCertificate& c) {
    auto result = validateBootstrapTarget(t);
    if (result.status != Status::Certified) return result;
    if (!p.input || !validContext(*p.input))
        return fail(Status::MissingExactModulusContext, "input.context", "Resolve the actual input before certification");

    for (std::size_t i = 0; i < c.gates.size(); ++i) {
        const auto& gate = c.gates[i];
        if (!gate.verified || gate.provenance.empty() || (bootstrapGateRequiresNumericalBound(static_cast<BootstrapGate>(i))
            && gate.requiredBounds.empty()))
            return fail(gateStatuses[i], gateNames[i], gate.provenance.empty() ? "Required gate evidence unavailable" : gate.provenance);
        for (const auto& b : gate.requiredBounds) {
            if (!known(b, c.failureEvents)) return fail(Status::RequiredBoundUnavailable, gateNames[i], b.provenance.empty() ? "Required local bound unavailable" : b.provenance);
        }
    }
    if (!known(c.outputError, c.failureEvents)) return fail(Status::RequiredBoundUnavailable, "final.error", "Output bound unavailable");
    if (c.outputError.upperBound > t.targetAbsoluteError)
        return fail(Status::ErrorBudgetExceeded, "final.error", c.outputError.provenance);
    const auto failure = bootstrapFailureLog2UpperBound(c.failureEvents);
    if (!failure)
        return fail(Status::FailureProbabilityExceeded, "final.failure", "Invalid primary events or unsupported floating-point rounding mode");
    if (*failure > t.maxFailureProbabilityLog2)
        return fail(Status::FailureProbabilityExceeded, "final.failure", "Conservative union bound exceeds budget (v9 section 6.4)");
    if (p.levelsUsed > p.input->raisedPrimes.size() - 1)
        return fail(Status::InsufficientLevels, "levels", "Plan consumes more than available raised levels");
    if (!c.minimumSecurityBits || *c.minimumSecurityBits < t.targetSecurityBits)
        return fail(Status::SecurityBudgetExceeded, "security", "Minimum over all public RLWE families unavailable or insufficient");
    if ((p.latencyMs && (!std::isfinite(*p.latencyMs) || *p.latencyMs < 0))
        || (t.maxEvaluationKeyBytes && (!p.evaluationKeyBytes || *p.evaluationKeyBytes > t.maxEvaluationKeyBytes))
        || (t.maxLatencyMs && (!p.latencyMs || *p.latencyMs > t.maxLatencyMs)))
        return fail(Status::ResourceBudgetExceeded, "resources", "Invalid or exceeded resource estimate");
    return {Status::Certified, {}, {}};
}
} // namespace m2424
