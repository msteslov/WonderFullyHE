#include "m2424/bootstrap_contract.hpp"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <stdexcept>

namespace m2424 {
namespace {
using Status = BootstrapCertificationStatus;
constexpr double negativeInfinity = -std::numeric_limits<double>::infinity();
BootstrapContractResult fail(Status status, std::string gate, std::string detail) {
    return {status, std::move(gate), std::move(detail)};
}
bool known(const BootstrapBound& b) {
    if (!std::isfinite(b.upperBound) || b.upperBound < 0 || b.provenance.empty()) return false;
    if (b.kind == BootstrapBoundKind::Deterministic)
        return b.log2FailureProbability == negativeInfinity;
    return b.kind == BootstrapBoundKind::Probabilistic
        && std::isfinite(b.log2FailureProbability) && b.log2FailureProbability <= 0;
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
std::optional<double> bootstrapFailureLog2UpperBound(const std::vector<BootstrapBound>& bounds) {
    double maximum = negativeInfinity;
    std::size_t count = 0;
    for (const auto& b : bounds) {
        if (!known(b)) return std::nullopt;
        if (b.kind == BootstrapBoundKind::Probabilistic) {
            maximum = std::max(maximum, b.log2FailureProbability);
            ++count;
        }
    }
    if (!count) return negativeInfinity;
    unsigned ceilLog2 = 0;
    for (auto n = count - 1; n; n >>= 1) ++ceilLog2;
    if (!ceilLog2) return maximum;
    // Error-free TwoSum detects whether binary64 addition rounded downward.
    const double sum = maximum + ceilLog2;
    const double b = sum - maximum;
    const double residual = (maximum - (sum - b)) + (ceilLog2 - b);
    return std::min(0.0, residual > 0
        ? std::nextafter(sum, std::numeric_limits<double>::infinity()) : sum);
}
BootstrapContractResult validateBootstrapCertificate(const BootstrapTarget& t,
    const BootstrapPlanMetadata& p, const BootstrapCertificate& c) {
    auto result = validateBootstrapTarget(t);
    if (result.status != Status::Certified) return result;
    if (!p.input || !validContext(*p.input))
        return fail(Status::MissingExactModulusContext, "input.context", "Resolve the actual input before certification");
    std::vector<BootstrapBound> bounds;
    for (std::size_t i = 0; i < c.gates.size(); ++i) {
        const auto& gate = c.gates[i];
        if (!gate.verified || gate.provenance.empty() || gate.requiredBounds.empty())
            return fail(gateStatuses[i], gateNames[i], gate.provenance.empty() ? "Required gate evidence unavailable" : gate.provenance);
        for (const auto& b : gate.requiredBounds) {
            if (!known(b)) return fail(Status::RequiredBoundUnavailable, gateNames[i], b.provenance.empty() ? "Required local bound unavailable" : b.provenance);
            bounds.push_back(b);
        }
    }
    if (!known(c.outputError)) return fail(Status::RequiredBoundUnavailable, "final.error", "Output bound unavailable");
    bounds.push_back(c.outputError);
    if (c.outputError.upperBound > t.targetAbsoluteError)
        return fail(Status::ErrorBudgetExceeded, "final.error", c.outputError.provenance);
    if (*bootstrapFailureLog2UpperBound(bounds) > t.maxFailureProbabilityLog2)
        return fail(Status::FailureProbabilityExceeded, "final.failure", "Conservative union bound exceeds budget (v9 section 6.4)");
    if (p.levelsUsed > p.input->raisedPrimes.size() - 1)
        return fail(Status::InsufficientLevels, "levels", "Plan consumes more than available raised levels");
    if (!c.minimumSecurityBits || *c.minimumSecurityBits < t.targetSecurityBits)
        return fail(Status::SecurityBudgetExceeded, "security", "Minimum over all public RLWE families unavailable or insufficient");
    if (!std::isfinite(p.latencyMs) || p.latencyMs < 0
        || (t.maxEvaluationKeyBytes && p.evaluationKeyBytes > t.maxEvaluationKeyBytes)
        || (t.maxLatencyMs && p.latencyMs > t.maxLatencyMs))
        return fail(Status::ResourceBudgetExceeded, "resources", "Invalid or exceeded resource estimate");
    return {Status::Certified, {}, {}};
}
} // namespace m2424
