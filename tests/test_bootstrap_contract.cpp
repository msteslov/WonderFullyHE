#include "m2424/bootstrap_contract.hpp"
#include <cmath>
#include <cfenv>
#include <cstdio>
#include <stdexcept>

using namespace m2424;
using Status = BootstrapCertificationStatus;
namespace {
void check(bool condition, const char* detail) {
    if (!condition) throw std::runtime_error(detail);
}
BootstrapBound deterministic(double value = 1e-12) {
    return {value, BootstrapBoundKind::Deterministic,
        "synthetic contract fixture, not a backend proof", {}};
}
}
int main() {
    try {
        BootstrapTarget target;
        check(validateBootstrapTarget(target).status == Status::Certified, "default target");
        for (double bad : {0.0, -1.0, std::numeric_limits<double>::infinity(), std::nan("")}) {
            auto t = target; t.targetAbsoluteError = bad;
            check(validateBootstrapTarget(t).status == Status::InvalidInput, "invalid error target");
        }
        auto invalid = target; invalid.targetSecurityBits = 0;
        check(validateBootstrapTarget(invalid).status == Status::InvalidInput, "security target");
        invalid = target; invalid.maxFailureProbabilityLog2 = std::nan("");
        check(validateBootstrapTarget(invalid).status == Status::InvalidInput, "failure target");
        invalid = target; invalid.maxLatencyMs = -1;
        check(validateBootstrapTarget(invalid).status == Status::InvalidInput, "latency target");
        SealAdapter empty;
        Cipher emptyCipher;
        check(resolveBootstrapInput(empty, emptyCipher).result.status == Status::MissingExactModulusContext, "missing context");
        auto adapter = SealAdapter::create({8192, {50, 40, 40, 50}, std::exp2(40), 4096});
        adapter.generateKeys(false, false);
        const auto cipher = adapter.encrypt(adapter.encode({0.25}));
        const auto before = adapter.saveCipher(cipher);
        auto resolved = resolveBootstrapInput(adapter, cipher);
        check(before == adapter.saveCipher(cipher), "resolver must not mutate ciphertext");
        check(resolved.context.has_value(), "resolve context");
        check(resolved.context->sourcePrimes == adapter.coeffModulusValues(cipher), "exact primes");
        check(resolveBootstrapInput(adapter, cipher, resolved.context->scaleBinary64Bits).context.has_value(), "exact scale");
        check(resolveBootstrapInput(adapter, cipher, resolved.context->scaleBinary64Bits + 1).result.status == Status::InputScaleMismatch, "one ulp mismatch");
        const auto lowered = adapter.rescaleToNext(adapter.multiplyPlain(cipher, adapter.encodeScalarFor(1, cipher)));
        const auto low = resolveBootstrapInput(adapter, lowered);
        check(low.context && low.context->sourcePrimes.size() + 1 == resolved.context->sourcePrimes.size(), "active source primes");
        check(low.context->raisedPrimes == resolved.context->raisedPrimes, "raised chain distinct from source");
        BootstrapPlanMetadata plan; plan.input = resolved.context;
        BootstrapCertificate cert;
        check(validateBootstrapCertificate(target, plan, cert).status != Status::Certified, "default certificate");
        for (std::size_t i = 0; i < cert.gates.size(); ++i) {
            cert.gates[i] = {true, "synthetic gate fixture", {}};
            if (bootstrapGateRequiresNumericalBound(static_cast<BootstrapGate>(i)))
                cert.gates[i].requiredBounds = {deterministic()};
        }
        for (auto gate : {BootstrapGate::ScaleSchedule, BootstrapGate::EvaluationKeys, BootstrapGate::Security})
            check(cert.gates[static_cast<std::size_t>(gate)].requiredBounds.empty(), "logical gate has no synthetic bound");
        cert.outputError = deterministic(1e-11); cert.minimumSecurityBits = 128;
        check(validateBootstrapCertificate(target, plan, cert).status == Status::Certified, "complete contract fixture");
        for (std::size_t i = 0; i < cert.gates.size(); ++i) {
            auto broken = cert; broken.gates[i].verified = false;
            auto r = validateBootstrapCertificate(target, plan, broken);
            check(r.status != Status::Certified && !r.gate.empty() && r.provenance == "synthetic gate fixture", "gate provenance");
            broken = cert; broken.gates[i].requiredBounds.clear();
            check((validateBootstrapCertificate(target, plan, broken).status != Status::Certified)
                == bootstrapGateRequiresNumericalBound(static_cast<BootstrapGate>(i)), "omitted bounds follow gate semantics");
            broken = cert; broken.gates[i].requiredBounds.push_back({});
            check(validateBootstrapCertificate(target, plan, broken).status == Status::RequiredBoundUnavailable, "unknown required bound");
        }
        for (double bad : {-1.0, std::numeric_limits<double>::infinity(), std::nan("")}) {
            auto malformedBound = cert; malformedBound.outputError.upperBound = bad;
            check(validateBootstrapCertificate(target, plan, malformedBound).status == Status::RequiredBoundUnavailable, "invalid bound value");
        }
        auto noProvenance = cert; noProvenance.gates[0].provenance.clear();
        check(validateBootstrapCertificate(target, plan, noProvenance).status != Status::Certified, "missing gate provenance");
        noProvenance = cert; noProvenance.outputError.provenance.clear();
        check(validateBootstrapCertificate(target, plan, noProvenance).status == Status::RequiredBoundUnavailable, "missing bound provenance");
        auto broken = cert; broken.gates[1].verified = false; broken.gates[2].verified = false;
        check(validateBootstrapCertificate(target, plan, broken).gate == "sparse_secret", "first failing gate");
        auto missing = plan; missing.input.reset();
        check(validateBootstrapCertificate(target, missing, cert).status == Status::MissingExactModulusContext, "missing plan context");
        broken = cert; broken.outputError.upperBound = 1e-9;
        check(validateBootstrapCertificate(target, plan, broken).status == Status::ErrorBudgetExceeded, "error budget");
        broken = cert; broken.minimumSecurityBits = 127;
        check(validateBootstrapCertificate(target, plan, broken).status == Status::SecurityBudgetExceeded, "security budget");
        auto costly = plan; costly.levelsUsed = 100;
        check(validateBootstrapCertificate(target, costly, cert).status == Status::InsufficientLevels, "levels");
        costly = plan; costly.evaluationKeyBytes = 2;
        auto limited = target; limited.maxEvaluationKeyBytes = 1;
        check(validateBootstrapCertificate(limited, costly, cert).status == Status::ResourceBudgetExceeded, "key memory");
        costly = plan; costly.latencyMs = 2; limited = target; limited.maxLatencyMs = 1;
        check(validateBootstrapCertificate(limited, costly, cert).status == Status::ResourceBudgetExceeded, "latency");
        limited = target; limited.maxLatencyMs = 1;
        check(validateBootstrapCertificate(limited, plan, cert).status == Status::ResourceBudgetExceeded, "unknown latency with limit");
        costly = plan; costly.latencyMs = 0;
        check(validateBootstrapCertificate(limited, costly, cert).status == Status::Certified, "measured zero latency");
        limited = target; limited.maxEvaluationKeyBytes = 1;
        check(validateBootstrapCertificate(limited, plan, cert).status == Status::ResourceBudgetExceeded, "unknown key memory with limit");
        costly = plan; costly.evaluationKeyBytes = 0;
        check(validateBootstrapCertificate(limited, costly, cert).status == Status::Certified, "measured zero key memory");
        for (double bad : {-1.0, std::numeric_limits<double>::infinity(), std::nan("")}) {
            costly = plan; costly.latencyMs = bad;
            check(validateBootstrapCertificate(target, costly, cert).status == Status::ResourceBudgetExceeded, "invalid measured latency without limit");
        }
        check(validateBootstrapCertificate(target, plan, cert).status == Status::Certified, "unknown resources without limits");
        BootstrapFailureEvent first{"noise", -129, "primary noise tail certificate"};
        BootstrapFailureEvent second{"switch", -129, "primary switching tail certificate"};
        check(*bootstrapFailureLog2UpperBound({first, second}) == -128, "union equality");
        check(*bootstrapFailureLog2UpperBound({first, first}) == -129, "duplicate event counted once");
        auto probabilistic = deterministic(); probabilistic.kind = BootstrapBoundKind::Probabilistic;
        probabilistic.failureEventIds = {first.id};
        broken = cert;
        broken.failureEvents = {first, first};
        broken.gates[0].requiredBounds = {probabilistic, probabilistic};
        broken.gates[1].requiredBounds = {probabilistic};
        broken.outputError = probabilistic;
        limited = target; limited.maxFailureProbabilityLog2 = -129;
        check(validateBootstrapCertificate(limited, plan, broken).status == Status::Certified, "derived bounds reuse one primary assumption");
        broken.failureEvents.push_back(second);
        check(validateBootstrapCertificate(limited, plan, broken).status == Status::FailureProbabilityExceeded, "distinct primary event counted");
        check(validateBootstrapCertificate(target, plan, broken).status == Status::Certified, "failure budget equality");
        broken.failureEvents.clear();
        check(validateBootstrapCertificate(target, plan, broken).status == Status::RequiredBoundUnavailable, "unresolved event reference");
        broken = cert; broken.outputError = probabilistic;
        check(validateBootstrapCertificate(target, plan, broken).status == Status::RequiredBoundUnavailable, "output unresolved event reference");
        broken = cert; broken.outputError.failureEventIds = {first.id}; broken.failureEvents = {first};
        check(validateBootstrapCertificate(target, plan, broken).status == Status::RequiredBoundUnavailable, "deterministic bound cannot depend on event");
        broken = cert; probabilistic.failureEventIds.clear(); broken.outputError = probabilistic;
        check(validateBootstrapCertificate(target, plan, broken).status == Status::RequiredBoundUnavailable, "probabilistic bound must name assumptions");
        first.log2FailureProbability = -2001; second.log2FailureProbability = -2001;
        check(*bootstrapFailureLog2UpperBound({first, second}) == -2000, "no underflow");
        check(*bootstrapFailureLog2UpperBound({}) == -std::numeric_limits<double>::infinity(), "deterministic zero failure");
        limited = target; limited.maxFailureProbabilityLog2 = -std::numeric_limits<double>::infinity();
        check(validateBootstrapCertificate(limited, plan, cert).status == Status::Certified, "deterministic-only certificate with zero failure budget");
        const int originalRounding = std::fegetround();
        check(std::fesetround(FE_DOWNWARD) == 0, "set unsupported rounding mode");
        const auto unsupportedRounding = bootstrapFailureLog2UpperBound({first});
        check(std::fesetround(originalRounding) == 0, "restore rounding mode");
        check(!unsupportedRounding, "unsupported rounding must fail closed");
        check(!bootstrapFailureLog2UpperBound({BootstrapFailureEvent{}}), "unknown event");
        auto conflicting = first; conflicting.log2FailureProbability = -127;
        check(!bootstrapFailureLog2UpperBound({first, conflicting}), "conflicting duplicate probability");
        conflicting = first; conflicting.provenance = "different premise";
        check(!bootstrapFailureLog2UpperBound({first, conflicting}), "conflicting duplicate provenance");
        for (double bad : {1.0, std::nan(""), std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()}) {
            auto malformed = first; malformed.log2FailureProbability = bad;
            check(!bootstrapFailureLog2UpperBound({malformed}), "malformed event probability");
        }
        auto malformed = first; malformed.provenance.clear();
        check(!bootstrapFailureLog2UpperBound({malformed}), "event provenance required");
        broken = cert; broken.failureEvents = {malformed};
        check(validateBootstrapCertificate(target, plan, broken).status == Status::FailureProbabilityExceeded, "malformed event rejected by validator");
        first.log2FailureProbability = -128; second.log2FailureProbability = -140;
        const double unequal = *bootstrapFailureLog2UpperBound({first, second});
        // Independent mathematical enclosure of -128 + log2(1 + 2^-12).
        // Decimal endpoints rounded outward from a 100-digit reference.
        check(unequal >= -127.999647822519699 && unequal <= -127.99964782251965,
            "tight outward unequal union");
        check(unequal == *bootstrapFailureLog2UpperBound({second, first}), "event order independence");
        // Independent 100-digit Decimal oracle, rounded to binary64 from above.
        const double fractionalCases[][4] = {
            {-0x1.0000000000000p+7, -0x1.0100000000000p+7, -0x1.fce9edee4f2a8p+6, -0x1.fce9edee4f288p+6},
            {-0x1.0000000000000p-1, -0x1.499999999999ap+3, -0x1.fe57fd354cfd5p-2, -0x1.fe57fd354cfb5p-2},
            {-0x1.0b80000000000p+10, -0x1.0cccccccccccdp+10, -0x1.0b7d858cf6b87p+10, -0x1.0b7d858cf6b67p+10},
            {-0x1.f408000000000p+10, -0x1.f470000000000p+10, -0x1.f3ee1252087e7p+10, -0x1.f3ee1252087c7p+10},
        };
        for (const auto& row : fractionalCases) {
            auto a = first, b = second;
            a.log2FailureProbability = row[0]; b.log2FailureProbability = row[1];
            const double observed = *bootstrapFailureLog2UpperBound({a, b});
            check(observed >= row[2] && observed <= row[3], "fractional log probability oracle enclosure");
        }
        second.log2FailureProbability = -2000;
        check(*bootstrapFailureLog2UpperBound({first, second}) > -128, "tiny positive tail cannot disappear");
        first.log2FailureProbability = -127;
        broken = cert; broken.failureEvents = {first};
        check(validateBootstrapCertificate(target, plan, broken).status == Status::FailureProbabilityExceeded, "failure budget");
        std::puts("[test_bootstrap_contract] PASS; synthetic output bound=1e-11 <= 1e-10; no homomorphic stages executed by contract");
        return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "FAIL: %s\n", e.what()); return 1;
    }
}
