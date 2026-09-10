#include "m2424/bootstrap_contract.hpp"
#include <cmath>
#include <cstdio>
#include <stdexcept>

using namespace m2424;
using Status = BootstrapCertificationStatus;
namespace {
void check(bool condition, const char* detail) {
    if (!condition) throw std::runtime_error(detail);
}
BootstrapBound deterministic(double value = 0) {
    return {value, BootstrapBoundKind::Deterministic,
        -std::numeric_limits<double>::infinity(), "synthetic contract fixture, not a backend proof"};
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
        for (auto& g : cert.gates) g = {true, "synthetic gate fixture", {deterministic()}};
        cert.outputError = deterministic(1e-11); cert.minimumSecurityBits = 128;
        check(validateBootstrapCertificate(target, plan, cert).status == Status::Certified, "complete contract fixture");
        for (std::size_t i = 0; i < cert.gates.size(); ++i) {
            auto broken = cert; broken.gates[i].verified = false;
            auto r = validateBootstrapCertificate(target, plan, broken);
            check(r.status != Status::Certified && !r.gate.empty() && r.provenance == "synthetic gate fixture", "gate provenance");
            broken = cert; broken.gates[i].requiredBounds.clear();
            check(validateBootstrapCertificate(target, plan, broken).status != Status::Certified, "omitted bounds");
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
        auto probabilistic = deterministic(); probabilistic.kind = BootstrapBoundKind::Probabilistic;
        probabilistic.log2FailureProbability = -129;
        check(*bootstrapFailureLog2UpperBound({probabilistic, probabilistic}) == -128, "union equality");
        check(*bootstrapFailureLog2UpperBound({probabilistic, probabilistic, probabilistic}) == -127, "conservative three-event union");
        broken = cert; broken.gates[0].requiredBounds = {probabilistic, probabilistic};
        check(validateBootstrapCertificate(target, plan, broken).status == Status::Certified, "failure budget equality");
        probabilistic.log2FailureProbability = -2001;
        check(*bootstrapFailureLog2UpperBound({probabilistic, probabilistic}) == -2000, "no underflow");
        check(*bootstrapFailureLog2UpperBound({deterministic()}) == -std::numeric_limits<double>::infinity(), "deterministic zero failure");
        check(!bootstrapFailureLog2UpperBound({BootstrapBound{}}), "unknown failure");
        auto malformed = deterministic(); malformed.log2FailureProbability = -128;
        check(!bootstrapFailureLog2UpperBound({malformed}), "malformed deterministic");
        malformed = probabilistic; malformed.log2FailureProbability = std::nan("");
        check(!bootstrapFailureLog2UpperBound({malformed}), "nan failure");
        broken = cert; probabilistic.log2FailureProbability = -127;
        broken.gates[0].requiredBounds = {probabilistic};
        check(validateBootstrapCertificate(target, plan, broken).status == Status::FailureProbabilityExceeded, "failure budget");
        std::puts("[test_bootstrap_contract] PASS; synthetic output bound=1e-11 <= 1e-10; no homomorphic stages executed by contract");
        return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "FAIL: %s\n", e.what()); return 1;
    }
}
