#include "m2424/m2424.hpp"
#include "m2424/evalround_plus_coeff_to_slot.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <stdexcept>

using namespace m2424;
namespace {
void check(bool value, const char* reason) { if (!value) throw std::runtime_error(reason); }
std::uint64_t bits(double value) { std::uint64_t b; std::memcpy(&b, &value, sizeof(b)); return b; }
template<class F> void rejects(F operation, const char* reason) {
    try { operation(); } catch (const std::exception&) { return; }
    throw std::runtime_error(reason);
}
void exactScalarTests() {
    BootstrapInputContext source;
    source.scaleBinary64Bits = bits(1.0);
    source.sourcePrimes = {2};
    const auto half = CoeffToSlotPrefactor::sourceNormalization(source);
    const double tiny = std::numeric_limits<double>::denorm_min();
    check(half.multiplyRounded(tiny) == 0, "tie rounds to even zero");
    check(half.multiplyRounded(3 * tiny) == 2 * tiny, "subnormal tie rounds to even");
    check(half.multiplyRounded(-3 * tiny) == -2 * tiny, "negative subnormal rounding");
    source.sourcePrimes.assign(30, 1152921504606830593ULL);
    source.scaleBinary64Bits = bits(std::ldexp(1.0, 900));
    const auto huge = CoeffToSlotPrefactor::sourceNormalization(source);
    // qSource > binary64 range, but Delta0*operand/qSource remains near 1.
    check(huge.multiplyRounded(std::ldexp(1.0, 900)) > 1.0, "exact product beyond binary64 range");
    source.sourcePrimes.assign(40, 1152921504606830593ULL);
    const auto underflowingAlpha = CoeffToSlotPrefactor::sourceNormalization(source);
    check(underflowingAlpha.multiplyRounded(1) == 0
        && underflowingAlpha.multiplyRounded(std::ldexp(1.0, 900)) > 0,
        "no premature alpha rounding before multiplication");
    source.scaleBinary64Bits = bits(std::ldexp(1.0, 59) * std::sqrt(2.0));
    source.sourcePrimes = {1152921504606830593ULL, 1152921504606748673ULL};
    const auto alpha = CoeffToSlotPrefactor::sourceNormalization(source);
    check(alpha.denominatorFactors() == source.sourcePrimes && alpha.numeratorScaleBits() == source.scaleBinary64Bits,
        "preserve exact factored modulus and dyadic scale");
    // Exact Fraction oracle vectors are appended below.
    const double oracle[][2] = {
        {0x1.0000000000000p+0, 0x1.6a09e667f3e47p-61},
        {-0x1.f9add3739635fp-4, -0x1.6591ad9eb212dp-64},
        {0x1.921fb54442d18p+1, 0x1.1c5831add64d6p-59},
        {0x1.0000000000000p-400, 0x1.6a09e667f3e47p-461},
        {0x1.0000000000000p+500, 0x1.6a09e667f3e47p+439},
    };
    for (const auto& row : oracle) check(alpha.multiplyRounded(row[0]) == row[1], "single-round exact rational oracle");
    source.sourcePrimes.clear();
    rejects([&] { CoeffToSlotPrefactor::sourceNormalization(source); }, "missing modulus");
}
void plainTests() {
    BootstrapInputContext source; source.sourcePrimes = {17, 97}; source.scaleBinary64Bits = bits(1024);
    const auto alpha = CoeffToSlotPrefactor::sourceNormalization(source);
    for (std::size_t n : {8, 16, 32}) for (std::size_t depth : {1, 2, 3}) {
        CoeffToSlotPlan plan(n, depth);
        // Include each basis vector, in particular every second-half phase column.
        for (std::size_t trial = 0; trial <= n; ++trial) {
            std::vector<double> coefficients(n);
            for (std::size_t i = 0; i < n; ++i)
                coefficients[i] = trial == n ? (static_cast<int>((i * 7) % 13) - 6) / 16.0 : (i == trial ? 1 : 0);
            const auto slots = coeffToSlotReference(coefficients);
            const auto hp = plan.applyPlain(slots);
            const auto lp = plan.applyPlain(slots, alpha);
            for (std::size_t i = 0; i < n / 2; ++i) {
                check(std::abs(hp.first[i] - coefficients[i]) <= 1e-12, "plaintext HP first reference");
                check(std::abs(hp.second[i] - coefficients[i + n / 2]) <= 1e-12, "plaintext HP second reference");
                // Small q=1649 is exactly representable: independent reference division.
                check(std::abs(lp.first[i] - coefficients[i] * 1024 / 1649) <= 1e-12, "plaintext LP first reference");
                check(std::abs(lp.second[i] - coefficients[i + n / 2] * 1024 / 1649) <= 1e-12, "plaintext LP second phase reference");
            }
        }
    }
}
void backendTest() {
    constexpr std::size_t n = 16384, slots = n / 2;
    constexpr double diagnosticLimit = 2e-10; // Same as existing CoeffToSlot backend test.
    const CkksProfile profile{n, {60, 60, 60, 60, 60, 60, 60}, std::exp2(59.5), slots};
    EvalRoundPlusCoeffToSlot transform(n);
    const auto requirements = transform.requirements();
    auto steps = requirements.rotationSteps; steps.push_back(0);
    auto adapter = SealAdapter::create(profile);
    adapter.generateKeys(steps, false);
    std::vector<double> values(slots);
    for (std::size_t i = 0; i < slots; ++i) values[i] = (static_cast<int>(i % 17) - 8) / 32.0;
    auto lowered = adapter.encrypt(adapter.encode(values));
    while (adapter.chainIndex(lowered)) lowered = adapter.rescaleToNext(adapter.multiplyPlain(
        lowered, adapter.encodeScalarAtScaleFor(1, std::exp2(60), lowered)));
    const auto source = *resolveBootstrapInput(adapter, lowered).context;
    const auto alpha = CoeffToSlotPrefactor::sourceNormalization(source);
    auto raised = adapter.modRaiseToTop(lowered);
    const auto uOverScale = adapter.decryptRaisedCoefficientsAtRaisedModulus(raised);
    const auto centeredOverScale = adapter.decryptRaisedCoefficientsAtSourceModulus(raised);
    const CoeffToSlotContract contract{"evalround_plus_pr1", slots, n, 59.5, 59.5, 0.25, diagnosticLimit};
    auto prepared = transform.prepare(adapter, raised, source, contract, contract);
    auto wrongSource = source; ++wrongSource.scaleBinary64Bits;
    check(transform.preflight(adapter, raised, wrongSource, contract, contract).status == BootstrapCertificationStatus::InputScaleMismatch,
        "exact one-ulp source mismatch");
    wrongSource = source; --wrongSource.sourcePrimes[0];
    check(transform.preflight(adapter, raised, wrongSource, contract, contract).status == BootstrapCertificationStatus::MissingExactModulusContext,
        "source modulus mismatch");
    wrongSource = source; ++wrongSource.contextFingerprint[0];
    check(transform.preflight(adapter, raised, wrongSource, contract, contract).status == BootstrapCertificationStatus::MissingExactModulusContext,
        "source context fingerprint mismatch");
    auto wrongContract = contract; wrongContract.inputScaleLog2 = 55;
    rejects([&] { transform.prepare(adapter, raised, source, contract, wrongContract); }, "LP contract mismatch must fail preparation");
    auto identityPrepared = transform.plan().prepare(adapter, raised, contract);
    auto lpPrepared = transform.plan().prepare(adapter, raised, contract, alpha);
    check(!transform.plan().isPreparedFor(lpPrepared, adapter, raised, contract), "legacy API rejects LP prefactor");
    check(!transform.plan().isPreparedFor(identityPrepared, adapter, raised, contract, alpha), "LP rejects identity prepared plan");
    check(transform.plan().isPreparedFor(lpPrepared, adapter, raised, contract, alpha), "explicit LP prefactor binding");
    auto differentDepth = EvalRoundPlusCoeffToSlot(n, 3);
    rejects([&] { differentDepth.apply(adapter, raised, prepared); }, "prepared factorization binding");
    // Analytic public upper bound: centered source c0,c1 have coefficients <=q/2,
    // and SEAL ternary secret coefficients have |s_i|<=1. Thus |sigma(u/Delta)|
    // <= N*(N+1)*q/(2*Delta). No secret-key observation enters execution.
    // q has ONE 60-bit prime here; long double is only a test-input bound, not alpha.
    const long double magnitude = static_cast<long double>(n) * (n + 1) * source.sourcePrimes[0]
        / (2 * static_cast<long double>(adapter.scale(lowered)));
    const BootstrapBound inputMagnitude{std::nextafter(static_cast<double>(magnitude), INFINITY),
        BootstrapBoundKind::Deterministic, "centered source components and ternary coefficient support", {}};
    const auto result = transform.apply(adapter, raised, prepared, inputMagnitude);
    const auto hp0 = adapter.decodeComplex(adapter.decrypt(result.hpFirst));
    const auto hp1 = adapter.decodeComplex(adapter.decrypt(result.hpSecond));
    const auto lp0 = adapter.decodeComplex(adapter.decrypt(result.lpFirst));
    const auto lp1 = adapter.decodeComplex(adapter.decrypt(result.lpSecond));
    double hpError = 0, lpError = 0;
    bool integerObserved = false;
    for (std::size_t i = 0; i < n; ++i) {
        // Existing debug oracles return coefficients already divided by Delta0.
        // Recover integer I from the difference of raised/source representatives.
        const double integer = std::round(alpha.multiplyRounded(uOverScale[i] - centeredOverScale[i]));
        integerObserved = integerObserved || integer != 0;
        const double expectedLP = integer + alpha.multiplyRounded(centeredOverScale[i]);
        check(std::abs(alpha.multiplyRounded(uOverScale[i]) - expectedLP) <= 1e-12,
            "debug identity u/q = I + (m+nu_b)/q");
        hpError = std::max(hpError, std::abs((i < slots ? hp0[i] : hp1[i - slots]) - uOverScale[i]));
        lpError = std::max(lpError, std::abs((i < slots ? lp0[i] : lp1[i - slots]) - expectedLP));
    }
    check(integerObserved, "nontrivial ModRaise integer component");
    check(hpError <= diagnosticLimit && lpError <= diagnosticLimit, "CKKS diagnostic regression tolerance");
    for (const auto* trace : {&result.hpTrace, &result.lpTrace}) {
        check(trace->levelsConsumed == transform.plan().depth(), "no standalone alpha level");
        check(trace->rescaleOperations == 2 * transform.plan().depth(), "only factor rescales on two halves");
        check(trace->certificate.status == BootstrapCertificationStatus::RequiredBoundUnavailable && !trace->evidence.verified,
            "no empirical branch certification");
        check(trace->outputError.kind == BootstrapBoundKind::Unknown && !trace->outputError.provenance.empty(), "unknown branch error provenance");
        for (const auto& half : trace->halves) {
            check(half.size() == transform.plan().depth() + 1 && half.back().operation == "conjugate_add", "factor and projection trace");
            for (std::size_t i = 0; i < half.size(); ++i) {
                const auto& factor = half[i];
                check(factor.idealMagnitude.kind == BootstrapBoundKind::Deterministic && factor.idealMagnitude.upperBound > 0,
                    "analytic ideal magnitude recurrence");
                check(factor.propagatedSemanticError.kind == BootstrapBoundKind::Unknown
                    && factor.localAddedError.kind == BootstrapBoundKind::Unknown, "local unknown propagates");
                check(factor.activePrimes.size() == factor.outputState.coeffModulusSize, "exact modulus trace");
                check(factor.outputState.chainIndex + std::min(i + 1, transform.plan().depth()) == adapter.info(raised).chainIndex,
                    "one rescale per factor and none on projection");
            }
        }
    }
    check(result.hpTrace.gate == BootstrapGate::CoeffToSlotHP && result.lpTrace.gate == BootstrapGate::CoeffToSlotLP,
        "only HP and LP gates");
    check(result.hpTrace.prefactor.isIdentity() && result.lpTrace.prefactor == alpha, "trace exact branch scalars");
    check(adapter.chainIndex(result.hpFirst) == adapter.chainIndex(result.lpFirst), "equal HP LP depth");
    // Same state, prepared identity path: new HP must be bitwise the old result.
    CoeffToSlot legacy(n);
    auto legacyRaised = adapter.modRaiseToTop(lowered);
    auto oldResult = legacy.apply(adapter, std::move(legacyRaised), contract, identityPrepared);
    check(adapter.saveCipher(oldResult.slotCipherFirst) == adapter.saveCipher(result.hpFirst)
        && adapter.saveCipher(oldResult.slotCipherSecond) == adapter.saveCipher(result.hpSecond), "legacy HP numerical behavior unchanged");
    // Missing keys are checked before execution; transfer only the public key.
    auto missingKeys = SealAdapter::create(profile);
    missingKeys.loadPublicKey(adapter.savePublicKey());
    const auto missingRaised = missingKeys.modRaiseToTop(missingKeys.loadCipher(adapter.saveCipher(lowered)));
    check(transform.preflight(missingKeys, missingRaised, source, contract, contract).provenance == "missing_rotation_keys", "missing rotations preflight");
    rejects([&] { transform.apply(missingKeys, missingRaised, prepared); }, "missing rotation apply preflight");
    auto noConjugation = SealAdapter::create(profile);
    noConjugation.generateKeys(requirements.rotationSteps, false);
    const auto noConjRaised = noConjugation.modRaiseToTop(noConjugation.loadCipher(adapter.saveCipher(lowered)));
    check(transform.preflight(noConjugation, noConjRaised, source, contract, contract).provenance == "missing_conjugation_key", "missing conjugation preflight");
    rejects([&] { transform.apply(noConjugation, noConjRaised, prepared); }, "missing conjugation apply preflight");
    auto foreign = SealAdapter::create({n, {60,60,60,60,60,59,60}, std::exp2(59.5), slots});
    rejects([&] { transform.apply(foreign, raised, prepared); }, "prepared context fingerprint");
    // One-ulp altered input must be rejected, independently of permissive legacy tolerance.
    auto wrongScale = adapter.normalizeScale(lowered, std::nextafter(adapter.scale(lowered), INFINITY));
    auto wrongRaised = adapter.modRaiseToTop(wrongScale);
    rejects([&] { transform.apply(adapter, wrongRaised, prepared); }, "prepared exact scale binding");
    std::printf("[test_evalround_plus_coeff_to_slot] HP=%.6e LP=%.6e certified_HP=Unknown certified_LP=Unknown "
        "status=RequiredBoundUnavailable levels=%zu rescales_per_branch=%zu PASS\n", hpError, lpError,
        result.hpTrace.levelsConsumed, result.hpTrace.rescaleOperations);
}
}
int main() {
    try { exactScalarTests(); plainTests(); backendTest(); return 0; }
    catch (const std::exception& e) { std::fprintf(stderr, "FAIL: %s\n", e.what()); return 1; }
}
