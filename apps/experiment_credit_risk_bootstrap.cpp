#include "m2424/m2424.hpp"
#include "m2424/canonical_embedding_reference.hpp"
#include "bootstrap_fixture.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include <gmpxx.h>

namespace {

using namespace m2424;
using Clock = std::chrono::steady_clock;

BootstrapBound deterministicBound(double value, const char* provenance) {
    return {value, BootstrapBoundKind::Deterministic, provenance, {}};
}

double elapsedMs(Clock::time_point start, Clock::time_point finish) {
    return std::chrono::duration<double, std::milli>(finish - start).count();
}

Cipher sumFirstSlots(SealAdapter& adapter, const Cipher& input, std::size_t count) {
    auto acc = input;
    for (std::size_t step = 1; step < count; step <<= 1) {
        acc = adapter.add(acc, adapter.rotate(acc, static_cast<int>(step)));
    }
    return acc;
}

} // namespace

int main() {
    try {
        constexpr std::size_t N = 16;
        constexpr std::size_t loans = N / 2;
        constexpr double inputNormalization = 4096.0;
        constexpr double publicWeightScale = std::ldexp(1.0, 20);
        const CkksProfile profile{N, std::vector<int>(15, 50), std::ldexp(1.0, 49), loans};

        auto adapter = test::BootstrapFixture::create(profile);
        Bootstrapper bootstrap(N, 2, 2);

        auto requiredKeys = bootstrap.rotationKeys();
        requiredKeys.insert(requiredKeys.end(), {0, 1, 2, 4});
        std::sort(requiredKeys.begin(), requiredKeys.end());
        requiredKeys.erase(std::unique(requiredKeys.begin(), requiredKeys.end()), requiredKeys.end());
        adapter.generateKeys(requiredKeys, true);

        if (adapter.securityLevelBits() != 0) {
            throw std::runtime_error("Financial bootstrap experiment must remain an explicit security-none fixture");
        }

        const std::vector<double> pd{
            0.004, 0.008, 0.012, 0.016,
            0.020, 0.024, 0.028, 0.032
        };
        const std::vector<double> lossWeights{
            42000.0, 56000.0, 67500.0, 72000.0,
            82500.0, 90000.0, 104000.0, 117000.0
        };

        ComplexVector normalizedPd;
        normalizedPd.reserve(loans);
        std::vector<double> normalizedLossWeights;
        normalizedLossWeights.reserve(loans);
        std::vector<std::complex<double>> weightsComplex;
        weightsComplex.reserve(loans);
        for (std::size_t i = 0; i < loans; ++i) {
            normalizedPd.emplace_back(pd[i] / inputNormalization, 0.0);
            const double normalizedWeight = lossWeights[i] * inputNormalization;
            normalizedLossWeights.push_back(normalizedWeight);
            weightsComplex.emplace_back(normalizedWeight, 0.0);
        }

        // Construct the same exact coefficient-domain fixture used by the passing
        // K=1 end-to-end integration path. The inverse reference embedding maps the
        // desired financial slot vector to real polynomial coefficients; rounding at
        // Delta_0 produces exact integer coefficients for RNS encoding.
        const auto normalizedCoefficients = slotToCoeffReference(normalizedPd);
        std::vector<double> integerCoefficients(N, 0.0);
        double messageMagnitude = 0.0;
        for (std::size_t i = 0; i < N; ++i) {
            integerCoefficients[i] = std::nearbyint(normalizedCoefficients[i] * profile.scale);
            messageMagnitude = std::max(messageMagnitude, std::abs(integerCoefficients[i]));
        }

        std::vector<double> roundedNormalizedCoefficients(N, 0.0);
        for (std::size_t i = 0; i < N; ++i) {
            roundedNormalizedCoefficients[i] = integerCoefficients[i] / profile.scale;
        }
        const auto bootstrapInputReference = coeffToSlotReference(roundedNormalizedCoefficients);

        double sourceEncodingError = 0.0;
        double expectedLossTarget = 0.0;
        double expectedLossPlain = 0.0;
        for (std::size_t i = 0; i < loans; ++i) {
            sourceEncodingError = std::max(
                sourceEncodingError, std::abs(bootstrapInputReference[i] - normalizedPd[i]));
            expectedLossTarget += pd[i] * lossWeights[i];
            expectedLossPlain += bootstrapInputReference[i].real() * normalizedLossWeights[i];
        }

        std::vector<std::uint64_t> residues;
        auto primes = adapter.dataModulusValues();
        primes.push_back(adapter.specialKeyModulusValue());
        residues.reserve(primes.size() * N);
        for (auto prime : primes) {
            const mpz_class modulus(std::to_string(prime));
            for (double coefficient : integerCoefficients) {
                const auto integer = static_cast<long>(coefficient);
                mpz_class value(integer), reduced;
                mpz_mod(reduced.get_mpz_t(), value.get_mpz_t(), modulus.get_mpz_t());
                residues.push_back(std::stoull(reduced.get_str()));
            }
        }

        auto top = adapter.encrypt(adapter.encode({0.0}));
        auto bottom = adapter.modSwitchToChainIndex(top, 0);
        auto encoded = adapter.encodePolynomialRnsAtKeyScale(residues, profile.scale);
        auto sourcePlain = adapter.modSwitchPlainTo(encoded, bottom);
        auto input = test::BootstrapFixture::input(adapter, sourcePlain);

        // Baseline: the application step requires a rescale, but the source ciphertext
        // is already at the last chain level. This check records whether the same
        // financial continuation is executable without refreshing the ciphertext.
        bool noBootstrapContinuation = false;
        std::string noBootstrapFailure;
        try {
            auto baselinePlain = adapter.encodeComplexAtScaleFor(
                weightsComplex, publicWeightScale, input);
            auto baselineWeighted = adapter.multiplyPlain(input, baselinePlain);
            baselineWeighted = adapter.rescaleToNext(baselineWeighted);
            noBootstrapContinuation = true;
        } catch (const std::exception& error) {
            noBootstrapFailure = error.what();
        }

        BootstrapRequest request;
        request.target.targetAbsoluteError = 1e-6;
        request.upstream.messageMagnitude = deterministicBound(
            messageMagnitude,
            "Exact maximum rounded coefficient magnitude of normalized financial fixture");
        request.upstream.sourceNoiseMagnitude = deterministicBound(
            1.0,
            "Synthetic c=(m,1): source noise is ternary secret with coefficient support one");
        request.upstream.raisedMagnitude = deterministicBound(
            static_cast<double>(N) * (messageMagnitude + 1.0) / profile.scale,
            "N*(M_m+1)/Delta canonical raised-magnitude bound for exact fixture coefficients");
        request.lift = {
            1,
            BootstrapLiftEvidence::TestFixtureAssumption,
            "Explicit K=1 integration-fixture assumption"
        };

        const auto prepareStart = Clock::now();
        const auto plan = bootstrap.prepare(adapter, input, request);
        const auto prepareFinish = Clock::now();
        if (plan.executionReadiness().status != BootstrapCertificationStatus::Certified) {
            throw std::runtime_error(
                plan.executionReadiness().gate + ": " + plan.executionReadiness().provenance);
        }

        const auto bootstrapStart = Clock::now();
        const auto refreshed = bootstrap.apply(adapter, input, plan);
        const auto bootstrapFinish = Clock::now();
        if (!refreshed.output) {
            throw std::runtime_error(
                "Bootstrap produced no output: " + refreshed.trace.result.gate + ": "
                + refreshed.trace.result.provenance);
        }

        const auto refreshedPd = adapter.decodeComplex(adapter.decrypt(*refreshed.output));
        double bootstrapObservedError = 0.0;
        for (std::size_t i = 0; i < loans; ++i) {
            bootstrapObservedError = std::max(
                bootstrapObservedError, std::abs(refreshedPd[i] - bootstrapInputReference[i]));
        }

        const auto financialStart = Clock::now();
        auto weightPlain = adapter.encodeComplexAtScaleFor(
            weightsComplex, publicWeightScale, *refreshed.output);
        auto weighted = adapter.multiplyPlain(*refreshed.output, weightPlain);
        weighted = adapter.rescaleToNext(weighted);

        const auto weightedSlots = adapter.decodeComplex(adapter.decrypt(weighted));
        auto totalCipher = sumFirstSlots(adapter, weighted, loans);
        const auto financialFinish = Clock::now();

        const auto decoded = adapter.decodeComplex(adapter.decrypt(totalCipher));
        const double expectedLossEncrypted = decoded.front().real();
        const double absError = std::abs(expectedLossEncrypted - expectedLossPlain);
        const double relativeError = absError / std::max(1.0, std::abs(expectedLossPlain));

        std::printf("financial_credit_risk_bootstrap_experiment\n");
        std::printf("mode=executable_k1_fixture\n");
        std::printf("security_bits=%d\n", adapter.securityLevelBits());
        std::printf("loans=%zu\n", loans);
        std::printf("input_normalization=%.0f\n", inputNormalization);
        std::printf("source_encoding_max_error=%.12e\n", sourceEncodingError);
        std::printf("message_magnitude_bound=%.12e\n", messageMagnitude);
        std::printf("raised_magnitude_bound=%.12e\n", request.upstream.raisedMagnitude.upperBound);
        std::printf("without_bootstrap_continuation=%s\n", noBootstrapContinuation ? "PASS" : "BLOCKED");
        if (!noBootstrapFailure.empty()) {
            std::printf("without_bootstrap_reason=%s\n", noBootstrapFailure.c_str());
        }
        std::printf("prepare_ms=%.6f\n", elapsedMs(prepareStart, prepareFinish));
        std::printf("bootstrap_ms=%.6f\n", elapsedMs(bootstrapStart, bootstrapFinish));
        std::printf("financial_stage_ms=%.6f\n", elapsedMs(financialStart, financialFinish));
        std::printf("bootstrap_global_status=%d\n", static_cast<int>(refreshed.trace.result.status));
        std::printf("bootstrap_first_global_gate=%s\n", refreshed.trace.result.gate.c_str());
        std::printf("bootstrap_certified_error=%.12e\n",
                    refreshed.trace.bounds.at("E_boot").upperBound);
        std::printf("bootstrap_observed_pd_error=%.12e\n", bootstrapObservedError);
        std::printf("post_bootstrap_chain_index=%zu\n", adapter.chainIndex(*refreshed.output));
        std::printf("final_chain_index=%zu\n", adapter.chainIndex(totalCipher));
        for (std::size_t i = 0; i < loans; ++i) {
            const double targetContribution = pd[i] * lossWeights[i];
            const double plainContribution = bootstrapInputReference[i].real() * normalizedLossWeights[i];
            std::printf("loan_%zu_target_el=%.12f\n", i + 1, targetContribution);
            std::printf("loan_%zu_plain_el=%.12f\n", i + 1, plainContribution);
            std::printf("loan_%zu_ckks_el=%.12f\n", i + 1, weightedSlots[i].real());
        }
        std::printf("expected_loss_target=%.12f\n", expectedLossTarget);
        std::printf("expected_loss_plain=%.12f\n", expectedLossPlain);
        std::printf("expected_loss_encrypted=%.12f\n", expectedLossEncrypted);
        std::printf("expected_loss_abs_error=%.12e\n", absError);
        std::printf("expected_loss_relative_error=%.12e\n", relativeError);

        const bool ok = !noBootstrapContinuation
            && sourceEncodingError <= 1e-12
            && bootstrapObservedError <= request.target.targetAbsoluteError
            && std::isfinite(expectedLossEncrypted)
            && relativeError <= 1e-4;
        std::printf("status=%s\n", ok ? "PASS" : "FAIL");
        return ok ? 0 : 1;
    } catch (const std::exception& error) {
        std::fprintf(stderr, "financial_credit_risk_bootstrap_experiment: %s\n", error.what());
        return 1;
    }
}
