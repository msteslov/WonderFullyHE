#include "m2424/m2424.hpp"
#include "bootstrap_fixture.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <vector>

namespace {

using namespace m2424;

BootstrapBound deterministicBound(double value, const char* provenance) {
    return {value, BootstrapBoundKind::Deterministic, provenance, {}};
}

double maxAbs(const std::vector<double>& values) {
    double result = 0.0;
    for (double value : values) {
        result = std::max(result, std::abs(value));
    }
    return result;
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
        const CkksProfile profile{N, std::vector<int>(15, 50), std::ldexp(1.0, 49), loans};

        // Research-only fixture: deliberately insecure tiny-N context and a synthetic
        // c=(m,1) input. It is used to exercise the currently executable K=1 bootstrap
        // path and must not be interpreted as production security evidence.
        auto adapter = test::BootstrapFixture::create(profile);
        Bootstrapper bootstrap(N, 2, 2);

        auto requiredKeys = bootstrap.rotationKeys();
        requiredKeys.push_back(0);
        adapter.generateKeys(requiredKeys, true);

        if (adapter.securityLevelBits() != 0) {
            throw std::runtime_error("Financial bootstrap experiment must remain an explicit security-none fixture");
        }

        // Synthetic confidential probabilities of default. In a complete application
        // these would be the output of an encrypted score/sigmoid stage.
        const std::vector<double> pd{
            0.004, 0.008, 0.012, 0.016,
            0.020, 0.024, 0.028, 0.032
        };

        // EAD * LGD * stress multiplier, dollars per unit probability of default.
        // These are public model parameters and therefore stay plaintext.
        const std::vector<double> lossWeights{
            42000.0, 56000.0, 67500.0, 72000.0,
            82500.0, 90000.0, 104000.0, 117000.0
        };

        double expectedLossPlain = 0.0;
        for (std::size_t i = 0; i < loans; ++i) {
            expectedLossPlain += pd[i] * lossWeights[i];
        }

        // Build a low-level source ciphertext for the research fixture. The plaintext
        // itself is ordinary CKKS slot encoding; only the synthetic c=(m,1) wrapper is
        // test-only.
        auto topZero = adapter.encrypt(adapter.encode({0.0}));
        auto bottomZero = adapter.modSwitchToChainIndex(topZero, 0);
        auto sourcePlain = adapter.modSwitchPlainTo(adapter.encode(pd), bottomZero);
        auto input = test::BootstrapFixture::input(adapter, sourcePlain);

        BootstrapRequest request;
        request.target.targetAbsoluteError = 1e-6;

        // Conservative coefficient-domain envelope for the CKKS encoding. For this
        // tiny fixture N * ||slots||_inf * Delta is comfortably below q_src/2.
        const double messageMagnitude = static_cast<double>(N) * maxAbs(pd) * profile.scale;
        request.upstream.messageMagnitude = deterministicBound(
            messageMagnitude,
            "Conservative N*||PD||_inf*Delta coefficient envelope for financial fixture");
        request.upstream.sourceNoiseMagnitude = deterministicBound(
            1.0,
            "Synthetic c=(m,1): source noise is ternary secret with coefficient support one");
        request.upstream.raisedMagnitude = deterministicBound(
            static_cast<double>(N) * (messageMagnitude + 1.0) / profile.scale,
            "Conservative canonical raised-magnitude envelope for financial fixture");
        request.lift = {
            1,
            BootstrapLiftEvidence::TestFixtureAssumption,
            "Explicit K=1 research-fixture assumption; never production evidence"
        };

        const auto plan = bootstrap.prepare(adapter, input, request);
        if (plan.executionReadiness().status != BootstrapCertificationStatus::Certified) {
            throw std::runtime_error(
                plan.executionReadiness().gate + ": " + plan.executionReadiness().provenance);
        }

        const auto refreshed = bootstrap.apply(adapter, input, plan);
        if (!refreshed.output) {
            throw std::runtime_error(
                "Bootstrap produced no output: " + refreshed.trace.result.gate + ": "
                + refreshed.trace.result.provenance);
        }

        // Continue the financial computation *after* refresh. One plaintext-vector
        // multiplication applies EAD*LGD*stress factors. A power-of-two plaintext
        // scale gives enough precision while consuming only one remaining level.
        std::vector<std::complex<double>> weightsComplex;
        weightsComplex.reserve(lossWeights.size());
        for (double value : lossWeights) {
            weightsComplex.emplace_back(value, 0.0);
        }

        const double publicWeightScale = std::ldexp(1.0, 20);
        auto weightPlain = adapter.encodeComplexAtScaleFor(
            weightsComplex, publicWeightScale, *refreshed.output);
        auto weighted = adapter.multiplyPlain(*refreshed.output, weightPlain);
        weighted = adapter.rescaleToNext(weighted);

        auto totalCipher = sumFirstSlots(adapter, weighted, loans);
        const auto decoded = adapter.decode(adapter.decrypt(totalCipher));
        const double expectedLossEncrypted = decoded.front();

        const double absError = std::abs(expectedLossEncrypted - expectedLossPlain);
        const double relativeError = absError / std::max(1.0, std::abs(expectedLossPlain));

        const auto refreshedPd = adapter.decode(adapter.decrypt(*refreshed.output));
        double bootstrapObservedError = 0.0;
        for (std::size_t i = 0; i < loans; ++i) {
            bootstrapObservedError = std::max(
                bootstrapObservedError, std::abs(refreshedPd[i] - pd[i]));
        }

        std::printf("financial_credit_risk_bootstrap_experiment\n");
        std::printf("mode=research_fixture_k1\n");
        std::printf("security_bits=%d\n", adapter.securityLevelBits());
        std::printf("loans=%zu\n", loans);
        std::printf("bootstrap_global_status=%d\n", static_cast<int>(refreshed.trace.result.status));
        std::printf("bootstrap_first_global_gate=%s\n", refreshed.trace.result.gate.c_str());
        std::printf("bootstrap_certified_error=%.12e\n",
                    refreshed.trace.bounds.at("E_boot").upperBound);
        std::printf("bootstrap_observed_pd_error=%.12e\n", bootstrapObservedError);
        std::printf("post_bootstrap_chain_index=%zu\n", adapter.chainIndex(*refreshed.output));
        std::printf("final_chain_index=%zu\n", adapter.chainIndex(totalCipher));
        std::printf("expected_loss_plain=%.12f\n", expectedLossPlain);
        std::printf("expected_loss_encrypted=%.12f\n", expectedLossEncrypted);
        std::printf("expected_loss_abs_error=%.12e\n", absError);
        std::printf("expected_loss_relative_error=%.12e\n", relativeError);
        std::printf("status=%s\n", std::isfinite(expectedLossEncrypted) ? "PASS" : "FAIL");

        return std::isfinite(expectedLossEncrypted) ? 0 : 1;
    } catch (const std::exception& error) {
        std::fprintf(stderr, "financial_credit_risk_bootstrap_experiment: %s\n", error.what());
        return 1;
    }
}
