#include "m2424/m2424.hpp"
#include "m2424/profiles.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdio>
#include <stdexcept>
#include <vector>

namespace {

using namespace m2424;
using Clock = std::chrono::steady_clock;

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
        constexpr std::size_t loans = 8;
        const std::vector<double> pd{
            0.004, 0.008, 0.012, 0.016,
            0.020, 0.024, 0.028, 0.032
        };
        const std::vector<double> lossWeights{
            42000.0, 56000.0, 67500.0, 72000.0,
            82500.0, 90000.0, 104000.0, 117000.0
        };

        double expectedLossPlain = 0.0;
        for (std::size_t i = 0; i < loans; ++i) {
            expectedLossPlain += pd[i] * lossWeights[i];
        }

        const auto contextStart = Clock::now();
        auto adapter = SealAdapter::create(profiles::basic_ckks());
        const auto contextFinish = Clock::now();

        const auto keygenStart = Clock::now();
        adapter.generateKeys(std::vector<int>{1, 2, 4}, true);
        const auto keygenFinish = Clock::now();

        if (adapter.securityLevelBits() < 128) {
            throw std::runtime_error("basic_ckks did not validate at the expected security level");
        }

        std::vector<double> pdSlots(adapter.slotCount(), 0.0);
        std::vector<double> weightSlots(adapter.slotCount(), 0.0);
        for (std::size_t i = 0; i < loans; ++i) {
            pdSlots[i] = pd[i];
            weightSlots[i] = lossWeights[i];
        }

        const auto encryptStart = Clock::now();
        auto encryptedPd = adapter.encrypt(adapter.encode(pdSlots));
        const auto encryptFinish = Clock::now();
        const auto initialInfo = adapter.info(encryptedPd);

        const auto financialStart = Clock::now();
        auto weightPlain = adapter.encodeFor(weightSlots, encryptedPd);
        auto weighted = adapter.multiplyPlain(encryptedPd, weightPlain);
        weighted = adapter.rescaleToNext(weighted);
        auto totalCipher = sumFirstSlots(adapter, weighted, loans);
        const auto financialFinish = Clock::now();

        const auto decryptStart = Clock::now();
        const auto decoded = adapter.decode(adapter.decrypt(totalCipher));
        const auto decryptFinish = Clock::now();
        const double expectedLossEncrypted = decoded.front();

        const double absError = std::abs(expectedLossEncrypted - expectedLossPlain);
        const double relativeError = absError / std::max(1.0, std::abs(expectedLossPlain));

        std::printf("financial_credit_risk_ckks_experiment\n");
        std::printf("profile=basic_ckks\n");
        std::printf("security_bits=%d\n", adapter.securityLevelBits());
        std::printf("poly_modulus_degree=%zu\n", profiles::basic_ckks().polyModulusDegree);
        std::printf("loans=%zu\n", loans);
        std::printf("context_ms=%.6f\n", elapsedMs(contextStart, contextFinish));
        std::printf("keygen_ms=%.6f\n", elapsedMs(keygenStart, keygenFinish));
        std::printf("encrypt_ms=%.6f\n", elapsedMs(encryptStart, encryptFinish));
        std::printf("financial_stage_ms=%.6f\n", elapsedMs(financialStart, financialFinish));
        std::printf("decrypt_ms=%.6f\n", elapsedMs(decryptStart, decryptFinish));
        std::printf("initial_chain_index=%zu\n", initialInfo.chainIndex);
        std::printf("final_chain_index=%zu\n", adapter.chainIndex(totalCipher));
        std::printf("ciphertext_bytes=%zu\n", adapter.serializedSize(encryptedPd));
        std::printf("public_key_bytes=%zu\n", adapter.publicKeySize());
        std::printf("relin_key_bytes=%zu\n", adapter.relinKeysSize());
        std::printf("galois_key_bytes=%zu\n", adapter.galoisKeysSize());
        std::printf("expected_loss_plain=%.12f\n", expectedLossPlain);
        std::printf("expected_loss_encrypted=%.12f\n", expectedLossEncrypted);
        std::printf("expected_loss_abs_error=%.12e\n", absError);
        std::printf("expected_loss_relative_error=%.12e\n", relativeError);

        const bool ok = std::isfinite(expectedLossEncrypted) && relativeError <= 1e-6;
        std::printf("status=%s\n", ok ? "PASS" : "FAIL");
        return ok ? 0 : 1;
    } catch (const std::exception& error) {
        std::fprintf(stderr, "financial_credit_risk_ckks_experiment: %s\n", error.what());
        return 1;
    }
}
