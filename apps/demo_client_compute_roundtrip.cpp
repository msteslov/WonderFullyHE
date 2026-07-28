#include "m2424/accuracy.hpp"
#include "m2424/abft.hpp"
#include "m2424/linear_transform.hpp"
#include "m2424/profiles.hpp"
#include "m2424/seal_adapter.hpp"

#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <vector>

namespace {

bool cannot_decrypt_without_secret(m2424::SealAdapter& adapter, const m2424::Cipher& cipher) {
    try {
        (void)adapter.decrypt(cipher);
    } catch (const std::runtime_error&) {
        return true;
    }
    return false;
}

} // namespace

int main() {
    const auto profile = m2424::profiles::basic_ckks();
    const std::size_t payload_size = 16;

    auto client = m2424::SealAdapter::create(profile);
    const auto rotationSteps = m2424::sum_slots_rotation_steps(payload_size);
    client.generateKeys(rotationSteps, false);

    std::vector<double> data(client.slotCount(), 0.0);
    for (std::size_t i = 0; i < payload_size; ++i) {
        data[i] = 4.0 + 0.25 * static_cast<double>(i) + std::sin(static_cast<double>(i) / 5.0);
    }

    const auto public_key = client.savePublicKey();
    const auto secret_key = client.saveSecretKey();
    const auto galois_keys = client.saveGaloisKeys();
    const auto encrypted_payload = client.saveCipher(client.encrypt(client.encode(data)));

    auto compute = m2424::SealAdapter::create(profile);
    compute.loadPublicKey(public_key);
    compute.loadGaloisKeys(galois_keys);
    auto encrypted_on_compute = compute.loadCipher(encrypted_payload);

    const bool secret_is_absent = cannot_decrypt_without_secret(compute, encrypted_on_compute);

    for (std::size_t step = 1; step < payload_size; step <<= 1) {
        encrypted_on_compute = compute.add(encrypted_on_compute,
                                           compute.rotate(encrypted_on_compute, static_cast<int>(step)));
    }
    const auto encrypted_result = compute.saveCipher(encrypted_on_compute);

    auto result_reader = m2424::SealAdapter::create(profile);
    result_reader.loadSecretKey(secret_key);
    auto result_cipher = result_reader.loadCipher(encrypted_result);
    const auto result = result_reader.decode(result_reader.decrypt(result_cipher));

    const double expected_sum = m2424::abft::checksum(std::vector<double>(data.begin(), data.begin() + payload_size));
    const double expected_mean = expected_sum / static_cast<double>(payload_size);
    const double encrypted_sum = result.front();
    const double encrypted_mean = encrypted_sum / static_cast<double>(payload_size);
    const auto accuracy = m2424::compare({expected_sum, expected_mean}, {encrypted_sum, encrypted_mean}, 1e-5);

    std::printf("client_compute_roundtrip\n");
    std::printf("payload_size=%zu\n", payload_size);
    std::printf("public_key_bytes=%zu\n", public_key.size());
    std::printf("secret_key_bytes=%zu\n", secret_key.size());
    std::printf("galois_keys_bytes=%zu\n", galois_keys.size());
    std::printf("ciphertext_in_bytes=%zu\n", encrypted_payload.size());
    std::printf("ciphertext_out_bytes=%zu\n", encrypted_result.size());
    std::printf("compute_context_has_secret_key=%s\n", secret_is_absent ? "NO" : "YES");
    std::printf("sum_cpu=%.12e\n", expected_sum);
    std::printf("sum_encrypted=%.12e\n", encrypted_sum);
    std::printf("mean_cpu=%.12e\n", expected_mean);
    std::printf("mean_encrypted=%.12e\n", encrypted_mean);
    std::printf("max_abs_error=%.6e\n", accuracy.max_abs_error);
    std::printf("status=%s\n", secret_is_absent && accuracy.ok ? "PASS" : "FAIL");

    return secret_is_absent && accuracy.ok ? 0 : 1;
}
