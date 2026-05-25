#include "m2424/accuracy.hpp"
#include "m2424/eval_mod.hpp"
#include "m2424/profiles.hpp"
#include "m2424/seal_adapter.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

namespace {

std::vector<double> head(std::vector<double> values, std::size_t n) {
    values.resize(std::min(values.size(), n));
    return values;
}

} // namespace

int main() {
    constexpr std::size_t payload_size = 64;
    constexpr double sine_tolerance = 1e-12;
    constexpr double cipher_tolerance = 1e-5;

    m2424::EvalModPolynomial eval_mod;
    std::vector<double> input;
    input.reserve(payload_size);
    for (std::size_t i = 0; i < payload_size; ++i) {
        const double t = -1.0 + 2.0 * static_cast<double>(i) / static_cast<double>(payload_size - 1);
        input.push_back(t * m2424::EvalModPolynomial::approximation_bound);
    }

    const auto polynomial_ref = eval_mod.evaluate_plain(input);
    const auto sine_ref = eval_mod.sine_reference(input);
    const auto plain_accuracy = m2424::compare(sine_ref, polynomial_ref, sine_tolerance);

    auto adapter = m2424::SealAdapter::create(m2424::profiles::boot_ckks());
    adapter.keygen(true, false);
    auto encrypted = adapter.encrypt(adapter.encode(input));
    auto evaluated = eval_mod.evaluate(adapter, encrypted);
    auto decoded = head(adapter.decode(adapter.decrypt(evaluated)), payload_size);
    const auto cipher_accuracy = m2424::compare(polynomial_ref, decoded, cipher_tolerance);
    const auto info = adapter.info(evaluated);

    std::printf("metric,value,status\n");
    std::printf("payload_size,%zu,PASS\n", payload_size);
    std::printf("approximation_bound,%.12e,PASS\n", m2424::EvalModPolynomial::approximation_bound);
    std::printf("plain_vs_sine_max_error,%.6e,%s\n", plain_accuracy.max_abs_error, plain_accuracy.ok ? "PASS" : "FAIL");
    std::printf("cipher_vs_plain_max_error,%.6e,%s\n", cipher_accuracy.max_abs_error, cipher_accuracy.ok ? "PASS" : "FAIL");
    std::printf("chain_index_after,%zu,PASS\n", info.chain_index);
    std::printf("coeff_modulus_size_after,%zu,PASS\n", info.coeff_modulus_size);
    std::printf("ciphertext_size_after,%zu,PASS\n", info.ciphertext_size);
    std::printf("serialized_bytes_after,%zu,PASS\n", adapter.serialized_size(evaluated));

    return plain_accuracy.ok && cipher_accuracy.ok ? 0 : 1;
}
