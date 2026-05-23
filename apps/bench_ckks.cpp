#include "m2424/accuracy.hpp"
#include "m2424/seal_adapter.hpp"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <functional>
#include <algorithm>
#include <string>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

double elapsed_ms(const std::function<void()>& fn) {
    const auto start = Clock::now();
    fn();
    const auto finish = Clock::now();
    return std::chrono::duration<double, std::milli>(finish - start).count();
}

void print_row(const char* operation, std::size_t poly_degree, std::size_t payload_size,
               double time_ms, double max_error, double mean_error, std::size_t bytes) {
    std::printf("%s,%zu,%zu,%.6f,%.6e,%.6e,%zu\n",
                operation, poly_degree, payload_size, time_ms, max_error, mean_error, bytes);
}

std::vector<double> payload_head(const std::vector<double>& values, std::size_t payload_size) {
    return std::vector<double>(values.begin(), values.begin() + std::min(values.size(), payload_size));
}

} // namespace

int main() {
    const std::size_t poly_degree = 8192;
    const std::size_t payload_size = 64;
    m2424::CkksProfile profile{poly_degree, {60, 40, 40, 60}, std::pow(2.0, 40), poly_degree / 2};

    auto adapter = m2424::SealAdapter::create(profile);
    adapter.keygen(true, true);

    std::vector<double> input;
    input.reserve(payload_size);
    for (std::size_t i = 0; i < payload_size; ++i) {
        input.push_back(std::sin(static_cast<double>(i) / 10.0));
    }

    m2424::Plain plain;
    m2424::Cipher cipher;
    double encode_ms = elapsed_ms([&] { plain = adapter.encode(input); });
    double encrypt_ms = elapsed_ms([&] { cipher = adapter.encrypt(plain); });

    m2424::Cipher added;
    double add_ms = elapsed_ms([&] { added = adapter.add(cipher, cipher); });
    auto add_decoded = payload_head(adapter.decode(adapter.decrypt(added)), payload_size);
    std::vector<double> add_ref;
    add_ref.reserve(payload_size);
    for (double value : input) add_ref.push_back(value + value);
    auto add_accuracy = m2424::compare(add_ref, add_decoded, 1e-5);

    m2424::Cipher multiplied;
    double mul_ms = elapsed_ms([&] { multiplied = adapter.mul_relin_rescale(cipher, cipher); });
    auto mul_decoded = payload_head(adapter.decode(adapter.decrypt(multiplied)), payload_size);
    std::vector<double> mul_ref;
    mul_ref.reserve(payload_size);
    for (double value : input) mul_ref.push_back(value * value);
    auto mul_accuracy = m2424::compare(mul_ref, mul_decoded, 1e-5);

    m2424::Cipher rotated;
    double rotate_ms = elapsed_ms([&] { rotated = adapter.rotate(cipher, 1); });
    auto rotated_decoded = payload_head(adapter.decode(adapter.decrypt(rotated)), payload_size);
    std::vector<double> rotate_ref;
    rotate_ref.reserve(payload_size);
    for (std::size_t i = 0; i < payload_size; ++i) {
        rotate_ref.push_back(input[(i + 1) % payload_size]);
    }
    auto rotate_accuracy = m2424::compare(rotate_ref, rotated_decoded, 1e-5);

    m2424::Plain decrypted;
    double decrypt_ms = elapsed_ms([&] { decrypted = adapter.decrypt(cipher); });
    const double decode_ms = elapsed_ms([&] { (void)adapter.decode(decrypted); });

    std::printf("operation,poly_modulus_degree,payload_size,time_ms,max_abs_error,mean_abs_error,serialized_bytes\n");
    print_row("encode", poly_degree, payload_size, encode_ms, 0.0, 0.0, 0);
    print_row("encrypt", poly_degree, payload_size, encrypt_ms, 0.0, 0.0, adapter.serialized_size(cipher));
    print_row("decrypt", poly_degree, payload_size, decrypt_ms, 0.0, 0.0, 0);
    print_row("decode", poly_degree, payload_size, decode_ms, 0.0, 0.0, 0);
    print_row("add", poly_degree, payload_size, add_ms, add_accuracy.max_abs_error,
              add_accuracy.mean_abs_error, adapter.serialized_size(added));
    print_row("mul_relin_rescale", poly_degree, payload_size, mul_ms, mul_accuracy.max_abs_error,
              mul_accuracy.mean_abs_error, adapter.serialized_size(multiplied));
    print_row("rotate", poly_degree, payload_size, rotate_ms, rotate_accuracy.max_abs_error,
              rotate_accuracy.mean_abs_error, adapter.serialized_size(rotated));
    print_row("public_key", poly_degree, payload_size, 0.0, 0.0, 0.0, adapter.public_key_size());
    print_row("relin_keys", poly_degree, payload_size, 0.0, 0.0, 0.0, adapter.relin_keys_size());
    print_row("galois_keys", poly_degree, payload_size, 0.0, 0.0, 0.0, adapter.galois_keys_size());

    return 0;
}
