#include "m2424/accuracy.hpp"
#include "m2424/profiles.hpp"
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

void print_row(const char* profile_name, const char* operation, std::size_t poly_degree, std::size_t payload_size,
               double time_ms, double max_error, double mean_error, std::size_t bytes) {
    std::printf("%s,%s,%zu,%zu,%.6f,%.6e,%.6e,%zu\n",
                profile_name,
                operation, poly_degree, payload_size, time_ms, max_error, mean_error, bytes);
}

std::vector<double> payload_head(const std::vector<double>& values, std::size_t payload_size) {
    return std::vector<double>(values.begin(), values.begin() + std::min(values.size(), payload_size));
}

} // namespace

int main(int argc, char** argv) {
    const std::size_t payload_size = 64;
    const std::string profile_name = argc > 1 ? argv[1] : "basic_ckks";
    const auto profile = m2424::profiles::by_name(profile_name);
    const std::size_t poly_degree = profile.polyModulusDegree;

    auto adapter = m2424::SealAdapter::create(profile);
    adapter.generateKeys(true, true);

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
    double mul_ms = elapsed_ms([&] { multiplied = adapter.rescaleToNext(adapter.relinearize(adapter.multiply(cipher, cipher))); });
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

    std::printf("profile,operation,polyModulusDegree,payload_size,time_ms,max_abs_error,mean_abs_error,serialized_bytes\n");
    print_row(profile_name.c_str(), "encode", poly_degree, payload_size, encode_ms, 0.0, 0.0, 0);
    print_row(profile_name.c_str(), "encrypt", poly_degree, payload_size, encrypt_ms, 0.0, 0.0, adapter.serializedSize(cipher));
    print_row(profile_name.c_str(), "decrypt", poly_degree, payload_size, decrypt_ms, 0.0, 0.0, 0);
    print_row(profile_name.c_str(), "decode", poly_degree, payload_size, decode_ms, 0.0, 0.0, 0);
    print_row(profile_name.c_str(), "add", poly_degree, payload_size, add_ms, add_accuracy.max_abs_error,
              add_accuracy.mean_abs_error, adapter.serializedSize(added));
    print_row(profile_name.c_str(), "multiplyRelinearizeAndRescale", poly_degree, payload_size, mul_ms, mul_accuracy.max_abs_error,
              mul_accuracy.mean_abs_error, adapter.serializedSize(multiplied));
    print_row(profile_name.c_str(), "rotate", poly_degree, payload_size, rotate_ms, rotate_accuracy.max_abs_error,
              rotate_accuracy.mean_abs_error, adapter.serializedSize(rotated));
    print_row(profile_name.c_str(), "public_key", poly_degree, payload_size, 0.0, 0.0, 0.0, adapter.publicKeySize());
    print_row(profile_name.c_str(), "relin_keys", poly_degree, payload_size, 0.0, 0.0, 0.0, adapter.relinKeysSize());
    print_row(profile_name.c_str(), "galois_keys", poly_degree, payload_size, 0.0, 0.0, 0.0, adapter.galoisKeysSize());

    return 0;
}
