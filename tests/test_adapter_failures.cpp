#include "m2424/profiles.hpp"
#include "m2424/seal_adapter.hpp"

#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

bool throws_invalid_argument(void (*fn)()) {
    try {
        fn();
    } catch (const std::invalid_argument&) {
        return true;
    }
    return false;
}

bool throws_runtime_error(void (*fn)()) {
    try {
        fn();
    } catch (const std::runtime_error&) {
        return true;
    }
    return false;
}

bool throws_exception(void (*fn)()) {
    try {
        fn();
    } catch (const std::exception&) {
        return true;
    }
    return false;
}

void encode_nan_case() {
    auto adapter = m2424::SealAdapter::create(m2424::profiles::basic_ckks());
    (void)adapter.encode({1.0, std::numeric_limits<double>::quiet_NaN()});
}

void encode_too_many_slots_case() {
    auto adapter = m2424::SealAdapter::create(m2424::profiles::fast_demo_ckks());
    std::vector<double> values(adapter.slot_count() + 1, 1.0);
    (void)adapter.encode(values);
}

void encrypt_without_key_case() {
    auto adapter = m2424::SealAdapter::create(m2424::profiles::basic_ckks());
    (void)adapter.encrypt(adapter.encode({1.0}));
}

void rotate_without_key_case() {
    auto adapter = m2424::SealAdapter::create(m2424::profiles::basic_ckks());
    adapter.keygen(true, false);
    auto encrypted = adapter.encrypt(adapter.encode({1.0, 2.0}));
    (void)adapter.rotate(encrypted, 1);
}

void load_empty_cipher_case() {
    auto adapter = m2424::SealAdapter::create(m2424::profiles::basic_ckks());
    (void)adapter.load_cipher({});
}

void load_corrupt_public_key_case() {
    auto adapter = m2424::SealAdapter::create(m2424::profiles::basic_ckks());
    (void)adapter.load_public_key({1, 2, 3, 4, 5});
}

void unsafe_scale_invalid_multiplier_case() {
    auto adapter = m2424::SealAdapter::create(m2424::profiles::basic_ckks());
    adapter.keygen(true, true);
    auto encrypted = adapter.encrypt(adapter.encode({1.0}));
    (void)adapter.unsafe_reinterpret_scale_for_diagnostics(encrypted, 0.0);
}

} // namespace

int main() {
    const bool ok = throws_invalid_argument(encode_nan_case)
        && throws_invalid_argument(encode_too_many_slots_case)
        && throws_runtime_error(encrypt_without_key_case)
        && throws_runtime_error(rotate_without_key_case)
        && throws_invalid_argument(load_empty_cipher_case)
        && throws_exception(load_corrupt_public_key_case)
        && throws_invalid_argument(unsafe_scale_invalid_multiplier_case);

    std::printf("[test_adapter_failures] negative_cases=%s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
