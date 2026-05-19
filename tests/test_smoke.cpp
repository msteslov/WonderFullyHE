#include <cstdio>
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include "m2424/abft.hpp"
#include "m2424/accuracy.hpp"
#include "m2424/seal_adapter.hpp"

static std::vector<double> head(const std::vector<double>& values, std::size_t n) {
    return std::vector<double>(values.begin(), values.begin() + std::min(values.size(), n));
}

static bool close_enough(const std::vector<double>& expected, const std::vector<double>& actual, double threshold) {
    return m2424::compare(expected, actual, threshold).ok;
}

static bool expect_runtime_error(void (*fn)()) {
    try {
        fn();
    } catch (const std::runtime_error&) {
        return true;
    }
    return false;
}

static bool expect_invalid_argument(void (*fn)()) {
    try {
        fn();
    } catch (const std::invalid_argument&) {
        return true;
    }
    return false;
}

static void invalid_profile_case() {
    (void)m2424::SealAdapter::create(m2424::CkksProfile{0, {60, 40, 60}, std::pow(2.0, 40), 0});
}

static void empty_encode_case() {
    auto adapter = m2424::SealAdapter::create(m2424::CkksProfile{8192, {60, 40, 40, 60}, std::pow(2.0, 40), 4096});
    (void)adapter.encode({});
}

static void encrypt_without_keys_case() {
    auto adapter = m2424::SealAdapter::create(m2424::CkksProfile{8192, {60, 40, 40, 60}, std::pow(2.0, 40), 4096});
    auto p = adapter.encode({1.0, 2.0});
    (void)adapter.encrypt(p);
}

static void mul_without_relin_case() {
    auto adapter = m2424::SealAdapter::create(m2424::CkksProfile{8192, {60, 40, 40, 60}, std::pow(2.0, 40), 4096});
    adapter.keygen(false, false);
    auto ct = adapter.encrypt(adapter.encode({1.0, 2.0}));
    (void)adapter.mul_relin_rescale(ct, ct);
}

static void rotate_without_galois_case() {
    auto adapter = m2424::SealAdapter::create(m2424::CkksProfile{8192, {60, 40, 40, 60}, std::pow(2.0, 40), 4096});
    adapter.keygen(false, false);
    auto ct = adapter.encrypt(adapter.encode({1.0, 2.0}));
    (void)adapter.rotate(ct, 1);
}

int main() {
    const std::size_t poly_modulus_degree = 8192;
    const std::size_t slots = poly_modulus_degree / 2;
    m2424::CkksProfile prof{poly_modulus_degree, {60, 40, 40, 60}, std::pow(2.0, 40), slots};
    auto adapter = m2424::SealAdapter::create(prof);
    adapter.keygen(true, true);

    const std::size_t N = 64;
    std::vector<double> input; input.reserve(N);
    for (std::size_t i = 0; i < N; ++i) input.push_back(std::sin(static_cast<double>(i) / 10.0));

    auto p = adapter.encode(input);
    auto ct = adapter.encrypt(p);
    auto ct2 = adapter.mul_relin_rescale(ct, ct);
    auto p2 = adapter.decrypt(ct2);
    auto out = adapter.decode(p2);

    std::vector<double> ref; ref.reserve(N);
    for (double x : input) ref.push_back(x * x);
    std::vector<double> out_head(out.begin(), out.begin() + std::min(out.size(), N));

    auto mul_accuracy = m2424::compare(ref, out_head, 1e-3);

    auto ct_add = adapter.add(ct, ct);
    std::vector<double> add_ref; add_ref.reserve(N);
    for (double x : input) add_ref.push_back(x + x);
    auto add_out = head(adapter.decode(adapter.decrypt(ct_add)), N);

    auto ct_sub = adapter.sub(ct_add, ct);
    auto sub_out = head(adapter.decode(adapter.decrypt(ct_sub)), N);

    auto ct_rot = adapter.rotate(ct, 1);
    auto rot_out = head(adapter.decode(adapter.decrypt(ct_rot)), N);
    std::vector<double> rot_ref; rot_ref.reserve(N);
    for (std::size_t i = 0; i < N; ++i) rot_ref.push_back(input[(i + 1) % N]);

    auto abft_payload = m2424::abft::append_checksum(input);
    auto abft_ct = adapter.add(adapter.encrypt(adapter.encode(abft_payload)), adapter.encrypt(adapter.encode(abft_payload)));
    auto abft_out = adapter.decode(adapter.decrypt(abft_ct));
    auto abft_check = m2424::abft::verify_appended_checksum(abft_out, N, 1e-5);

    bool ok = std::isfinite(mul_accuracy.max_abs_error) && std::isfinite(mul_accuracy.mean_abs_error) && mul_accuracy.ok
        && close_enough(add_ref, add_out, 1e-5)
        && close_enough(input, sub_out, 1e-5)
        && close_enough(rot_ref, rot_out, 1e-5)
        && abft_check.ok
        && adapter.slot_count() == slots
        && adapter.serialized_size(ct) > 0
        && adapter.scale(ct) > 0.0
        && adapter.coeff_modulus_size(ct) > 0
        && adapter.public_key_size() > 0
        && adapter.relin_keys_size() > 0
        && adapter.galois_keys_size() > 0
        && expect_invalid_argument(invalid_profile_case)
        && expect_invalid_argument(empty_encode_case)
        && expect_runtime_error(encrypt_without_keys_case)
        && expect_runtime_error(mul_without_relin_case)
        && expect_runtime_error(rotate_without_galois_case);
    std::printf("[test_smoke] max=%.6e mean=%.6e threshold=1e-3 => %s\n",
               mul_accuracy.max_abs_error, mul_accuracy.mean_abs_error, ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
