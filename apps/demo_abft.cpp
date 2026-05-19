#include "m2424/abft.hpp"
#include "m2424/seal_adapter.hpp"

#include <cmath>
#include <cstdio>
#include <vector>

namespace {

void print_check(const char* name, const m2424::abft::ChecksumResult& check) {
    std::printf("%s expected=%.12e observed=%.12e abs_error=%.6e => %s\n",
                name, check.expected, check.observed, check.abs_error, check.ok ? "PASS" : "FAIL");
}

std::vector<double> head(const std::vector<double>& values, std::size_t size) {
    return std::vector<double>(values.begin(), values.begin() + size);
}

} // namespace

int main() {
    const std::size_t poly_degree = 8192;
    const std::size_t payload_size = 16;
    m2424::CkksProfile profile{poly_degree, {60, 40, 40, 60}, std::pow(2.0, 40), poly_degree / 2};

    auto adapter = m2424::SealAdapter::create(profile);
    adapter.keygen(true, true);

    std::vector<double> lhs;
    std::vector<double> rhs;
    lhs.reserve(payload_size);
    rhs.reserve(payload_size);
    for (std::size_t i = 0; i < payload_size; ++i) {
        lhs.push_back(std::sin(static_cast<double>(i) / 8.0));
        rhs.push_back(std::cos(static_cast<double>(i) / 9.0));
    }

    auto lhs_ct = adapter.encrypt(adapter.encode(m2424::abft::append_checksum(lhs)));
    auto rhs_ct = adapter.encrypt(adapter.encode(m2424::abft::append_checksum(rhs)));

    auto sum_ct = adapter.add(lhs_ct, rhs_ct);
    auto add_check = m2424::abft::verify_appended_checksum(adapter.decode(adapter.decrypt(sum_ct)), payload_size, 1e-5);

    auto diff_ct = adapter.sub(lhs_ct, rhs_ct);
    auto sub_check = m2424::abft::verify_appended_checksum(adapter.decode(adapter.decrypt(diff_ct)), payload_size, 1e-5);

    std::vector<double> product_ref;
    product_ref.reserve(payload_size);
    for (std::size_t i = 0; i < payload_size; ++i) {
        product_ref.push_back(lhs[i] * rhs[i]);
    }
    auto mul_ct = adapter.mul_relin_rescale(adapter.encrypt(adapter.encode(lhs)), adapter.encrypt(adapter.encode(rhs)));
    auto mul_check = m2424::abft::verify_checksum_value(head(adapter.decode(adapter.decrypt(mul_ct)), payload_size),
                                                       payload_size, m2424::abft::checksum(product_ref), 1e-5);

    std::vector<double> rotate_payload;
    rotate_payload.reserve(adapter.slot_count());
    for (std::size_t i = 0; i < adapter.slot_count(); ++i) {
        rotate_payload.push_back(0.001 * std::sin(static_cast<double>(i) / 17.0));
    }
    auto rotate_ct = adapter.rotate(adapter.encrypt(adapter.encode(rotate_payload)), 7);
    auto rotate_check = m2424::abft::verify_checksum_value(adapter.decode(adapter.decrypt(rotate_ct)),
                                                          adapter.slot_count(), m2424::abft::checksum(rotate_payload), 1e-4);

    print_check("abft_add_appended_checksum", add_check);
    print_check("abft_sub_appended_checksum", sub_check);
    print_check("abft_mul_reference_checksum", mul_check);
    print_check("abft_rotate_sum_invariant", rotate_check);

    return add_check.ok && sub_check.ok && mul_check.ok && rotate_check.ok ? 0 : 1;
}
