#include "m2424/abft.hpp"
#include "m2424/seal_adapter.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

namespace {

double max_abs_error(const std::vector<double>& a, const std::vector<double>& b) {
    const std::size_t n = std::min(a.size(), b.size());
    double result = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        result = std::max(result, std::fabs(a[i] - b[i]));
    }
    return result;
}

std::vector<double> encrypted_prefix_sum(m2424::SealAdapter& adapter, const m2424::Cipher& input,
                                         std::size_t payload_size) {
    auto acc = input;
    for (std::size_t step = 1; step < payload_size; step <<= 1) {
        acc = adapter.add(acc, adapter.rotate(acc, static_cast<int>(step)));
    }
    return adapter.decode(adapter.decrypt(acc));
}

} // namespace

int main() {
    const std::size_t poly_degree = 8192;
    const std::size_t payload_size = 16;
    m2424::CkksProfile profile{poly_degree, {60, 40, 40, 60}, std::pow(2.0, 40), poly_degree / 2};

    auto adapter = m2424::SealAdapter::create(profile);
    adapter.keygen(false, true);

    std::vector<double> data(adapter.slot_count(), 0.0);
    for (std::size_t i = 0; i < payload_size; ++i) {
        data[i] = 10.0 + 0.5 * static_cast<double>(i) + std::sin(static_cast<double>(i) / 3.0);
    }

    const double expected_sum = m2424::abft::checksum(std::vector<double>(data.begin(), data.begin() + payload_size));
    const double expected_mean = expected_sum / static_cast<double>(payload_size);

    auto encrypted = adapter.encrypt(adapter.encode(data));
    auto summed_slots = encrypted_prefix_sum(adapter, encrypted, payload_size);
    const double encrypted_sum = summed_slots[0];
    const double encrypted_mean = encrypted_sum / static_cast<double>(payload_size);

    std::vector<double> expected{expected_sum, expected_mean};
    std::vector<double> actual{encrypted_sum, encrypted_mean};
    const double error = max_abs_error(expected, actual);

    std::printf("secure_statistics\n");
    std::printf("payload_size=%zu\n", payload_size);
    std::printf("sum_cpu=%.12e\n", expected_sum);
    std::printf("sum_encrypted=%.12e\n", encrypted_sum);
    std::printf("mean_cpu=%.12e\n", expected_mean);
    std::printf("mean_encrypted=%.12e\n", encrypted_mean);
    std::printf("max_abs_error=%.6e\n", error);
    std::printf("status=%s\n", error < 1e-5 ? "PASS" : "FAIL");

    return error < 1e-5 ? 0 : 1;
}
