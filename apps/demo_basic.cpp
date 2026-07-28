#include "m2424/accuracy.hpp"
#include "m2424/m2424.hpp"
#include "m2424/seal_adapter.hpp"

#include <cstdio>
#include <vector>
#include <cmath>
#include <algorithm>

static std::vector<double> rotate_left(const std::vector<double>& v, int steps) {
    if (v.empty()) return {};
    const std::size_t n = v.size();
    const std::size_t s = static_cast<std::size_t>((steps % static_cast<int>(n) + static_cast<int>(n)) % static_cast<int>(n));
    std::vector<double> out(n);
    for (std::size_t i = 0; i < n; ++i) {
        out[i] = v[(i + s) % n];
    }
    return out;
}

int main() {
    const auto prof = m2424::profiles::basic_ckks();
    auto adapter = m2424::SealAdapter::create(prof);
    adapter.generateKeys(true, true);

    // Prepare input vector within slots, modest length for a quick demo
    const std::size_t N = 64;
    std::vector<double> input;
    input.reserve(N);
    for (std::size_t i = 0; i < N; ++i) {
        input.push_back(std::sin(static_cast<double>(i) / 10.0));
    }

    // Encrypt
    auto p = adapter.encode(input);
    auto ct = adapter.encrypt(p);

    // ct2 = ct * ct (with relin + rescale), ct3 = rotate(ct2, 1), ct4 = ct2 + ct3
    auto ct2 = adapter.rescaleToNext(adapter.relinearize(adapter.multiply(ct, ct)));
    auto ct3 = adapter.rotate(ct2, 1);
    auto ct4 = adapter.add(ct2, ct3);

    // Decrypt and decode result
    auto p_res = adapter.decrypt(ct4);
    auto out = adapter.decode(p_res);

    // Reference computation
    std::vector<double> ref2; ref2.reserve(N);
    for (double x : input) ref2.push_back(x * x);
    auto ref3 = rotate_left(ref2, 1);
    std::vector<double> ref4(N);
    for (std::size_t i = 0; i < N; ++i) ref4[i] = ref2[i] + ref3[i];

    // Compare only the first N entries from decoded output
    if (out.size() < N) {
        std::fprintf(stderr, "decoded vector smaller than expected: %zu < %zu\n", out.size(), N);
        return 1;
    }
    std::vector<double> out_head(out.begin(), out.begin() + N);

    auto accuracy = m2424::compare(ref4, out_head, 1e-5);
    std::printf("max_abs_error=%.6e\nmean_abs_error=%.6e\n", accuracy.max_abs_error, accuracy.mean_abs_error);
    (void)m2424::version();
    return 0;
}
