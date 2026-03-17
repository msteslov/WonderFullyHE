#include <cstdio>
#include <vector>
#include <cmath>
#include <algorithm>
#include "m2424/seal_adapter.hpp"

static double max_abs_error(const std::vector<double>& a, const std::vector<double>& b) {
    const std::size_t n = std::min(a.size(), b.size());
    double m = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        double e = std::fabs(a[i] - b[i]);
        if (e > m) m = e;
    }
    return m;
}

static double mean_abs_error(const std::vector<double>& a, const std::vector<double>& b) {
    const std::size_t n = std::min(a.size(), b.size());
    if (n == 0) return 0.0;
    long double acc = 0.0L;
    for (std::size_t i = 0; i < n; ++i) acc += std::fabsl(static_cast<long double>(a[i] - b[i]));
    return static_cast<double>(acc / n);
}

int main() {
    const std::size_t poly_modulus_degree = 8192;
    const std::size_t slots = poly_modulus_degree / 2;
    m2424::CkksProfile prof{poly_modulus_degree, {60, 40, 40, 60}, std::pow(2.0, 40), slots};
    auto adapter = m2424::SealAdapter::create(prof);
    adapter.keygen(true, false); // relin needed, galois not needed

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

    double max_err = max_abs_error(ref, out_head);
    double mean_err = mean_abs_error(ref, out_head);

    bool ok = std::isfinite(max_err) && std::isfinite(mean_err) && max_err < 1e-3;
    std::printf("[test_smoke] max=%.6e mean=%.6e threshold=1e-3 => %s\n",
               max_err, mean_err, ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
