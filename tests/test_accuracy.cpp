#include "m2424/accuracy.hpp"

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

void size_mismatch_case() {
    (void)m2424::compare({1.0, 2.0}, {1.0}, 1e-6);
}

void negative_tolerance_case() {
    (void)m2424::compare({1.0}, {1.0}, -1.0);
}

void max_abs_nan_case() {
    (void)m2424::max_abs_error({1.0}, {std::numeric_limits<double>::quiet_NaN()});
}

} // namespace

int main() {
    const auto report = m2424::compare({1.0, 2.0, 4.0}, {1.0, 2.1, 3.9}, 0.11);
    const auto nan_report = m2424::compare(
        {1.0, 2.0},
        {1.0, std::numeric_limits<double>::quiet_NaN()},
        1e-6);
    const auto inf_report = m2424::compare(
        {1.0, std::numeric_limits<double>::infinity()},
        {1.0, 2.0},
        1e-6);

    const bool ok = report.ok
        && report.size_ok
        && report.finite_ok
        && report.worst_index == 1
        && report.max_abs_error > 0.09
        && report.mean_abs_error > 0.06
        && report.rms_error > 0.08
        && report.max_relative_error > 0.02
        && !nan_report.ok
        && nan_report.size_ok
        && !nan_report.finite_ok
        && !inf_report.ok
        && !inf_report.finite_ok
        && throws_invalid_argument(size_mismatch_case)
        && throws_invalid_argument(negative_tolerance_case)
        && throws_invalid_argument(max_abs_nan_case);

    std::printf("[test_accuracy] finite=%s nan_fail=%s inf_fail=%s => %s\n",
                report.finite_ok ? "PASS" : "FAIL",
                !nan_report.finite_ok ? "PASS" : "FAIL",
                !inf_report.finite_ok ? "PASS" : "FAIL",
                ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
