#include "m2424/accuracy.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace m2424 {
namespace {

static std::size_t checked_size(const std::vector<double>& expected, const std::vector<double>& actual) {
    if (expected.empty() || actual.empty()) {
        throw std::invalid_argument("accuracy comparison requires non-empty vectors");
    }
    if (expected.size() != actual.size()) {
        throw std::invalid_argument("accuracy comparison requires vectors of equal size");
    }
    return expected.size();
}

} // namespace

double max_abs_error(const std::vector<double>& expected, const std::vector<double>& actual) {
    const std::size_t n = checked_size(expected, actual);
    double result = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        if (!std::isfinite(expected[i]) || !std::isfinite(actual[i])) {
            throw std::invalid_argument("accuracy comparison requires finite values");
        }
        result = std::max(result, std::fabs(expected[i] - actual[i]));
    }
    return result;
}

double mean_abs_error(const std::vector<double>& expected, const std::vector<double>& actual) {
    const std::size_t n = checked_size(expected, actual);
    double total = 0.0;
    double correction = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        if (!std::isfinite(expected[i]) || !std::isfinite(actual[i])) {
            throw std::invalid_argument("accuracy comparison requires finite values");
        }
        const double value = std::fabs(expected[i] - actual[i]);
        const double y = value - correction;
        const double t = total + y;
        correction = (t - total) - y;
        total = t;
    }
    return total / static_cast<double>(n);
}

AccuracyReport compare(const std::vector<double>& expected, const std::vector<double>& actual, double tolerance) {
    if (!std::isfinite(tolerance) || tolerance < 0.0) {
        throw std::invalid_argument("tolerance must be non-negative");
    }
    const std::size_t n = checked_size(expected, actual);

    AccuracyReport report;
    report.tolerance = tolerance;
    report.size_ok = true;
    report.finite_ok = true;

    double total_abs = 0.0;
    double total_square = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        if (!std::isfinite(expected[i]) || !std::isfinite(actual[i])) {
            report.finite_ok = false;
            report.ok = false;
            report.worst_index = i;
            return report;
        }
        const double abs_error = std::fabs(expected[i] - actual[i]);
        const double relative_denominator = std::max(std::fabs(expected[i]), tolerance);
        const double relative_error = abs_error / relative_denominator;
        total_abs += abs_error;
        total_square += abs_error * abs_error;
        if (abs_error > report.max_abs_error) {
            report.max_abs_error = abs_error;
            report.worst_index = i;
        }
        report.max_relative_error = std::max(report.max_relative_error, relative_error);
    }
    report.mean_abs_error = total_abs / static_cast<double>(n);
    report.rms_error = std::sqrt(total_square / static_cast<double>(n));
    report.ok = report.finite_ok && report.max_abs_error <= tolerance;
    return report;
}

} // namespace m2424
