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
        result = std::max(result, std::fabs(expected[i] - actual[i]));
    }
    return result;
}

double mean_abs_error(const std::vector<double>& expected, const std::vector<double>& actual) {
    const std::size_t n = checked_size(expected, actual);
    double total = 0.0;
    double correction = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double value = std::fabs(expected[i] - actual[i]);
        const double y = value - correction;
        const double t = total + y;
        correction = (t - total) - y;
        total = t;
    }
    return total / static_cast<double>(n);
}

AccuracyReport compare(const std::vector<double>& expected, const std::vector<double>& actual, double tolerance) {
    if (tolerance < 0.0) {
        throw std::invalid_argument("tolerance must be non-negative");
    }
    const double max_error = max_abs_error(expected, actual);
    const double mean_error = mean_abs_error(expected, actual);
    return AccuracyReport{max_error, mean_error, tolerance, max_error <= tolerance};
}

} // namespace m2424
