#include "m2424/abft.hpp"

#include <cmath>
#include <stdexcept>

namespace m2424::abft {
namespace {

double compensated_sum(const std::vector<double>& values, std::size_t count) {
    double total = 0.0;
    double correction = 0.0;
    for (std::size_t i = 0; i < count; ++i) {
        const double y = values[i] - correction;
        const double t = total + y;
        correction = (t - total) - y;
        total = t;
    }
    return total;
}

} // namespace

double checksum(const std::vector<double>& values) {
    return compensated_sum(values, values.size());
}

std::vector<double> append_checksum(const std::vector<double>& values) {
    if (values.empty()) {
        throw std::invalid_argument("payload must not be empty");
    }
    std::vector<double> encoded = values;
    encoded.push_back(checksum(values));
    return encoded;
}

ChecksumResult verify_appended_checksum(const std::vector<double>& values, std::size_t payload_size, double tolerance) {
    if (payload_size == 0) {
        throw std::invalid_argument("payload_size must be positive");
    }
    if (values.size() <= payload_size) {
        throw std::invalid_argument("values do not contain an appended checksum slot");
    }
    if (tolerance < 0.0) {
        throw std::invalid_argument("tolerance must be non-negative");
    }

    const double expected = compensated_sum(values, payload_size);
    const double observed = values[payload_size];
    const double abs_error = std::fabs(expected - observed);
    return ChecksumResult{expected, observed, abs_error, abs_error <= tolerance};
}

ChecksumResult verify_checksum_value(const std::vector<double>& values, std::size_t payload_size,
                                     double expected, double tolerance) {
    if (payload_size == 0) {
        throw std::invalid_argument("payload_size must be positive");
    }
    if (values.size() < payload_size) {
        throw std::invalid_argument("values are smaller than payload_size");
    }
    if (tolerance < 0.0) {
        throw std::invalid_argument("tolerance must be non-negative");
    }

    const double observed = compensated_sum(values, payload_size);
    const double abs_error = std::fabs(expected - observed);
    return ChecksumResult{expected, observed, abs_error, abs_error <= tolerance};
}

} // namespace m2424::abft
