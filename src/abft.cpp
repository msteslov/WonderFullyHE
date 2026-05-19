#include "m2424/abft.hpp"

#include <cmath>
#include <stdexcept>

namespace m2424::abft {

double checksum(const std::vector<double>& values) {
    double total = 0.0;
    for (double value : values) {
        total += value;
    }
    return total;
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

    double expected = 0.0;
    for (std::size_t i = 0; i < payload_size; ++i) {
        expected += values[i];
    }
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

    double observed = 0.0;
    for (std::size_t i = 0; i < payload_size; ++i) {
        observed += values[i];
    }
    const double abs_error = std::fabs(expected - observed);
    return ChecksumResult{expected, observed, abs_error, abs_error <= tolerance};
}

} // namespace m2424::abft
