#pragma once

#include <cstddef>
#include <vector>

namespace m2424::abft {

struct ChecksumResult {
    double expected{};
    double observed{};
    double abs_error{};
    bool ok{};
};

double checksum(const std::vector<double>& values);
std::vector<double> append_checksum(const std::vector<double>& values);
ChecksumResult verify_appended_checksum(const std::vector<double>& values, std::size_t payload_size, double tolerance);
ChecksumResult verify_checksum_value(const std::vector<double>& values, std::size_t payload_size,
                                     double expected, double tolerance);

} // namespace m2424::abft
