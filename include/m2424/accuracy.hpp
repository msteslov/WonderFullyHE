#pragma once

#include <cstddef>
#include <vector>

namespace m2424 {

/// Целевой абсолютный предел ошибки для всех проверяемых CKKS-результатов.
inline constexpr double kTargetAbsoluteError = 1e-9;

struct AccuracyReport {
    double max_abs_error{};
    double mean_abs_error{};
    double rms_error{};
    double max_relative_error{};
    double tolerance{};
    std::size_t worst_index{};
    bool size_ok{};
    bool finite_ok{};
    bool ok{};
};

double max_abs_error(const std::vector<double>& expected, const std::vector<double>& actual);
double mean_abs_error(const std::vector<double>& expected, const std::vector<double>& actual);
AccuracyReport compare(const std::vector<double>& expected, const std::vector<double>& actual, double tolerance);

} // namespace m2424
