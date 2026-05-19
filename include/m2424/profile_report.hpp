#pragma once

#include "m2424/seal_adapter.hpp"

#include <cstddef>
#include <string>

namespace m2424 {

struct ProfileReport {
    std::string name;
    std::size_t poly_modulus_degree{};
    std::size_t slot_count{};
    std::string coeff_modulus_bits;
    int total_coeff_modulus_bits{};
    double scale{};
    double scale_log2{};
    std::size_t estimated_mul_depth{};
};

ProfileReport describe_profile(const std::string& name, const CkksProfile& profile);
std::string profile_report_csv_header();
std::string to_csv_row(const ProfileReport& report);

} // namespace m2424
