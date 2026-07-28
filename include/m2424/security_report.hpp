#pragma once

#include "m2424/seal_adapter.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace m2424 {

enum class SecurityLevel {
    None,
    TC128,
    TC192,
    TC256
};

struct SecurityReport {
    std::string profile_name;
    std::size_t polyModulusDegree{};
    int total_coeff_modulus_bits{};
    int tc128_limit{};
    int tc192_limit{};
    int tc256_limit{};
    bool passes_tc128{};
    bool passes_tc192{};
    bool passes_tc256{};
    SecurityLevel effective_level{SecurityLevel::None};
};

int coeff_modulus_max_bit_count(std::size_t polyModulusDegree, SecurityLevel level);
int total_coeff_modulus_bits(const CkksProfile& profile);
SecurityReport analyze_security(const std::string& profile_name, const CkksProfile& profile);
SecurityLevel project_minimum_security(const std::vector<SecurityReport>& reports);
const char* to_string(SecurityLevel level) noexcept;
std::string security_report_csv_header();
std::string to_csv_row(const SecurityReport& report);

} // namespace m2424
