#include "m2424/security_report.hpp"

#include <algorithm>
#include <sstream>
#include <stdexcept>

namespace m2424 {

const char* to_string(SecurityLevel level) noexcept {
    switch (level) {
        case SecurityLevel::None:
            return "none";
        case SecurityLevel::TC128:
            return "tc128";
        case SecurityLevel::TC192:
            return "tc192";
        case SecurityLevel::TC256:
            return "tc256";
    }
    return "unknown";
}

int coeff_modulus_max_bit_count(std::size_t poly_modulus_degree, SecurityLevel level) {
    switch (poly_modulus_degree) {
        case 1024:
            if (level == SecurityLevel::TC128) return 27;
            break;
        case 2048:
            if (level == SecurityLevel::TC128) return 54;
            break;
        case 4096:
            if (level == SecurityLevel::TC128) return 109;
            if (level == SecurityLevel::TC192) return 75;
            if (level == SecurityLevel::TC256) return 58;
            break;
        case 8192:
            if (level == SecurityLevel::TC128) return 218;
            if (level == SecurityLevel::TC192) return 152;
            if (level == SecurityLevel::TC256) return 118;
            break;
        case 16384:
            if (level == SecurityLevel::TC128) return 438;
            if (level == SecurityLevel::TC192) return 305;
            if (level == SecurityLevel::TC256) return 237;
            break;
        case 32768:
            if (level == SecurityLevel::TC128) return 881;
            if (level == SecurityLevel::TC192) return 611;
            if (level == SecurityLevel::TC256) return 476;
            break;
    }
    return 0;
}

int total_coeff_modulus_bits(const CkksProfile& profile) {
    int total = 0;
    for (int bits : profile.coeff_modulus_bits) {
        if (bits <= 0) {
            throw std::invalid_argument("coeff_modulus_bits entries must be positive");
        }
        total += bits;
    }
    return total;
}

static SecurityLevel effective_level(bool tc128, bool tc192, bool tc256) {
    if (tc256) return SecurityLevel::TC256;
    if (tc192) return SecurityLevel::TC192;
    if (tc128) return SecurityLevel::TC128;
    return SecurityLevel::None;
}

SecurityReport analyze_security(const std::string& profile_name, const CkksProfile& profile) {
    if (profile.poly_modulus_degree == 0) {
        throw std::invalid_argument("poly_modulus_degree must be positive");
    }
    const int total = total_coeff_modulus_bits(profile);
    const int tc128 = coeff_modulus_max_bit_count(profile.poly_modulus_degree, SecurityLevel::TC128);
    const int tc192 = coeff_modulus_max_bit_count(profile.poly_modulus_degree, SecurityLevel::TC192);
    const int tc256 = coeff_modulus_max_bit_count(profile.poly_modulus_degree, SecurityLevel::TC256);

    const bool pass128 = tc128 > 0 && total <= tc128;
    const bool pass192 = tc192 > 0 && total <= tc192;
    const bool pass256 = tc256 > 0 && total <= tc256;

    return SecurityReport{
        profile_name,
        profile.poly_modulus_degree,
        total,
        tc128,
        tc192,
        tc256,
        pass128,
        pass192,
        pass256,
        effective_level(pass128, pass192, pass256)
    };
}

SecurityLevel project_minimum_security(const std::vector<SecurityReport>& reports) {
    if (reports.empty()) {
        throw std::invalid_argument("security reports must not be empty");
    }
    auto rank = [](SecurityLevel level) {
        switch (level) {
            case SecurityLevel::None: return 0;
            case SecurityLevel::TC128: return 1;
            case SecurityLevel::TC192: return 2;
            case SecurityLevel::TC256: return 3;
        }
        return 0;
    };

    SecurityLevel result = reports.front().effective_level;
    for (const auto& report : reports) {
        if (rank(report.effective_level) < rank(result)) {
            result = report.effective_level;
        }
    }
    return result;
}

std::string security_report_csv_header() {
    return "profile,poly_modulus_degree,total_coeff_modulus_bits,tc128_limit,tc192_limit,tc256_limit,passes_tc128,passes_tc192,passes_tc256,effective_security_level";
}

std::string to_csv_row(const SecurityReport& report) {
    std::ostringstream out;
    out << report.profile_name << ","
        << report.poly_modulus_degree << ","
        << report.total_coeff_modulus_bits << ","
        << report.tc128_limit << ","
        << report.tc192_limit << ","
        << report.tc256_limit << ","
        << (report.passes_tc128 ? "PASS" : "FAIL") << ","
        << (report.passes_tc192 ? "PASS" : "FAIL") << ","
        << (report.passes_tc256 ? "PASS" : "FAIL") << ","
        << to_string(report.effective_level);
    return out.str();
}

} // namespace m2424
