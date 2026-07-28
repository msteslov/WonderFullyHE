#include "m2424/profile_report.hpp"

#include <cmath>
#include <sstream>
#include <stdexcept>

namespace m2424 {

static std::string join_bits(const std::vector<int>& bits) {
    std::ostringstream out;
    for (std::size_t i = 0; i < bits.size(); ++i) {
        if (i != 0) out << "-";
        out << bits[i];
    }
    return out.str();
}

static int total_bits(const std::vector<int>& bits) {
    int total = 0;
    for (int value : bits) total += value;
    return total;
}

static std::size_t estimated_mul_depth(const CkksProfile& profile) {
    if (profile.coeffModulusBits.size() <= 2) {
        return 0;
    }
    const double log_scale = std::log2(profile.scale);
    std::size_t usable = 0;
    for (std::size_t i = 1; i + 1 < profile.coeffModulusBits.size(); ++i) {
        if (static_cast<double>(profile.coeffModulusBits[i]) + 1.0 >= log_scale) {
            ++usable;
        }
    }
    return usable;
}

ProfileReport describe_profile(const std::string& name, const CkksProfile& profile) {
    if (profile.polyModulusDegree == 0) {
        throw std::invalid_argument("polyModulusDegree must be positive");
    }
    if (profile.scale <= 0.0 || !std::isfinite(profile.scale)) {
        throw std::invalid_argument("scale must be positive and finite");
    }

    return ProfileReport{
        name,
        profile.polyModulusDegree,
        profile.polyModulusDegree / 2,
        join_bits(profile.coeffModulusBits),
        total_bits(profile.coeffModulusBits),
        profile.scale,
        std::log2(profile.scale),
        estimated_mul_depth(profile)
    };
}

std::string profile_report_csv_header() {
    return "profile,polyModulusDegree,slotCount,coeffModulusBits,total_coeff_modulus_bits,scale,scale_log2,estimated_mul_depth";
}

std::string to_csv_row(const ProfileReport& report) {
    std::ostringstream out;
    out << report.name << ","
        << report.polyModulusDegree << ","
        << report.slotCount << ","
        << report.coeffModulusBits << ","
        << report.total_coeff_modulus_bits << ","
        << report.scale << ","
        << report.scale_log2 << ","
        << report.estimated_mul_depth;
    return out.str();
}

} // namespace m2424
