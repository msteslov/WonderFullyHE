#include "m2424/profiles.hpp"

#include <cmath>
#include <stdexcept>

namespace m2424::profiles {

CkksProfile fast_demo_ckks() {
    return CkksProfile{4096, {40, 30, 30}, std::pow(2.0, 30), 2048};
}

CkksProfile basic_ckks() {
    return CkksProfile{8192, {60, 40, 40, 60}, std::pow(2.0, 40), 4096};
}

CkksProfile balanced_ckks() {
    return CkksProfile{8192, {50, 40, 40, 40, 48}, std::pow(2.0, 40), 4096};
}

CkksProfile depth_ckks() {
    return CkksProfile{16384, {60, 40, 40, 40, 40, 60}, std::pow(2.0, 40), 8192};
}

CkksProfile high_precision_ckks() {
    return CkksProfile{16384, {60, 50, 50, 50, 60}, std::pow(2.0, 50), 8192};
}

CkksProfile boot_ckks() {
    return CkksProfile{16384, {60, 40, 40, 40, 40, 40, 40, 40, 60}, std::pow(2.0, 40), 8192};
}

CkksProfile boot_deep_ckks() {
    return CkksProfile{32768,
                       {60,
                        40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
                        40, 40, 40, 40, 40, 40, 40, 40, 40,
                        60},
                       std::pow(2.0, 40),
                       16384};
}

CkksProfile precision_boot_deep_ckks() {
    return CkksProfile{32768,
                       {60,
                        50, 50, 50, 50, 50, 50, 50, 50,
                        50, 50, 50, 50, 50, 50, 50,
                        60},
                       std::pow(2.0, 50),
                       16384};
}

CkksProfile precision_boot_ultra_ckks_59() {
    return CkksProfile{32768,
                       {60, 60, 60, 60, 60, 60, 60,
                        60, 60, 60, 60, 60, 60, 60},
                       std::pow(2.0, 59),
                       16384};
}

std::vector<std::pair<std::string, CkksProfile>> all() {
    return {
        {"fast_demo_ckks", fast_demo_ckks()},
        {"basic_ckks", basic_ckks()},
        {"balanced_ckks", balanced_ckks()},
        {"depth_ckks", depth_ckks()},
        {"high_precision_ckks", high_precision_ckks()},
        {"boot_ckks", boot_ckks()},
        {"boot_deep_ckks", boot_deep_ckks()},
        {"precision_boot_deep_ckks", precision_boot_deep_ckks()},
        {"precision_boot_ultra_ckks_59", precision_boot_ultra_ckks_59()}
    };
}

CkksProfile by_name(const std::string& name) {
    for (const auto& entry : all()) {
        if (entry.first == name) {
            return entry.second;
        }
    }
    throw std::invalid_argument("unknown CKKS profile: " + name);
}

} // namespace m2424::profiles
