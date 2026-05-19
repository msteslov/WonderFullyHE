#include "m2424/profile_report.hpp"

#include <cmath>
#include <iostream>
#include <vector>

int main() {
    const std::vector<std::pair<std::string, m2424::CkksProfile>> profiles = {
        {"basic_ckks", {8192, {60, 40, 40, 60}, std::pow(2.0, 40), 4096}},
        {"depth_ckks", {16384, {60, 40, 40, 40, 40, 60}, std::pow(2.0, 40), 8192}}
    };

    std::cout << m2424::profile_report_csv_header() << '\n';
    for (const auto& entry : profiles) {
        std::cout << m2424::to_csv_row(m2424::describe_profile(entry.first, entry.second)) << '\n';
    }

    return 0;
}
