#include "m2424/security_report.hpp"

#include <cmath>
#include <iostream>
#include <utility>
#include <vector>

int main() {
    const std::vector<std::pair<std::string, m2424::CkksProfile>> profiles = {
        {"basic_ckks", {8192, {60, 40, 40, 60}, std::pow(2.0, 40), 4096}},
        {"depth_ckks", {16384, {60, 40, 40, 40, 40, 60}, std::pow(2.0, 40), 8192}}
    };

    std::vector<m2424::SecurityReport> reports;
    reports.reserve(profiles.size());

    std::cout << m2424::security_report_csv_header() << '\n';
    for (const auto& profile : profiles) {
        reports.push_back(m2424::analyze_security(profile.first, profile.second));
        std::cout << m2424::to_csv_row(reports.back()) << '\n';
    }

    std::cout << "project_minimum_security,"
              << m2424::to_string(m2424::project_minimum_security(reports))
              << '\n';

    return 0;
}
