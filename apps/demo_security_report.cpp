#include "m2424/profiles.hpp"
#include "m2424/security_report.hpp"

#include <iostream>
#include <vector>

int main() {
    std::vector<m2424::SecurityReport> reports;
    const auto profiles = m2424::profiles::all();
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
