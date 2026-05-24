#include "m2424/profile_report.hpp"
#include "m2424/profiles.hpp"

#include <iostream>

int main() {
    std::cout << m2424::profile_report_csv_header() << '\n';
    for (const auto& entry : m2424::profiles::all()) {
        std::cout << m2424::to_csv_row(m2424::describe_profile(entry.first, entry.second)) << '\n';
    }

    return 0;
}
