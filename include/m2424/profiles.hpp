#pragma once

#include "m2424/seal_adapter.hpp"

#include <string>
#include <utility>
#include <vector>

namespace m2424::profiles {

CkksProfile fast_demo_ckks();
CkksProfile basic_ckks();
CkksProfile balanced_ckks();
CkksProfile depth_ckks();
CkksProfile high_precision_ckks();
CkksProfile boot_ckks();

std::vector<std::pair<std::string, CkksProfile>> all();
CkksProfile by_name(const std::string& name);

} // namespace m2424::profiles
