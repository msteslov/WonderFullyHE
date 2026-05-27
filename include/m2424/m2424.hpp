#pragma once

#include "m2424/abft.hpp"
#include "m2424/accuracy.hpp"
#include "m2424/bootstrap.hpp"
#include "m2424/bootstrap_plan.hpp"
#include "m2424/bootstrap_prototype.hpp"
#include "m2424/bootstrap_scaling.hpp"
#include "m2424/checked_evaluator.hpp"
#include "m2424/diagonal_transform.hpp"
#include "m2424/eval_mod.hpp"
#include "m2424/linear_transform.hpp"
#include "m2424/polynomial.hpp"
#include "m2424/profiles.hpp"
#include "m2424/profile_report.hpp"
#include "m2424/seal_adapter.hpp"
#include "m2424/security_report.hpp"

namespace m2424 {

// Returns a semantic-ish version string of the library.
const char* version() noexcept;

} // namespace m2424
