#pragma once

#include "m2424/abft.hpp"
#include "m2424/accuracy.hpp"
#include "m2424/bootstrap_candidates.hpp"
#include "m2424/canonical_embedding_reference.hpp"
#include "m2424/profiles.hpp"
#include "m2424/profile_report.hpp"
#include "m2424/seal_adapter.hpp"
#include "m2424/security_report.hpp"

namespace m2424 {

/// Возвращает строку версии библиотеки.
const char* version() noexcept;

} // namespace m2424
