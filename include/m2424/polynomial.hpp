#pragma once

#include "m2424/seal_adapter.hpp"

#include <cstddef>
#include <vector>

namespace m2424 {

struct PolynomialTerm {
    std::size_t degree{};
    double coefficient{};
};

class PolynomialEvaluator {
public:
    explicit PolynomialEvaluator(std::vector<PolynomialTerm> terms);

    const std::vector<PolynomialTerm>& terms() const noexcept;
    Cipher evaluate(SealAdapter& adapter, const Cipher& input) const;

private:
    std::vector<PolynomialTerm> terms_;
};

} // namespace m2424
