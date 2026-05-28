#pragma once

#include "m2424/seal_adapter.hpp"

#include <complex>
#include <cstddef>
#include <string>
#include <unordered_map>
#include <vector>

namespace m2424 {

using Complex = std::complex<double>;
using ComplexVector = std::vector<Complex>;
using ComplexMatrix = std::vector<ComplexVector>;

struct DiagonalTerm {
    int rotation{};
    ComplexVector diagonal;
};

class DiagonalLinearTransform {
public:
    explicit DiagonalLinearTransform(std::vector<DiagonalTerm> terms);

    static DiagonalLinearTransform from_matrix(const ComplexMatrix& matrix, double zero_tolerance = 0.0);

    const std::vector<DiagonalTerm>& terms() const noexcept;
    std::vector<int> rotation_steps() const;
    ComplexVector apply_plain(const ComplexVector& input) const;
    Cipher apply(SealAdapter& adapter, const Cipher& input) const;
    Cipher apply_at_plain_scale(SealAdapter& adapter, const Cipher& input, double plain_scale) const;

private:
    std::vector<DiagonalTerm> terms_;
    std::size_t dimension_{};
    mutable std::unordered_map<std::string, std::vector<Plain>> encoded_diagonal_cache_;
};

ComplexMatrix canonical_embedding_matrix(std::size_t slots);
ComplexMatrix invert_matrix(const ComplexMatrix& matrix);

} // namespace m2424
