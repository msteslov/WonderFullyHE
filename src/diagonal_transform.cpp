#include "m2424/diagonal_transform.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace m2424 {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

bool is_power_of_two(std::size_t value) {
    return value != 0 && (value & (value - 1)) == 0;
}

std::uint64_t pow_mod(std::uint64_t base, std::uint64_t exp, std::uint64_t mod) {
    std::uint64_t result = 1;
    base %= mod;
    while (exp != 0) {
        if (exp & 1U) {
            result = (result * base) % mod;
        }
        base = (base * base) % mod;
        exp >>= 1U;
    }
    return result;
}

void validate_square_matrix(const ComplexMatrix& matrix) {
    if (matrix.empty()) {
        throw std::invalid_argument("matrix must not be empty");
    }
    const std::size_t n = matrix.size();
    for (const auto& row : matrix) {
        if (row.size() != n) {
            throw std::invalid_argument("matrix must be square");
        }
        for (const auto& value : row) {
            if (!std::isfinite(value.real()) || !std::isfinite(value.imag())) {
                throw std::invalid_argument("matrix values must be finite");
            }
        }
    }
}

bool has_nonzero(const ComplexVector& values, double tolerance) {
    return std::any_of(values.begin(), values.end(), [tolerance](const Complex& value) {
        return std::abs(value) > tolerance;
    });
}

std::size_t baby_step_count(std::size_t dimension) {
    return static_cast<std::size_t>(std::ceil(std::sqrt(static_cast<double>(dimension))));
}

ComplexVector place_diagonal_for_giant_step(const ComplexVector& diagonal,
                                            std::size_t giant_step,
                                            std::size_t slot_count) {
    if (giant_step + diagonal.size() > slot_count) {
        throw std::invalid_argument("giant step exceeds available CKKS slots");
    }
    ComplexVector shifted(slot_count, Complex{0.0, 0.0});
    for (std::size_t i = 0; i < diagonal.size(); ++i) {
        shifted[i + giant_step] = diagonal[i];
    }
    return shifted;
}

Cipher pairwise_add(SealAdapter& adapter, std::vector<Cipher> terms) {
    if (terms.empty()) {
        throw std::invalid_argument("diagonal transform has no non-zero terms");
    }
    while (terms.size() > 1) {
        std::vector<Cipher> next;
        next.reserve((terms.size() + 1) / 2);
        for (std::size_t i = 0; i < terms.size(); i += 2) {
            if (i + 1 < terms.size()) {
                next.push_back(adapter.add(terms[i], terms[i + 1]));
            } else {
                next.push_back(std::move(terms[i]));
            }
        }
        terms = std::move(next);
    }
    return std::move(terms.front());
}

} // namespace

DiagonalLinearTransform::DiagonalLinearTransform(std::vector<DiagonalTerm> terms)
    : terms_(std::move(terms)) {
    if (terms_.empty()) {
        throw std::invalid_argument("diagonal transform must contain at least one term");
    }
    dimension_ = terms_.front().diagonal.size();
    if (dimension_ == 0) {
        throw std::invalid_argument("diagonal transform terms must not be empty");
    }
    for (const auto& term : terms_) {
        if (term.rotation < 0) {
            throw std::invalid_argument("diagonal transform rotations must be non-negative");
        }
        if (term.diagonal.size() != dimension_) {
            throw std::invalid_argument("diagonal transform terms must have equal size");
        }
        for (const auto& value : term.diagonal) {
            if (!std::isfinite(value.real()) || !std::isfinite(value.imag())) {
                throw std::invalid_argument("diagonal transform coefficients must be finite");
            }
        }
    }
}

DiagonalLinearTransform DiagonalLinearTransform::from_matrix(const ComplexMatrix& matrix, double zero_tolerance) {
    validate_square_matrix(matrix);
    if (zero_tolerance < 0.0 || !std::isfinite(zero_tolerance)) {
        throw std::invalid_argument("zero_tolerance must be a non-negative finite value");
    }

    const std::size_t n = matrix.size();
    std::vector<DiagonalTerm> terms;
    terms.reserve(n);
    for (std::size_t rotation = 0; rotation < n; ++rotation) {
        ComplexVector diagonal;
        diagonal.reserve(n);
        for (std::size_t row = 0; row < n; ++row) {
            diagonal.push_back(matrix[row][(row + rotation) % n]);
        }
        if (has_nonzero(diagonal, zero_tolerance)) {
            terms.push_back(DiagonalTerm{static_cast<int>(rotation), std::move(diagonal)});
        }
    }
    return DiagonalLinearTransform(std::move(terms));
}

const std::vector<DiagonalTerm>& DiagonalLinearTransform::terms() const noexcept {
    return terms_;
}

std::vector<int> DiagonalLinearTransform::rotation_steps() const {
    std::vector<int> steps;
    const std::size_t baby = baby_step_count(dimension_);
    for (const auto& term : terms_) {
        const std::size_t rotation = static_cast<std::size_t>(term.rotation);
        const std::size_t baby_rotation = rotation % baby;
        const std::size_t giant_rotation = rotation - baby_rotation;
        if (baby_rotation != 0) {
            steps.push_back(static_cast<int>(baby_rotation));
        }
        if (giant_rotation != 0) {
            steps.push_back(static_cast<int>(giant_rotation));
        }
    }
    std::sort(steps.begin(), steps.end());
    steps.erase(std::unique(steps.begin(), steps.end()), steps.end());
    return steps;
}

ComplexVector DiagonalLinearTransform::apply_plain(const ComplexVector& input) const {
    if (input.size() != dimension_) {
        throw std::invalid_argument("input size must match diagonal transform dimension");
    }
    ComplexVector result(dimension_, Complex{0.0, 0.0});
    for (const auto& term : terms_) {
        const std::size_t rotation = static_cast<std::size_t>(term.rotation);
        for (std::size_t i = 0; i < dimension_; ++i) {
            result[i] += term.diagonal[i] * input[(i + rotation) % dimension_];
        }
    }
    return result;
}

Cipher DiagonalLinearTransform::apply(SealAdapter& adapter, const Cipher& input) const {
    const std::size_t baby = baby_step_count(dimension_);
    std::vector<Cipher> baby_rotations(baby);
    std::vector<bool> baby_ready(baby, false);
    baby_rotations[0] = input;
    baby_ready[0] = true;

    for (const auto& term : terms_) {
        const std::size_t baby_rotation = static_cast<std::size_t>(term.rotation) % baby;
        if (!baby_ready[baby_rotation]) {
            baby_rotations[baby_rotation] = adapter.rotate(input, static_cast<int>(baby_rotation));
            baby_ready[baby_rotation] = true;
        }
    }

    const std::size_t giant_count = (dimension_ + baby - 1) / baby;
    std::vector<std::vector<const DiagonalTerm*>> groups(giant_count);
    for (const auto& term : terms_) {
        groups[static_cast<std::size_t>(term.rotation) / baby].push_back(&term);
    }

    std::vector<Cipher> giant_terms;
    giant_terms.reserve(giant_count);
    for (std::size_t giant_index = 0; giant_index < groups.size(); ++giant_index) {
        if (groups[giant_index].empty()) {
            continue;
        }

        const std::size_t giant_rotation = giant_index * baby;
        std::vector<Cipher> inner_terms;
        inner_terms.reserve(groups[giant_index].size());
        for (const auto* term : groups[giant_index]) {
            const std::size_t baby_rotation = static_cast<std::size_t>(term->rotation) % baby;
            const auto shifted_diagonal = place_diagonal_for_giant_step(term->diagonal,
                                                                        giant_rotation,
                                                                        adapter.slot_count());
            Plain encoded_diagonal = adapter.encode_complex_like(shifted_diagonal, baby_rotations[baby_rotation]);
            inner_terms.push_back(adapter.mul_plain_rescale(baby_rotations[baby_rotation], encoded_diagonal));
        }

        Cipher inner = pairwise_add(adapter, std::move(inner_terms));
        if (giant_rotation == 0) {
            giant_terms.push_back(std::move(inner));
        } else {
            giant_terms.push_back(adapter.rotate(inner, static_cast<int>(giant_rotation)));
        }
    }

    return pairwise_add(adapter, std::move(giant_terms));
}

ComplexMatrix canonical_embedding_matrix(std::size_t slots) {
    if (!is_power_of_two(slots)) {
        throw std::invalid_argument("slots must be a non-zero power of two");
    }

    const std::uint64_t polynomial_degree = static_cast<std::uint64_t>(2 * slots);
    const std::uint64_t cyclotomic_order = 2 * polynomial_degree;
    const Complex zeta = std::exp(Complex{0.0, 2.0 * kPi / static_cast<double>(cyclotomic_order)});

    ComplexMatrix matrix(slots, ComplexVector(slots));
    for (std::size_t i = 0; i < slots; ++i) {
        const std::uint64_t g = pow_mod(5, static_cast<std::uint64_t>(i), cyclotomic_order);
        for (std::size_t j = 0; j < slots; ++j) {
            matrix[i][j] = std::pow(zeta, static_cast<double>(g * j));
        }
    }
    return matrix;
}

ComplexMatrix invert_matrix(const ComplexMatrix& matrix) {
    validate_square_matrix(matrix);
    const std::size_t n = matrix.size();
    ComplexMatrix a(n, ComplexVector(2 * n, Complex{0.0, 0.0}));

    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            a[i][j] = matrix[i][j];
        }
        a[i][n + i] = Complex{1.0, 0.0};
    }

    for (std::size_t col = 0; col < n; ++col) {
        std::size_t pivot = col;
        for (std::size_t row = col + 1; row < n; ++row) {
            if (std::abs(a[row][col]) > std::abs(a[pivot][col])) {
                pivot = row;
            }
        }
        if (std::abs(a[pivot][col]) < 1e-14) {
            throw std::invalid_argument("matrix is singular");
        }
        if (pivot != col) {
            std::swap(a[pivot], a[col]);
        }

        const Complex divisor = a[col][col];
        for (std::size_t j = 0; j < 2 * n; ++j) {
            a[col][j] /= divisor;
        }

        for (std::size_t row = 0; row < n; ++row) {
            if (row == col) {
                continue;
            }
            const Complex factor = a[row][col];
            if (std::abs(factor) == 0.0) {
                continue;
            }
            for (std::size_t j = 0; j < 2 * n; ++j) {
                a[row][j] -= factor * a[col][j];
            }
        }
    }

    ComplexMatrix inverse(n, ComplexVector(n));
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            inverse[i][j] = a[i][n + j];
        }
    }
    return inverse;
}

} // namespace m2424
