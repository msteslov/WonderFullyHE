#include "m2424/polynomial.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace m2424 {

static bool is_zero(double value) {
    return std::abs(value) < 1e-15;
}

PolynomialEvaluator::PolynomialEvaluator(std::vector<PolynomialTerm> terms) : terms_(std::move(terms)) {
    if (terms_.empty()) {
        throw std::invalid_argument("polynomial must contain at least one term");
    }
    std::sort(terms_.begin(), terms_.end(), [](const PolynomialTerm& a, const PolynomialTerm& b) {
        return a.degree < b.degree;
    });
    for (std::size_t i = 0; i < terms_.size(); ++i) {
        if (terms_[i].degree == 0) {
            throw std::invalid_argument("constant-only terms are not supported without an input ciphertext term");
        }
        if (i > 0 && terms_[i - 1].degree == terms_[i].degree) {
            throw std::invalid_argument("polynomial degrees must be unique");
        }
        if (!std::isfinite(terms_[i].coefficient)) {
            throw std::invalid_argument("polynomial coefficients must be finite");
        }
    }
}

const std::vector<PolynomialTerm>& PolynomialEvaluator::terms() const noexcept {
    return terms_;
}

Cipher PolynomialEvaluator::evaluate(SealAdapter& adapter, const Cipher& input) const {
    bool has_result = false;
    Cipher result;
    Cipher power = input;
    std::size_t current_degree = 1;

    for (const auto& term : terms_) {
        // A zero coefficient has no observable contribution. Skipping it before
        // exponentiation avoids needless ciphertext multiplications and rescales.
        if (is_zero(term.coefficient)) {
            continue;
        }

        while (current_degree < term.degree) {
            Cipher base = adapter.modSwitchTo(input, power);
            power = adapter.rescaleToNext(adapter.relinearize(adapter.multiply(power, base)));
            ++current_degree;
        }

        Plain coefficient = adapter.encodeScalarFor(term.coefficient, power);
        Cipher weighted = adapter.rescaleToNext(adapter.multiplyPlain(power, coefficient));

        if (!has_result) {
            result = std::move(weighted);
            has_result = true;
        } else {
            result = adapter.alignForAddition(result, weighted);
            result = adapter.add(result, weighted);
        }
    }

    if (!has_result) {
        throw std::invalid_argument("polynomial contains only zero coefficients");
    }
    return result;
}

} // namespace m2424
