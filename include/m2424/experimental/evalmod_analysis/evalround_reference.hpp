#pragma once

#include "m2424/evalround.hpp"
#include "m2424/experimental/evalmod_analysis/exact_modular_oracle.hpp"

namespace m2424::experimental {

struct EvalRoundExactComplex { mpq_class real, imag; };
/// Exact domain membership and piecewise digit target. No floating-point rounding
/// decides the center or endpoint membership; rho is its exact binary64 dyadic.
ExactInteger evalRoundReferenceCenter(const EvalRoundProblem&, const mpq_class& z);
std::vector<int> extractEvalRoundReferenceDigits(const EvalRoundProblem&, EvalRoundRadix, const mpq_class& z);
ExactInteger reconstructEvalRoundDigitsExact(const std::vector<int>&, std::uint32_t K, EvalRoundRadix);
mpq_class evalRoundBinaryCleanerExact(const mpq_class&);
EvalRoundExactComplex evalRoundTernaryCleanerExact(const EvalRoundExactComplex&);
EvalRoundExactComplex evalRoundRootReference(int trit, std::size_t precisionBits = 256);

struct EvalRoundReferenceTrace {
    ExactInteger center;
    std::vector<int> targetDigits;
    std::vector<std::vector<EvalRoundExactComplex>> digitStages;
    mpq_class reconstructed;
    std::size_t precisionBits{};
};
/// Exact rational polynomials, MPFR phase/root evaluation, and rounding to the
/// requested reference precision after each stage. These values are diagnostics,
/// NEVER analytical error bounds. Concrete external polynomials can be evaluated
/// as diagnostics; this evaluation does not prove their digit approximation.
EvalRoundReferenceTrace evaluateEvalRoundReference(const EvalRoundProblem&, const EvalRoundCandidate&,
    const mpq_class& z, const std::vector<std::size_t>& cleaningCounts, std::size_t precisionBits = 256);

} // namespace m2424::experimental
