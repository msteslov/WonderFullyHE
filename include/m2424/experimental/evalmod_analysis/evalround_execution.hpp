#pragma once
#include "m2424/evalround_execution.hpp"
namespace m2424::experimental {
struct EvalRoundBinaryDigitSearchResult;
BootstrapBound evalRoundBackendKeyNoiseSupport();
EvalRoundDigitPolynomial certifyEvalRoundDigitPolynomial(const EvalRoundProblem&,std::size_t digit,
    const EvalModPolynomial&,std::size_t subdivisions=64);
struct EvalRoundExecutionOptions {
    // Bound relative to a real z in D_{K,rho}; supplied by the upstream stage.
    // PR-1 Unknown bounds remain Unknown here. PR-3 accepts deterministic inputs.
    BootstrapBound inputSemanticError;
    BootstrapBound evaluationKeyNoiseCoefficientSupport{evalRoundBackendKeyNoiseSupport()};
};
class EvalRoundExecutionCompiler {
public:
    // Verifies polynomial whole-domain evidence and compiles actual coefficients.
    // Never trusts reference local zeros or reference-only extraction targets.
    static EvalRoundExecutionPlan compile(SealAdapter&, const Cipher&, const EvalRoundPlan&,
        const EvalRoundExecutionOptions&);
    // Backend-aware bounded planner for every already generated/certified
    // record in the fixed K=64 binary search result.
    static EvalRoundExecutionPlan compile(SealAdapter&, const Cipher&,
        const EvalRoundBinaryDigitSearchResult&, const EvalRoundExecutionOptions&);
};
}
