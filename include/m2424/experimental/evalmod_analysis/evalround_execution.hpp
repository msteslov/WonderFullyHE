#pragma once
#include "m2424/evalround_execution.hpp"
namespace m2424::experimental {
BootstrapBound evalRoundBackendKeyNoiseSupport();
struct EvalRoundExecutionOptions {
    // Bound relative to a real z in D_{K,rho}; supplied by the upstream stage.
    // PR-1 Unknown bounds remain Unknown here. PR-3 accepts deterministic inputs.
    BootstrapBound inputSemanticError;
    BootstrapBound evaluationKeyNoiseCoefficientSupport{evalRoundBackendKeyNoiseSupport()};
};
class EvalRoundExecutionCompiler {
public:
    // Rebuilds the known polynomial certificate. Never trusts reference local zeros.
    static EvalRoundExecutionPlan compile(SealAdapter&, const Cipher&, const EvalRoundPlan&,
        const EvalRoundExecutionOptions&);
};
}
