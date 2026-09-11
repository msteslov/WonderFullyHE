#pragma once
#include "m2424/bootstrap_contract.hpp"
namespace m2424 {
struct LinearTransformScale { std::string numerator,denominator; std::uint64_t binary64Bits{}; };
struct LinearTransformFactorBound {
    BootstrapBound kappa,delta,localArithmeticError;
};
struct LinearTransformPropagation {
    BootstrapContractResult result;
    BootstrapBound idealMagnitude,propagatedError,operatorPerturbationError,semanticError;
};
LinearTransformPropagation propagateLinearTransform(const BootstrapBound& magnitude,
    const BootstrapBound& error,const LinearTransformFactorBound&);
BootstrapBound linearTransformGain(const std::vector<LinearTransformFactorBound>&);
}
