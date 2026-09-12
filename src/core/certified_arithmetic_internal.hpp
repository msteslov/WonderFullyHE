#pragma once
#include "m2424/evalround_execution.hpp"
namespace m2424 {
Cipher executeCertifiedArithmetic(SealAdapter&,const std::vector<Cipher>&,const std::vector<EvalRoundExecutionNode>&,const std::vector<Plain>&,std::size_t,const std::function<void(std::size_t,const Cipher&)>&);
}
