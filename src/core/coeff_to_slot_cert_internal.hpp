#pragma once
#include "slot_to_coeff_internal.hpp"
#include "m2424/evalround_plus_coeff_to_slot.hpp"
namespace m2424 {
struct CertifiedEvalRoundPlusCoeffToSlot::Impl {
    BootstrapInputContext source;
    CoeffToSlotContract hpContract,lpContract;
    CoeffToSlotFactorization factorization;
    std::shared_ptr<const PreparedRootLinearTransform> branches[2];
    CoeffToSlotBranchTrace traces[2];
    CoeffToSlotDomainCertificate domain;
};
}
