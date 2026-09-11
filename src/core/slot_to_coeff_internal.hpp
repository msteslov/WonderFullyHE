#pragma once
#include "m2424/slot_to_coeff.hpp"
#include <map>
namespace m2424 {
// Every nonzero entry is exactly a power of zeta. -1 denotes exact zero.
using RootDiagonalMap=std::map<std::size_t,std::vector<int>>;
struct SlotToCoeffPlan::Impl {
    std::size_t degree{};
    SlotToCoeffFactorization factorization;
    std::vector<RootDiagonalMap> factors[2];
    std::vector<std::size_t> babySteps;
    SlotToCoeffRequirements keys;
    SlotToCoeffMetrics metrics;
};
struct StCPreparedTerm { int baby{}; Plain diagonal; BootstrapBound perturbation; };
struct StCPreparedGroup { int giant{}; std::vector<StCPreparedTerm> terms; };
struct PreparedSlotToCoeffPlan::Impl {
    std::size_t degree{},chainIndex{};
    std::size_t inputComponents[2]{};
    SlotToCoeffFactorization factorization;
    std::array<std::uint64_t,4> fingerprint{};
    std::vector<std::uint64_t> primes;
    std::vector<std::vector<StCPreparedGroup>> factors[2];
    SlotToCoeffCertificate certificate;
};
}
