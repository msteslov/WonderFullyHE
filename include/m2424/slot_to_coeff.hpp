#pragma once
#include "m2424/canonical_embedding_reference.hpp"
#include "m2424/linear_transform_contract.hpp"
#include <functional>
#include <memory>
namespace m2424 {
struct SlotToCoeffFactorization { std::vector<std::size_t> radices; };
struct SlotToCoeffRequirements { std::vector<int> rotations; bool relinearization{}; };
struct SlotToCoeffMetrics {
    std::size_t depth{},rotations{},rescales{},innerModDowns{},finalModDowns{},plaintexts{},relinearizations{};
};
struct SlotToCoeffContract {
    BootstrapBound inputMagnitude[2],inputError[2];
    BootstrapBound evaluationKeyNoiseSupport{21,BootstrapBoundKind::Deterministic,"SEAL CBD coefficient support 21",{}};
};
struct SlotToCoeffRuntimeStage {
    std::string operation;
    std::size_t branch{},factor{},chainIndex{};
    std::vector<std::uint64_t> activePrimes;
    LinearTransformScale inputScale,outputScale,arithmeticScale;
    BootstrapBound magnitude,semanticError,localError,scaleRepresentationError;
    std::string centeredHeadroomNumerator,headroomProvenance;
};
struct SlotToCoeffFactorTrace {
    std::size_t branch{},factor{};
    LinearTransformFactorBound bounds;
    LinearTransformPropagation propagation;
    std::vector<SlotToCoeffRuntimeStage> runtime;
    // Analytic headroom checks inside the atomic double-hoisted kernel.
    std::vector<std::string> internalHeadroomProofs;
};
struct SlotToCoeffCertificate {
    BootstrapContractResult result;
    std::vector<SlotToCoeffFactorTrace> factors;
    std::vector<SlotToCoeffRuntimeStage> inputStages;
    BootstrapBound outputError,gamma;
    SlotToCoeffMetrics metrics;
    SlotToCoeffRequirements keys;
    LinearTransformScale inputScales[2],outputScale;
};
class SlotToCoeffPlan;
class PreparedSlotToCoeffPlan {
public:
    PreparedSlotToCoeffPlan();
    const SlotToCoeffCertificate& certificate() const;
private:
    struct Impl; std::shared_ptr<const Impl> impl_;
    friend class SlotToCoeffPlan;
};
struct SlotToCoeffResult { Cipher ciphertext; SlotToCoeffCertificate certificate; };
class SlotToCoeffPlan {
public:
    explicit SlotToCoeffPlan(std::size_t degree,std::size_t depth=4);
    SlotToCoeffPlan(std::size_t degree,SlotToCoeffFactorization);
    std::size_t polyModulusDegree() const;
    const SlotToCoeffFactorization& factorization() const;
    SlotToCoeffRequirements requirements() const;
    SlotToCoeffMetrics metrics() const;
    ComplexVector applyPlain(const ComplexVector& first,const ComplexVector& second) const;
    // Exposes the same factorized mathematical path for per-factor diagnostics.
    std::vector<ComplexVector> applyPlainTrace(const ComplexVector&,std::size_t branch) const;
    // Rigorous preparation is implemented by optional m2424_evalmod_analysis.
    // Core plan/plaintext execution/preflight/apply have no GMP/MPFR dependency.
    PreparedSlotToCoeffPlan prepare(SealAdapter&,const Cipher&,const Cipher&,const SlotToCoeffContract&) const;
    BootstrapContractResult preflight(const SealAdapter&,const Cipher&,const Cipher&,const PreparedSlotToCoeffPlan&) const;
    SlotToCoeffResult apply(SealAdapter&,const Cipher&,const Cipher&,const PreparedSlotToCoeffPlan&,
        const std::function<void(const SlotToCoeffRuntimeStage&,const Cipher&)>& observer={}) const;
private:
    struct Impl; std::shared_ptr<const Impl> impl_;
};
// Same plan -> prepare -> preflight -> apply API under the operation name.
using SlotToCoeff = SlotToCoeffPlan;
}
