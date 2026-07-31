#pragma once

#include "m2424/canonical_embedding_reference.hpp"
#include "m2424/coeff_to_slot_contract.hpp"
#include "m2424/seal_adapter.hpp"

#include <cstddef>
#include <memory>
#include <vector>

namespace m2424 {

struct CoeffToSlotResult {
    Cipher slotCipherFirst;
    Cipher slotCipherSecond;
};

struct CoeffToSlotFactorization {
    /// Число исходных разреженных факторов, объединённых на каждом ciphertext-уровне.
    std::vector<std::size_t> radices;
};

struct CoeffToSlotPlanMetrics {
    std::size_t butterflyStages{};
    std::size_t permutationStages{};
    std::size_t depth{};
    std::vector<std::size_t> diagonalsPerStage;
    std::size_t plaintextMultiplicationsPerApply{};
    std::size_t rotationsPerApply{};
    std::size_t additionsPerApply{};
    std::size_t rescalesPerApply{};
    std::size_t uniqueEvaluationKeys{};
    std::size_t storedComplexValues{};
};

struct CoeffToSlotFactorizationEstimate {
    CoeffToSlotFactorization factorization;
    std::size_t estimatedDiagonals{};
    std::size_t estimatedRotationsPerOutput{};
    std::size_t estimatedUniqueRotations{};
};

class PreparedCoeffToSlotPlan {
public:
    PreparedCoeffToSlotPlan(PreparedCoeffToSlotPlan&&) noexcept;
    PreparedCoeffToSlotPlan& operator=(PreparedCoeffToSlotPlan&&) noexcept;
    ~PreparedCoeffToSlotPlan();
    PreparedCoeffToSlotPlan(const PreparedCoeffToSlotPlan&) = delete;
    PreparedCoeffToSlotPlan& operator=(const PreparedCoeffToSlotPlan&) = delete;

    std::size_t plaintextCount() const;
    std::size_t serializedPlaintextBytes(const SealAdapter& adapter) const;

private:
    struct Impl;
    explicit PreparedCoeffToSlotPlan(std::unique_ptr<Impl> implementation);
    std::unique_ptr<Impl> pimpl_;
    friend class CoeffToSlotPlan;
    friend class CoeffToSlot;
};

/**
 * Подготовленный математический план canonical inverse embedding в порядке
 * слотов SEAL. План хранит только разреженные факторизованные диагонали и
 * никогда не материализует плотную матрицу U_0^T/N.
 */
class CoeffToSlotPlan {
public:
    explicit CoeffToSlotPlan(std::size_t polyModulusDegree,
                             std::size_t targetDepth = 4);
    CoeffToSlotPlan(std::size_t polyModulusDegree,
                    CoeffToSlotFactorization factorization);
    ~CoeffToSlotPlan();
    CoeffToSlotPlan(CoeffToSlotPlan&&) noexcept;
    CoeffToSlotPlan& operator=(CoeffToSlotPlan&&) noexcept;
    CoeffToSlotPlan(const CoeffToSlotPlan&) = delete;
    CoeffToSlotPlan& operator=(const CoeffToSlotPlan&) = delete;

    std::size_t polyModulusDegree() const;
    std::size_t butterflyStageCount() const;
    std::size_t rawStageCount() const;
    std::size_t depth() const;
    const CoeffToSlotFactorization& factorization() const;
    std::vector<std::vector<int>> stageRotationSteps() const;
    CoeffToSlotPlanRequirements requirements() const;
    CoeffToSlotPlanMetrics metrics() const;

    /// Кодирует все stage-диагонали после обязательного preflight.
    PreparedCoeffToSlotPlan prepare(SealAdapter& adapter,
                                    const RaisedCipher& input,
                                    const CoeffToSlotContract& contract) const;
    /// Проверяет context fingerprint, parms_id, уровень, scale и контракт подготовки.
    bool isPreparedFor(const PreparedCoeffToSlotPlan& prepared,
                       const SealAdapter& adapter,
                       const RaisedCipher& input,
                       const CoeffToSlotContract& contract) const;

    /// Независимое plaintext-исполнение факторизации для математических тестов.
    std::pair<ComplexVector, ComplexVector> applyPlain(const ComplexVector& slots) const;

    static std::vector<CoeffToSlotFactorizationEstimate> rankFactorizations(
        std::size_t polyModulusDegree,
        std::size_t depth,
        std::size_t limit = 10);

private:
    struct Impl;
    std::unique_ptr<Impl> pimpl_;
    friend class CoeffToSlot;
};

/**
 * Canonical bootstrap CoeffToSlot для полного CKKS packing.
 *
 * Потребляет только RaisedCipher и возвращает две половины коэффициентов.
 * Каждый ciphertext-уровень исполняет один сгруппированный разреженный
 * FFT-фактор через rotations, plaintext multiplication и сложения.
 */
class CoeffToSlot {
public:
    explicit CoeffToSlot(std::size_t polyModulusDegree,
                         std::size_t targetDepth = 4);
    CoeffToSlot(std::size_t polyModulusDegree,
                CoeffToSlotFactorization factorization);
    ~CoeffToSlot();
    CoeffToSlot(CoeffToSlot&&) noexcept;
    CoeffToSlot& operator=(CoeffToSlot&&) noexcept;
    CoeffToSlot(const CoeffToSlot&) = delete;
    CoeffToSlot& operator=(const CoeffToSlot&) = delete;

    const CoeffToSlotPlan& plan() const;
    CoeffToSlotPlanRequirements requirements() const;
    PreparedCoeffToSlotPlan prepare(SealAdapter& adapter,
                                    const RaisedCipher& input,
                                    const CoeffToSlotContract& contract) const;
    CoeffToSlotResult apply(SealAdapter& adapter,
                            RaisedCipher&& input,
                            const CoeffToSlotContract& contract,
                            const PreparedCoeffToSlotPlan& prepared) const;

private:
    CoeffToSlotPlan plan_;
};

} // namespace m2424
