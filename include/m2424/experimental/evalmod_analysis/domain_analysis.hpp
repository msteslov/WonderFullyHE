#pragma once

#include <cstddef>
#include <cstdint>

namespace m2424::experimental {

/// Явный вероятностный контракт для целой части после ModRaise.
struct EvalModCiphertextModel {
    std::size_t coefficientCount{};
    double deterministicIntegerOffset{};
    double integerNoiseStddev{};
    double normalizedMessageAbsBound{};
    double normalizedEncodingErrorAbsBound{};
    double normalizedCoeffToSlotErrorAbsBound{};
    double normalizedScalePeriodErrorAbsBound{};
    double failureProbabilityLog2{-128.0};
};

struct EvalModDomain {
    std::size_t integerBound{};
    double normalizedResidualBound{};
    double discontinuityMargin{};
    double failureProbabilityLog2{};
};

/// Выводит K по union bound, а rho из всех консервативно нормированных bounds.
/// Нормировку на exact qSource вызывающая сторона должна выполнить без double.
/// Модель должна быть получена отдельно из параметров secret/error/ciphertext state.
EvalModDomain analyzeEvalModDomain(const EvalModCiphertextModel& model);

bool isEvalModDomainValid(const EvalModDomain& domain);

} // namespace m2424::experimental
