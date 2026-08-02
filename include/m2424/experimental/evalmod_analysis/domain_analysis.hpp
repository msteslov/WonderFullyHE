#pragma once

#include <cstddef>
#include <string>

namespace m2424::experimental {

enum class TailModel { Deterministic, Gaussian, Subgaussian };

struct TailModelProvenance {
    std::string derivation;
    std::string includedNoiseSources;
    std::string assumptions;
};

/// Параметр sigma является доказанным subgaussian parameter, а не standard deviation.
struct EvalModCiphertextModel {
    TailModel tailModel{TailModel::Deterministic};
    std::size_t coefficientCount{};
    double deterministicIntegerOffset{};
    double integerNoiseSubgaussianSigma{};
    double normalizedMessageAbsBound{};
    double normalizedEncodingErrorAbsBound{};
    double normalizedCoeffToSlotErrorAbsBound{};
    double relativePeriodMismatchAbsBound{};
    double additiveNormalizationErrorAbsBound{};
    double failureProbabilityLog2{-128.0};
    std::size_t analysisPrecisionBits{256};
    TailModelProvenance provenance;
};

struct EvalModDomainErrorBreakdown {
    double message{};
    double encoding{};
    double coeffToSlot{};
    double periodMismatch{};
    double additiveNormalization{};
};

struct EvalModDomain {
    std::size_t integerBound{};
    double normalizedResidualBound{};
    std::string normalizedResidualBoundDecimal;
    double failureProbabilityLog2{};
    EvalModDomainErrorBreakdown errors;

    /// Консервативный margin, округлённый только вниз.
    double discontinuityMargin() const;
};

/// MPFR estimate с направленным вверх округлением; proof зависит от корректности model sigma.
EvalModDomain estimateEvalModDomain(const EvalModCiphertextModel& model);

bool isEvalModDomainValid(const EvalModDomain& domain);

} // namespace m2424::experimental
