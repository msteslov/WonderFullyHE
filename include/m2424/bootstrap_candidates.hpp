#pragma once

#include "m2424/accuracy.hpp"

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace m2424 {

/// Ограничения на ресурсы одного полного bootstrap-цикла.
struct BootstrapResourceLimits {
    double maxLatencyMs{};
    std::size_t maxRamMiB{};
    std::size_t maxEvaluationKeyMiB{};
};

/// Консервативное распределение итоговой ошибки между частями будущего цикла.
struct BootstrapErrorBudget {
    double encode{};
    double modRaise{};
    double coeffToSlot{};
    double evalMod{};
    double slotToCoeff{};
    double reserve{};
};

/// Проектный кандидат конфигурации bootstrap, ещё не являющийся CKKS-профилем исполнения.
struct BootstrapCandidate {
    std::string id;
    std::string description;
    std::size_t activeSlots{};
    std::size_t polyModulusDegree{};
    double targetScaleLog2{};
    double inputAbsBound{};
    int minimumSecurityBits{};
    BootstrapResourceLimits resources;
    BootstrapErrorBudget errorBudget;
};

/// Подтверждённые верхние границы ошибок компонентов для одного кандидата.
struct BootstrapComponentBounds {
    std::optional<double> encode;
    std::optional<double> modRaise;
    std::optional<double> coeffToSlot;
    std::optional<double> evalMod;
    std::optional<double> slotToCoeff;
};

/// Результат консервативной оценки возможности достичь целевой точности.
struct BootstrapFeasibilityForecast {
    bool complete{};
    bool passes{};
    double predictedTotalError{};
    std::string blocker;
};

/// Возвращает каталог кандидатов для экспериментов; ни один из них не означает готовую реализацию.
const std::vector<BootstrapCandidate>& bootstrapCandidates();

/// Возвращает кандидата по идентификатору или генерирует исключение для неизвестного идентификатора.
const BootstrapCandidate& bootstrapCandidateById(const std::string& id);

/// Проверяет корректность лимитов и распределения бюджета кандидата.
bool isBootstrapCandidateValid(const BootstrapCandidate& candidate);

/// Оценивает кандидата только по подтверждённым верхним границам всех компонентных ошибок.
BootstrapFeasibilityForecast forecastBootstrapFeasibility(const BootstrapCandidate& candidate,
                                                          const BootstrapComponentBounds& bounds);

} // namespace m2424
