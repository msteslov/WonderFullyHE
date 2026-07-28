#pragma once

#include "m2424/seal_adapter.hpp"

#include <cstddef>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace m2424 {

/// Тип стадии исследовательского bootstrap-конвейера.
enum class BootstrapStageKind {
    ModUp,
    CoeffToSlot,
    EvalMod,
    SlotToCoeff
};

/// Требования стадии к evaluation keys.
struct BootstrapKeyRequirements {
    bool requiresRelin{};
    std::vector<int> rotationSteps;
};

/// Настройки сбора измерений при запуске pipeline.
struct BootstrapMeasurementConfig {
    bool measureDuration{true};
};

/// Метрики одной выполненной стадии.
struct BootstrapStageReport {
    BootstrapStageKind kind{};
    std::string name;
    CipherInfo input;
    CipherInfo output;
    double durationMs{};
};

/// Результат выполнения одной стадии.
struct BootstrapStageResult {
    Cipher ciphertext;
};

/**
 * Контракт взаимозаменяемой bootstrap-стадии.
 *
 * Реализация обязана объявить необходимые keys и не должна выполнять
 * неявные операции вне своей математической модели. До включения в рабочий
 * experiment runner она должна иметь plaintext/reference-тесты.
 */
class BootstrapStage {
public:
    virtual ~BootstrapStage() = default;

    /// Возвращает назначение стадии в конвейере.
    virtual BootstrapStageKind kind() const noexcept = 0;
    /// Возвращает стабильное имя реализации для отчётов.
    virtual std::string_view name() const noexcept = 0;
    /// Возвращает минимальный набор evaluation keys для стадии.
    virtual BootstrapKeyRequirements keyRequirements() const = 0;
    /// Применяет стадию к ciphertext без скрытого управления pipeline.
    virtual BootstrapStageResult execute(SealAdapter& adapter, const Cipher& input) const = 0;
};

/// Итог запуска pipeline вместе с метриками всех стадий.
struct BootstrapPipelineResult {
    Cipher ciphertext;
    std::vector<BootstrapStageReport> stages;
};

/**
 * Последовательно запускает выбранные стратегии bootstrap.
 *
 * Класс не выбирает алгоритмы и не меняет scale-policy: это обязанность
 * конкретных стадий и будущих planner-ов.
 */
class BootstrapPipeline {
public:
    explicit BootstrapPipeline(std::vector<std::unique_ptr<BootstrapStage>> stages);

    /// Объединяет требования всех стадий и удаляет дублирующиеся rotation steps.
    BootstrapKeyRequirements keyRequirements() const;
    /// Выполняет стадии в заданном порядке и формирует отчёт.
    BootstrapPipelineResult run(SealAdapter& adapter,
                                const Cipher& input,
                                const BootstrapMeasurementConfig& measurements = {}) const;

private:
    std::vector<std::unique_ptr<BootstrapStage>> stages_;
};

} // namespace m2424
