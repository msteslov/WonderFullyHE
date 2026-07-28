#pragma once

#include "m2424/seal_adapter.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace m2424 {

/**
 * @brief Тип стадии исследовательского bootstrap-конвейера.
 *
 * Значение описывает математическую роль стадии, а не конкретный алгоритм.
 * Например, одна и та же стадия CoeffToSlot может иметь реализации DenseDiagonal,
 * BSGS или FFT-like.
 */
enum class BootstrapStageKind {
    /// Расширение modulus chain ciphertext.
    ModUp,
    /// Переход из коэффициентного представления в слоты.
    CoeffToSlot,
    /// Аппроксимация модульной функции.
    EvalMod,
    /// Обратный переход из слотов в коэффициенты.
    SlotToCoeff
};

/**
 * @brief Минимальные требования стадии к evaluation keys.
 *
 * Pipeline объединяет требования всех стадий до запуска. По полученному
 * набору вызывающий код генерирует или загружает ключи в SealAdapter.
 */
struct BootstrapKeyRequirements {
    /// Нужны ли relinearization keys.
    bool requiresRelin{};
    /// Шаги ротации, для которых требуются Galois keys.
    std::vector<int> rotationSteps;
};

/**
 * @brief Настройки сбора измерений при запуске pipeline.
 *
 * Измерения не влияют на математику стадий. Они нужны только для сравнения
 * реализаций в experiment runner.
 */
struct BootstrapMeasurementConfig {
    /// Измерять ли wall-clock время выполнения каждой стадии.
    bool measureDuration{true};
};

/**
 * @brief Отчёт об одном вызове BootstrapStage.
 *
 * Отчёт содержит фактические метаданные ciphertext до и после стадии, а не
 * прогноз planner-а. Благодаря этому можно сопоставлять план с исполнением.
 */
struct BootstrapStageReport {
    /// Тип выполненной стадии.
    BootstrapStageKind kind{};
    /// Стабильное имя реализации стадии.
    std::string name;
    /// Метаданные ciphertext перед запуском стадии.
    CipherInfo input;
    /// Метаданные ciphertext после завершения стадии.
    CipherInfo output;
    /// Изменение индекса modulus chain: положительное значение означает подъём, отрицательное — расход уровней.
    std::int64_t chainIndexDelta{};
    /// Измеренное время выполнения; равно нулю при отключённом измерении.
    double durationMs{};
};

/**
 * @brief Результат выполнения одной стадии.
 *
 * Стадия владеет возвращаемым ciphertext. Pipeline передаёт его следующей
 * стадии без изменения уровня или scale от своего имени.
 */
struct BootstrapStageResult {
    /// Ciphertext, который передаётся следующей стадии.
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

    /// @brief Возвращает математическую роль стадии в конвейере.
    virtual BootstrapStageKind kind() const noexcept = 0;
    /// @brief Возвращает стабильное имя реализации для отчётов и benchmark-ов.
    virtual std::string_view name() const noexcept = 0;
    /// @brief Возвращает минимальный набор evaluation keys для стадии.
    /// @return Требования, которые должны быть удовлетворены до execute().
    virtual BootstrapKeyRequirements keyRequirements() const = 0;
    /**
     * @brief Применяет стадию к ciphertext.
     *
     * Реализация сама отвечает за допустимый входной уровень, scale и
     * математические преобразования. Pipeline не добавляет relinearize,
     * rescale или modulus switch неявно.
     * @throws std::runtime_error Если вход не удовлетворяет контракту стадии.
     */
    virtual BootstrapStageResult execute(SealAdapter& adapter, const Cipher& input) const = 0;
};

/**
 * @brief Итог запуска pipeline вместе с отчётами по всем стадиям.
 *
 * Порядок элементов stages совпадает с порядком, в котором стадии были
 * переданы конструктору BootstrapPipeline.
 */
struct BootstrapPipelineResult {
    /// Итоговый ciphertext после всех стадий.
    Cipher ciphertext;
    /// Отчёты в порядке выполнения стадий.
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
    /**
     * @brief Создаёт неизменяемую последовательность bootstrap-стадий.
     * @throws std::invalid_argument Если список пуст или содержит null-стадию.
     */
    explicit BootstrapPipeline(std::vector<std::unique_ptr<BootstrapStage>> stages);

    /**
     * @brief Объединяет требования к ключам всех стадий.
     *
     * Дублирующиеся rotation steps удаляются, а результат сортируется по
     * возрастанию для воспроизводимой генерации ключей.
     * @throws std::invalid_argument Если стадия запросила rotation step 0.
     */
    BootstrapKeyRequirements keyRequirements() const;
    /**
     * @brief Выполняет стадии в заданном порядке и формирует фактический отчёт.
     *
     * До вызова все requirements должны быть удовлетворены: pipeline умеет
     * их агрегировать, но не генерирует и не загружает ключи сам.
     * @return Итоговый ciphertext и один BootstrapStageReport на каждую стадию.
     */
    BootstrapPipelineResult run(SealAdapter& adapter,
                                const Cipher& input,
                                const BootstrapMeasurementConfig& measurements = {}) const;

private:
    std::vector<std::unique_ptr<BootstrapStage>> stages_;
};

} // namespace m2424
