# Архитектура WonderFullyHE

WonderFullyHE определяет C++17 API для приближённых защищённых вычислений на CKKS поверх Microsoft SEAL. Публичная поверхность библиотеки разделена на стабильный слой `m2424/m2424.hpp` и исследовательский слой `m2424/experimental.hpp`.

```text
Приложения и тесты
        |
        v
Публичный API m2424
        |
        +-- SealAdapter
        +-- CheckedEvaluator и OperationBudget
        +-- accuracy
        +-- abft
        +-- LinearTransform
        +-- PolynomialEvaluator
        +-- ParameterPlanner
        +-- ProfileReport и SecurityReport
        +-- Experimental bootstrap API
        |
        v
Microsoft SEAL
```

## Слои API

### `SealAdapter`

`SealAdapter` изолирует публичный API WonderFullyHE от типов Microsoft SEAL и предоставляет операции:

- создание CKKS-контекста из `CkksProfile`;
- генерация public/secret/relinearization/Galois-ключей;
- `encode`, `decode`, `encrypt`, `decrypt`;
- `add`, `sub`, `add_plain`, `sub_plain`;
- `mul_plain`, `mul_plain_rescale`, `mul_relin_rescale`;
- `rotate`;
- `mod_switch_to` и `match_level_and_scale`;
- сериализация ключей и ciphertext;
- получение параметров ciphertext: `scale`, `chain_index`, `coeff_modulus_size`, размер объекта.

Сериализация должна поддерживать разделение контекстов: один контекст шифрует данные, второй выполняет вычисления с evaluation keys, третий расшифровывает результат.

### `CheckedEvaluator` и `OperationBudget`

`CheckedEvaluator` выполняет операции через `SealAdapter`, собирает `CkksOperationBudget` и возвращает `CheckedResult` с численными метриками. Этот бюджет используется planner-ом для выбора параметров и проверки необходимости refresh перед следующим вычислительным блоком.

`OperationBudgetBuilder` нужен для явного описания вычисления, когда операции вызываются напрямую через `SealAdapter`.

### `accuracy`

Модуль `accuracy` задаёт единый способ сравнения результата с эталоном:

```text
max_abs_error
mean_abs_error
compare(expected, actual, tolerance)
```

Эти метрики являются общим критерием точности для тестов, демонстрационных программ и benchmark-ов.

### `abft`

ABFT-слой реализует checksum-проверки для:

- сложения;
- вычитания;
- покомпонентного умножения;
- ротации.

Эти проверки контролируют численную согласованность вычисления и не заменяют криптографическую аутентификацию результата.

### `LinearTransform`

`LinearTransform` вычисляет преобразования вида:

```text
sum_i a_i * rotate(ct, k_i)
```

Блок используется для rotation-based линейных преобразований и bootstrap-этапов `CoeffToSlot` / `SlotToCoeff`. Функция `sum_slots` сворачивает заданное число слотов в первый слот ciphertext через ротации по степеням двойки.

### `PolynomialEvaluator`

`PolynomialEvaluator` вычисляет полином от ciphertext по набору степеней и коэффициентов. В исследовательском bootstrap-контуре этот слой используется для EvalMod.

### `ParameterPlanner`

`plan_ckks_parameters` выбирает CKKS-профиль по требованиям:

```text
target_error
multiplicative_depth
slots
ops_profile или CkksOperationBudget
security_bits
```

Результат содержит выбранный `CkksProfile`, рабочую битность, масштаб, оценку ошибки и `SecurityReport`. Для нетривиальных программ должен использоваться `CkksOperationBudget`, потому что он учитывает разные типы операций: ciphertext multiplication, plaintext-rescale, rotations, linear transforms, EvalMod и refresh.

## CKKS-профили

`m2424::profiles` содержит профили:

- `fast_demo_ckks`
- `basic_ckks`
- `balanced_ckks`
- `depth_ckks`
- `high_precision_ckks`
- `boot_ckks`
- `boot_deep_ckks`
- `precision_boot_deep_ckks`
- `precision_boot_ultra_ckks_59`

Таблица профилей должна формироваться из исполняемого отчёта:

```bash
./build/demo_profile_report
```

## Bootstrap-слой

Исследовательский bootstrap API должен подключаться явно:

```cpp
#include "m2424/experimental.hpp"
```

Компоненты bootstrap-слоя:

- `Bootstrapper` — точка входа для анализа глубины и guarded refresh;
- `BootstrapPipelinePlan` — описание порядка стадий, backend-а трансформаций и активных проверок;
- `BootstrapRefreshPlanningResult` — решение, помещается ли следующий блок без refresh или нужен guarded refresh;
- `BootstrapPrototypeReport` — отчёт по стадиям prototype-refresh;
- `BootstrapScaleDesign` и layout-модели — проверка масштаба, периода, уровней и ёмкости EvalMod;
- `BootstrapStcReferencePlan` — plaintext/reference модель pre-EvalMod lattice form.

Backend-и линейных transform-стадий:

- `DenseDiagonal` — reference-путь для малых `slots`;
- `FftLike` — staged FFT-like путь для масштабируемых трансформаций.

## Исследовательский цикл

Изменяемые bootstrap-компоненты не должны выбираться в benchmark-ах через
локальные `bool` и набор setter-вызовов. Для этого research API предоставляет
`BootstrapExperimentConfig` и `run_bootstrap_experiment`:

```text
BootstrapExperimentConfig
        |
        +-- circuit order / transform backend / EvalMod
        +-- scaling, period и STC/ModUp параметры
        v
bootstrap_experiment_rotation_steps
        v
keygen только нужных Galois keys
        v
run_bootstrap_experiment(s)
        v
BootstrapExperimentResult { passed | blocked | failed, blocker, stage report }
```

Каждая конфигурация описывает один вариант, поэтому suite запускает одинаковый
input независимо для каждого варианта и возвращает сопоставимые отчёты. Runner
не создаёт ключи неявно: вызывающий код сначала объединяет rotation steps всех
конфигураций, создаёт adapter и ключи, затем запускает suite. Это исключает
скрытое влияние key-set и setup на сравнение backend-ов.

Пример исполняемого формата — research benchmark
`bench_bootstrap_experiments` (доступен с `M2424_BUILD_RESEARCH_APPS=ON`).

Scalable refresh-путь задаёт порядок:

```text
SlotsToCoeff -> ScaleDown/ModUp -> CoeffToSlot -> EvalMod P3 -> output_scale_repair
```

Набор Galois-ключей для этого пути должен строиться через:

```cpp
m2424::Bootstrapper::scalable_refresh_rotation_steps(slots)
```

## Контракты стадий

- `ModRaise` меняет RNS-базу ciphertext и не является восстановлением вычислительной глубины сам по себе.
- `CoeffToSlot` и `SlotToCoeff` должны сохранять значение в пределах заданного error budget.
- `EvalMod` применяется только к нормализованному входу, попадающему в допустимый интервал аппроксимации.
- `output_scale_repair` возвращает результат к рабочему CKKS-масштабу после EvalMod.
- Guarded refresh не запускает prototype-refresh, если следующий вычислительный блок помещается в доступные уровни.

## Границы применимости

- Полный production-level bootstrapping не входит в стабильную поверхность `m2424/m2424.hpp`.
- Scalable refresh не должен рассматриваться как полный цикл с гарантией неограниченного продолжения вычислений при `1e-9`.
- Корректный bootstrap требует совпадения encrypted `ScaleDown/ModUp -> CoeffToSlot` с reference lattice form `P * (k + gamma * m)` до запуска EvalMod.
- Увеличение степени EvalMod не исправляет mismatch pre-EvalMod lattice form.
- Для внешнего применения нужны отдельные протоколы управления ключами, форматы обмена данными и анализ side-channel рисков.

## Модули реализации

- `src/seal_adapter/seal_adapter.cpp` — адаптер Microsoft SEAL.
- `src/checked_evaluator.cpp` — проверяемый evaluator и сбор бюджета операций.
- `src/parameter_planner.cpp` — подбор CKKS-параметров.
- `src/profiles.cpp` — готовые CKKS-профили.
- `src/profile_report.cpp` и `src/security_report.cpp` — отчёты по параметрам.
- `src/linear_transform.cpp` и `src/diagonal_transform.cpp` — rotation-based линейные преобразования.
- `src/polynomial.cpp` и `src/eval_mod.cpp` — полиномы и EvalMod.
- `src/bootstrap.cpp` — публичный bootstrap facade.
- `src/bootstrap_plan.cpp` — планирование refresh и проверка уровней.
- `src/bootstrap_layout_v2.cpp` — layout-планирование bootstrap-профиля.
- `src/bootstrap_dft.cpp` — DFT-планы для `CoeffToSlot` / `SlotToCoeff`.
- `src/bootstrap_scaling.cpp` — управление масштабом bootstrap-стадий.
- `src/bootstrap_stc_modup.cpp` — STC-first ScaleDown/ModUp contract.
- `src/bootstrap_stc_reference.cpp` — plaintext/reference lattice model.
- `src/bootstrap_prototype*.cpp` — prototype-refresh реализации.
