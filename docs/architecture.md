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
- `encode`, `encodeFor`, `encodeScalarFor`, `decode`, `encrypt`, `decrypt`;
- `add`, `sub`, `addPlain`, `subPlain`;
- `multiplyPlain`, `multiplyPlainAndRescale`, `multiplyRelinearizeAndRescale`;
- `rotate`;
- `modSwitchTo`;
- сериализация ключей и ciphertext;
- получение параметров ciphertext: `scale`, `chain_index`, `coeff_modulus_size`, размер объекта.

Сериализация должна поддерживать разделение контекстов: один контекст шифрует данные, второй выполняет вычисления с evaluation keys, третий расшифровывает результат.

Имена публичного API используют `camelCase`. Суффикс `For` означает, что plaintext кодируется на уровне и с масштабом указанного ciphertext; это нужно для совместимости plaintext-операций CKKS. `SealAdapter` не содержит ручной ModUp или произвольную подмену scale. Единственное исключение — явно названный `alignForAddition`: он допускает только близкие scale (до 1%) после снижения уровня и предназначен исключительно для сложения.

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

Поддерживаются только исследовательские модели и planner-ы:

- `BootstrapPipelinePlan` и `BootstrapRefreshPlanningResult` — расчёт требований к refresh без запуска ciphertext-refresh;
- `BootstrapScaleDesign` и layout-модели — проверка масштаба, периода, уровней и ёмкости EvalMod;
- `BootstrapStcReferencePlan` — plaintext/reference модель pre-EvalMod lattice form.

Backend-и линейных transform-стадий:

- `DenseDiagonal` — reference-путь для малых `slots`;
- `FftLike` — staged FFT-like путь для масштабируемых трансформаций.

## Контракты стадий

- `ModRaise` меняет RNS-базу ciphertext и не является восстановлением вычислительной глубины сам по себе.
- `CoeffToSlot` и `SlotToCoeff` должны сохранять значение в пределах заданного error budget.
- `EvalMod` применяется только к нормализованному входу, попадающему в допустимый интервал аппроксимации.
- `output_scale_repair` возвращает результат к рабочему CKKS-масштабу после EvalMod.

## Границы применимости

- Полный production-level bootstrapping не входит в стабильную поверхность `m2424/m2424.hpp`.
- Scalable refresh не должен рассматриваться как полный цикл с гарантией неограниченного продолжения вычислений при `1e-9`.
- Корректный bootstrap требует совпадения encrypted `ScaleDown/ModUp -> CoeffToSlot` с reference lattice form `P * (k + gamma * m)` до запуска EvalMod.
- Увеличение степени EvalMod не исправляет mismatch pre-EvalMod lattice form.
- Для внешнего применения нужны отдельные протоколы управления ключами, форматы обмена данными и анализ side-channel рисков.

## Модули реализации

- `src/core/` — базовые типы, профили, метрики точности и ABFT.
- `src/ckks/` — адаптер Microsoft SEAL, проверяемый evaluator и бюджет операций.
- `src/math/` — полиномы и rotation-based линейные преобразования.
- `src/planning/` — подбор параметров и отчёты по профилю и безопасности.
- `src/research/` — исследовательские модели и планы bootstrapping, DFT, EvalMod и Mod1.
