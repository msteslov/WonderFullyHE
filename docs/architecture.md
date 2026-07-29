# Архитектура WonderFullyHE

WonderFullyHE определяет C++17 API для приближённых защищённых вычислений на CKKS поверх Microsoft SEAL.

```text
Приложения и тесты
        |
        v
Публичный API m2424
        |
        +-- SealAdapter
        +-- Bootstrap candidates
        +-- Canonical embedding reference
        +-- CoeffToSlot contract
        +-- CoeffToSlot FFT plan
        +-- CoeffToSlot BSGS plan
        +-- Homomorphic linear transform
        +-- accuracy
        +-- abft
        +-- ProfileReport и SecurityReport
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
- `multiplyPlain`, `multiply`, `relinearize`, `rescaleToNext`;
- `rotate`;
- `modSwitchTo`;
- сериализация ключей и ciphertext;
- получение параметров ciphertext: `scale`, `chain_index`, `coeff_modulus_size`, размер объекта.

Сериализация должна поддерживать разделение контекстов: один контекст шифрует данные, второй выполняет вычисления с evaluation keys, третий расшифровывает результат.

Имена публичного API используют `camelCase`. Суффикс `For` означает, что plaintext кодируется на уровне и с масштабом указанного ciphertext; это нужно для совместимости plaintext-операций CKKS. Адаптер не меняет scale ciphertext вручную.

### `accuracy`

Модуль `accuracy` задаёт единый способ сравнения результата с эталоном:

```text
max_abs_error
mean_abs_error
compare(expected, actual, tolerance)
```

Эти метрики являются общим критерием точности для тестов, демонстрационных программ и benchmark-ов.

### `BootstrapCandidate`

`BootstrapCandidate` — это каталог целевых конфигураций для будущих экспериментов,
а не реализация bootstrap. Каждый кандидат содержит число слотов, степень кольца,
целевой scale, диапазон входа, лимиты ресурсов и консервативное распределение
общей ошибки `1e-9`. `forecastBootstrapFeasibility` возвращает положительный
результат только при наличии подтверждённой верхней границы для каждого
компонента и при соблюдении всех частных бюджетов.

### `CanonicalEmbeddingReference`

`CanonicalEmbeddingReference` задаёт plaintext-эталон CoeffToSlot и SlotToCoeff
в порядке слотов Microsoft SEAL. Он намеренно имеет квадратичную сложность и
используется только для проверки будущих encrypted-стратегий, а не в runtime.

### `CoeffToSlotContract`

`CoeffToSlotContract` фиксирует семантику входного ciphertext после будущего
ModUp: нормализованные коэффициенты, input/output scale и частный бюджет ошибки.
Конкретный FFT- или BSGS-план передаёт минимальный запас уровней и необходимые
rotation keys; `preflightCoeffToSlot` отклоняет несовместимый ciphertext до
запуска вычисления.

### `CoeffToSlotFftPlan`

`CoeffToSlotFftPlan` — первый encrypted-кандидат CoeffToSlot для размеров 4, 8
и 16. План включает bit-reversal и все butterfly-слои, поэтому его стоимость
ротаций и уровней измеряется полностью, без бесплатных перестановок. Пока нет
ModUp, encrypted-тест использует synthetic slot fixture; это проверка исполнения
линейного FFT-плана, а не подтверждение полного bootstrap CoeffToSlot.

### `CoeffToSlotBsgsPlan`

`CoeffToSlotBsgsPlan` — второй encrypted-кандидат для той же DFT-матрицы.
Он переиспользует baby rotations внутри giant-групп и расходует один уровень;
результаты сравниваются с FFT только на одинаковом synthetic fixture. Пока ModUp
не реализован, это не подтверждение полного bootstrap CoeffToSlot.

Условия объективного сравнения обоих кандидатов зафиксированы в
[`coeff_to_slot_comparison_plan.md`](coeff_to_slot_comparison_plan.md).
Воспроизводимая матрица запусков и её результаты — в
[`coeff_to_slot_matrix.md`](coeff_to_slot_matrix.md).

### `HomomorphicLinearTransform`

`HomomorphicLinearTransform` исполняет проверяемую сумму plaintext-диагоналей
и rotations с одним контролируемым rescale. Это вычислительное ядро будущих FFT и BSGS,
но само по себе не является готовым CoeffToSlot.

### `abft`

ABFT-слой реализует checksum-проверки для:

- сложения;
- вычитания;
- покомпонентного умножения;
- ротации.

Эти проверки контролируют численную согласованность вычисления и не заменяют криптографическую аутентификацию результата.

## CKKS-профили

`m2424::profiles` содержит профили:

- `fast_demo_ckks`
- `basic_ckks`
- `balanced_ckks`
- `depth_ckks`
- `high_precision_ckks`

Таблица профилей должна формироваться из исполняемого отчёта:

```bash
./build/demo_profile_report
```

## Модули реализации

- `src/core/` — базовые типы, профили, метрики точности и ABFT.
- `src/ckks/` — адаптер Microsoft SEAL.
- `src/planning/` — отчёты по профилю и безопасности.

## Артефакты сборки

Корень `build/` используется только CMake и CTest. Исполняемые файлы располагаются в `build/bin/`, статические и динамические библиотеки — в `build/lib/`, а HTML-документация Doxygen — в `build/docs/html/`.

## Контракт полноценного CKKS bootstrapping

Bootstrap считается реализованным только как единый рабочий путь, а не как набор
разрозненных математических примеров:

```text
входной CKKS ciphertext на q
        → ModRaise(q → Q)
        → CoeffToSlot
        → EvalMod
        → SlotToCoeff
        → ModDown и нормализация
        → обновлённый CKKS ciphertext
```

### Обязательные блоки

1. **ModRaise.** Низкоуровневый RNS basis extension обоих компонентов ciphertext:
   перевод из NTT в coefficient form, centered lift из q в Q, обратный NTT и
   назначение верхнего `parms_id`. Операция не заменяется `modSwitchTo` и не
   должна требовать secret key.
2. **CoeffToSlot.** Линейное преобразование именно ciphertext после ModRaise,
   а не обычного CKKS slot fixture. Входной порядок коэффициентов, scale,
   уровни и rotation keys фиксируются контрактом ModRaise.
3. **EvalMod.** Одна подтверждённая аппроксимация modular reduction с явной
   степенью, схемой вычисления, relinearization/rescale-последовательностью и
   доказанной верхней границей ошибки.
4. **SlotToCoeff.** Обратное к CoeffToSlot преобразование с тем же порядком
   слотов, нормировкой и контролируемым расходом уровней.
5. **ModDown.** Финальное снижение modulus и нормализация scale до контракта
   выходного ciphertext.

### Общие требования к реализации

- один публичный `bootstrap` API и один тип конфигурации; пользователь не
  вручную не соединяет внутренние блоки;
- каждый блок имеет строгий входной и выходной контракт: representation,
  `chainIndex`, scale, необходимый набор ключей, расход уровней и максимум
  ошибки;
- selector выбирает только стратегии, подтверждённые для данного профиля; для
  непроверенной конфигурации возвращается ошибка, а не fallback;
- plaintext-константы и ключи готовятся вне рабочего вызова, а `bootstrap`
  не выполняет скрытых encode/allocations, не необходимых алгоритму;
- все имена публичного C++ API используют camelCase, а Doxygen-комментарии
  написаны по-русски.

### Критерии готовности

Реализация становится production-ready только при одновременном выполнении:

- весь маршрут выполняется без secret key;
- результат расшифровывается исходным ключом и соответствует входному сообщению;
- итоговая максимальная абсолютная ошибка не превышает `1e-9`, а каждый блок
  соблюдает свой частный бюджет;
- scale, уровень и размер ciphertext после bootstrap соответствуют контракту
  свежего ciphertext;
- тесты включают граничные и случайные комплексные входы, несколько независимых
  шифрований, отрицательные сценарии ключей/уровней/scale и полный цикл;
- benchmark фиксирует latency, память, размеры evaluation keys и уровень
  безопасности для каждой поддержанной конфигурации;
- нет synthetic fixture, непроверенных алгоритмов или экспериментальных
  fallback-ветвей в публичном рабочем пути.
