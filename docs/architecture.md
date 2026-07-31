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
        +-- CoeffToSlotPlan
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

`CoeffToSlotContract` фиксирует семантику входного `RaisedCipher` после
ModRaise: нормализованные коэффициенты, input/output scale и частный бюджет
ошибки. План передаёт минимальный запас уровней, rotation keys и требование
conjugation key; `preflightCoeffToSlot` отклоняет несовместимый RaisedCipher до
запуска вычисления.

### `CoeffToSlot`

`CoeffToSlot` принимает только `RaisedCipher` и вычисляет две половины
coefficient-side plaintext через factorized canonical inverse embedding
`U₀ᵀ/N`, rotations, plaintext multiplication и conjugated sums. Результат содержит
`slotCipherFirst` и `slotCipherSecond`. Независимый oracle расшифровывает
RaisedCipher, выполняет inverse NTT и centered CRT по поднятой базе Q, после
чего отдельно сравнивает обе половины. CRT по исходной базе q используется
только для проверки ModRaise residues и будущего результата EvalMod.

`CoeffToSlotPlan` является единственным рабочим планом этого преобразования.
Он строит radix-2 special FFT в порядке корней Microsoft SEAL, добавляет
разреженные стадии перестановки bit-reversal и группирует соседние факторы
согласно `factorization().radices`. План заранее хранит только диагонали этих
разреженных факторов, вычисляет точный набор rotations и реальную глубину.
Исходные radix-2 факторы освобождаются сразу после merge. Внутри каждого
сгруппированного фактора rotations исполняются по BSGS. Baby-step rotations
одного ciphertext используют single hoisting: DCRT-разложение второй
компоненты выполняется один раз на фактор, после чего переиспользуется для
всей пачки automorphisms. Giant-step rotations остаются обычными; их
совместное исполнение относится к отдельной double-hoisting оптимизации.
Все plaintext multiplication по-прежнему суммируются до единственного rescale.

Encoded plaintext для каждого уровня обязательно готовятся явным `prepare`
после preflight и до рабочего `apply`. `prepare` возвращает отдельный
неизменяемый `PreparedCoeffToSlotPlan`, привязанный к fingerprint
SEAL-контекста, точному входному `parms_id`, начальному уровню, scale и
контракту преобразования. `apply` повторяет preflight и отклоняет
несовместимый prepared-план; скрытого encode или mutable cache в рабочем
пути нет. Плотная матрица из 8192 диагоналей в runtime отсутствует.

`bench_coeff_to_slot` отдельно измеряет создание контекста, построение плана,
генерацию ключей, ModRaise, oracle, `prepare`, первый и повторный `apply`, а
также публикует число операций, объём prepared plaintext и лучшие варианты
radix-разбиения по cost model.

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
- `src/math/` — canonical reference и factorized `CoeffToSlotPlan`.
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
