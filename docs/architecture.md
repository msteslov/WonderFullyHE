# Архитектура

WonderFullyHE реализует библиотеку защищённых вычислений поверх Microsoft SEAL.

```text
Демонстрационные приложения и эксперименты
        |
        v
Публичный API WonderFullyHE
        |
        +-- SealAdapter: CKKS-контекст, ключи, encode/encrypt/evaluate/decrypt
        +-- accuracy: отчёты по максимальной и средней ошибке
        +-- abft: checksum-проверки корректности вычислений
        +-- LinearTransform: линейные преобразования через ротации
        +-- PolynomialEvaluator: вычисление полиномов от ciphertext
        +-- Bootstrapper: диагностика глубины и отчёт по bootstrapping-конвейеру
        +-- profiles: готовые CKKS-профили для типовых сценариев
        +-- parameter planning: целевой слой подбора scale/chain/N по требованиям
        +-- profile_report: воспроизводимые таблицы CKKS-параметров
        |
        v
Microsoft SEAL как CKKS-движок
```

## Слои публичного API

### SealAdapter

`SealAdapter` отделяет публичный API библиотеки от Microsoft SEAL. Прямые SEAL-типы скрыты за `Plain` и `Cipher`, а наружу вынесены CKKS-операции, используемые библиотекой:

- создание CKKS-контекста из `CkksProfile`;
- генерация ключей;
- encode/decode;
- encrypt/decrypt;
- сложение и вычитание;
- сложение и вычитание ciphertext с plaintext;
- умножение ciphertext на plaintext с rescale;
- умножение с relinearize и rescale;
- ротация вектора;
- генерация только нужных Galois-ключей для заданных ротаций;
- выравнивание уровня и масштаба ciphertext;
- сериализация и загрузка ключей/ciphertext;
- размеры сериализованных объектов;
- диагностика ciphertext: `scale`, `chain_index`, `coeff_modulus_size`.

Сериализация возвращает байтовые буферы для публичного ключа, секретного ключа, Relin/Galois-ключей и ciphertext. Это позволяет создавать отдельные контексты выполнения: один контекст шифрует данные, другой выполняет гомоморфные операции только с evaluation keys, третий расшифровывает результат.

Microsoft SEAL остаётся внутренней зависимостью. Текущая C++-реализация размещает публичные типы и функции в пространстве имён `m2424`.

### accuracy

Модуль `accuracy` задаёт общие критерии численной точности:

```text
max_abs_error
mean_abs_error
compare(expected, actual, tolerance)
```

Демонстрационные приложения и тесты используют этот модуль как единый критерий точности.

### operation_budget

`OperationBudgetBuilder` собирает `CkksOperationBudget` для planner-а. `CheckedEvaluator` уже накапливает бюджет для своих операций автоматически, поэтому типичный проверяемый pipeline может быть использован как вход в `plan_ckks_parameters`. Прямые вызовы `SealAdapter` остаются низкоуровневым API; для них budget нужно записывать явно.

### abft

Модуль `abft` реализует checksum-проверки корректности. Текущие проверки покрывают:

- добавленный checksum для сложения;
- добавленный checksum для вычитания;
- эталонный checksum для покомпонентного умножения;
- сохранение суммы при ротации.

ABFT-слой используется для проверки численной согласованности защищённых вычислений.

### LinearTransform

`LinearTransform` реализует преобразования вида:

```text
sum_i a_i * rotate(ct, k_i)
```

Этот блок нужен для этапов `CoeffToSlot` и `SlotToCoeff`, где вычисление строится из ротаций CKKS-слотов, умножений на plaintext-коэффициенты и сложений.

Отдельная функция `sum_slots` считает сумму заданного числа слотов и помещает результат в первый слот ciphertext. Для неё используется последовательность ротаций по степеням двойки.

### PolynomialEvaluator

`PolynomialEvaluator` вычисляет полином от ciphertext по заданным степеням и коэффициентам. Блок нужен для этапа `EvalMod`, где модульная редукция приближается полиномом.

### Bootstrapper

`Bootstrapper` является experimental точкой входа для bootstrapping-модуля. Текущая реализация предоставляет:

- диагностику мультипликативной глубины;
- отчёт по CKKS bootstrapping-конвейеру;
- явный `BootstrapPipelinePlan` с доменами значений, backend'ом трансформаций и active gate;
- явные этапы: `ModRaise`, `CoeffToSlot`, `eval_mod_normalization`, `EvalMod`, `SlotToCoeff`, `post_refresh_mod_raise`;
- параметры ciphertext на границе глубины: `scale`, `chain_index`, `coeff_modulus_size`, `serialized_bytes`;
- structural/scaling diagnostics для проверки `ModRaise -> CoeffToSlot` перед full refresh.

Full refresh не считается стабильным API, пока `bench_bootstrap_validation` не проходит end-to-end. `bench_bootstrap_scaling` теперь является первым gate: он поддерживает source-period diagnostics и decomposed plaintext scaling для tiny normalization scalar. Если нет real period-mode с EvalMod-ready magnitude and scale, `bench_bootstrap_validation` не запускает `EvalMod`. `bench_bootstrap_scale_strategy` сначала решает тот же scale/levels вопрос алгебраически по modulus-drop budget и capacity первого EvalMod multiplication, `BootstrapScaleDesign` превращает этот вывод в явный decision object со статусами `period_model_blocked`, `scale_strategy_blocked`, `evalmod_capacity_blocked` и `ready_for_evalmod_p3`, а `BootstrapPeriodFeasibilityWindow` проверяет, существует ли вообще период, который одновременно нормализует magnitude и оставляет P3 scale capacity. Для текущего dense/research path окно отрицательное: `boot_ckks` около `-52.5` бит, `boot_deep_ckks` около `-56.7` бит. `bench_bootstrap_prescaled_coeff_to_slot` проверяет попытку перенести часть normalization в `CoeffToSlot` diagonals; через `apply_at_plain_scale` prescaled diagonals становятся физически представимыми, но высокий plaintext scale поднимает scale результата и пока оставляет blocker `not_enough_levels_for_scale`. `bench_bootstrap_profile_budget` проверяет гипотетические modulus-chain layouts без ciphertext и выводит минимальный `period_offset_log2`, требуемый для P3-capacity при заданном target scale; для `target_scale_log2=60` текущий offset `44` невозможен, а нижняя граница равна `102`. `bench_bootstrap_period_model` проверяет period/scaling без EvalMod на ciphertext, а `bench_bootstrap_one_case` используется для stage trace вокруг normalization/EvalMod/denormalization. `bench_bootstrap_reference_path` добавляет E0 oracle-harness: для `slots=4/8/16` он держит CPU/reference path рядом с ciphertext path и печатает domain, scale, level, expected magnitude, error и blocker на каждом этапе. Diagnostic `boot_deep_ckks` работает на максимальном для `32768/tc128` modulus budget и нужен только для локализации blocker: если даже он не даёт P3-ready scale после value-preserving `ScaleSquash`, следующий слой должен менять scale-management, а не EvalMod coefficients или transforms.

Prescaled `CoeffToSlot` теперь проверяется через `apply_at_plain_scale`: это снимает transparent plaintext failure, но показывает следующий tradeoff - высокий plaintext scale делает transform физически представимым и одновременно поднимает ciphertext scale, из-за чего scale gate всё ещё блокируется по modulus-drop/capacity.
Физический diagnostic `CoeffToSlot -> normalization -> ScaleSquash` сейчас даёт `p3_ready_cases=0` на `boot_ckks`: ближайшие cases проходят scalar/interval, но не успевают снизить scale перед `EvalMod` при сохранении нужных levels.

Архитектурный target разделён на два backend'а:

- `DenseDiagonal`: текущий research backend для малых `slots`, где `CoeffToSlot` и `SlotToCoeff` строятся через полную canonical embedding matrix и diagonal decomposition.
- `FftLike`: целевой scalable backend для больших `slots`, где transforms должны быть staged FFT-like линейными преобразованиями вместо плотных таблиц. Первая часть реализована для `CoeffToSlot`: canonical embedding раскладывается через рекурсивный even/odd split, а dense matrix остаётся reference для проверки.

Контракты стадий фиксируют, что `ModRaise` является structural-only, `Scaling` является первым correctness gate, `EvalMod` проверяется отдельно от transform backend, а `SlotToCoeff` возвращается в correctness path только после прохождения предыдущих gates.

Целевой `BootstrapPlanner` должен стать обязательным входом в refresh. Он заранее считает для каждой стадии `chain_index`, `scale_log2`, остаток modulus, bound значения и error budget. Если план не проходит, refresh не запускается и возвращает явную причину: `blocked_by_levels`, `blocked_by_period_window`, `blocked_by_evalmod_capacity` или `blocked_by_error_budget`.

Первая часть этого слоя реализована как `plan_bootstrap_refresh`: она принимает текущий `CipherInfo` и `CkksOperationBudget` следующего вычислительного блока, проверяет хватает ли chain levels без refresh и возвращает `compute_fits_without_refresh`, `refresh_required` или `refresh_plan_blocked`.

`Bootstrapper::plan_refresh_for_budget` применяет этот gate к реальному ciphertext и должен стать входной точкой перед любым будущим guarded refresh.

`plan_bootstrap_refresh_scale_gate` является вторым gate для prototype-refresh: принимает факты после `ModRaise -> CoeffToSlot`, считает period/scaling design и блокирует запуск, если normalization scalar, остаток levels или capacity первого `EvalMod` multiplication не проходят.

`search_bootstrap_refresh_scale_gate` добавляет controlled sweep по period/plain-scale/target-scale и возвращает лучший design. Его задача - показать, существует ли feasible scale-management case при текущем profile/path, и если нет, какой blocker ближайший. `BootstrapScaleStrategyPlan` также хранит `missing_drop_log2`, `missing_scalar_levels` и `missing_total_levels`.

Search учитывает prescale вместе с plaintext scale transform-а. Это важно для dense `CoeffToSlot`: prescale может сделать magnitude algebraically ready, но если transform требует высокий plaintext scale, первый `EvalMod` multiplication снова блокируется по scale capacity.

Для `EvalMod` default-полиномом остаётся `P3`, пока нормализованный вход удовлетворяет `|u| <= 2^-10`. На этом интервале математическая ошибка аппроксимации около `1e-14`, поэтому для target `1e-9` bottleneck находится в CKKS scale/noise и линейных трансформах, а не в степени полинома.

### Parameter planning

Ручные профили не должны быть основным способом выбора параметров. Модуль `parameter_planner` строит профиль по требованиям:

```text
target_error
multiplicative_depth
slots
ops_profile
security_bits
```

Для нетривиальных программ planner должен получать явный `CkksOperationBudget`: число `add/sub`, plaintext ops, ciphertext multiplications, plaintext-rescale шагов, explicit rescale, mod-switch, rotations, linear transforms, EvalMod P3 и refresh-операций. Coarse `ops_profile` остаётся fallback для простых сценариев и быстрых demo.

Расчётная схема:

```text
required_result_bits = ceil(-log2(target_error))
work_bits = required_result_bits + calibrated_loss_bits(ops_profile)
scale_log2 ~= work_bits
work_levels >= multiplicative_depth
poly_modulus_degree = минимальный N, проходящий security
```

Для `target_error = 1e-9` текущая калибровка depth=2 даёт `work_bits ~= 44`, поэтому быстрый режим выбирает 45-битные рабочие модули, а conservative-режим — 50-битные. Первая реализация возвращает `CkksPlanningResult` с `CkksProfile`, выбранными битами, `estimated_abs_error_bound`, `passes_target_error` и `SecurityReport`.

### profiles

Модуль `profiles` содержит готовые наборы `CkksProfile` для типовых режимов:

- `fast_demo_ckks` — быстрые локальные проверки;
- `basic_ckks` — основные демонстрации и benchmark базовых операций;
- `balanced_ckks` — дополнительная глубина на `N = 8192`;
- `depth_ckks` — анализ глубины и bootstrapping-блоки;
- `high_precision_ckks` — профиль с повышенным масштабом для сценариев, где важнее точность.

Это снижает риск ручной ошибки при выборе `poly_modulus_degree`, `coeff_modulus_bits` и `scale`.

### profile_report

`profile_report` формирует воспроизводимый отчёт по CKKS-параметрам:

```bash
./build/demo_profile_report
```

Вывод содержит `N`, число слотов, цепочку модулей, суммарный размер коэффициентного модуля, масштаб и оценку мультипликативной глубины.

## Демонстрационные приложения

- Default CTest содержит только инвариантные tests: `test_smoke`, `test_accuracy`, `test_adapter_failures`. Demos и benchmarks не входят в default CTest и запускаются явно.
- `demo_secure_stats` показывает защищённую агрегацию: сумму и среднее над зашифрованными данными.
- `demo_client_compute_roundtrip` показывает передачу сериализованных ключей/ciphertext между отдельными контекстами выполнения.
- `demo_galois_key_optimization` сравнивает полный и ограниченный набор Galois-ключей для ротаций.
- `demo_abft` проверяет корректность гомоморфных операций через ABFT-инварианты.
- `demo_noise_growth` показывает расход глубины и остановку вычисления без bootstrapping.
- `demo_bootstrap_pipeline` печатает отчёт bootstrapping-модуля.
- `bench_bootstrap_scaling` проверяет представимость normalization scalar после `ModRaise -> CoeffToSlot`.
- `bench_bootstrap_reference_path` является oracle-harness для маленьких slots и stage-by-stage reference comparison.
- `demo_bootstrap_cipher_path` запускает experimental refresh-путь от существующего ciphertext.
- `demo_bootstrap_end_to_end` является historical experimental demo и не входит в default CTest.
- `bench_ckks` измеряет время операций, численную ошибку и размеры сериализованных объектов.
- `bench_chain_accuracy` калибрует связь между chain length, `scale_log2`, рабочей битностью и точностью.
- `bench_parameter_planner` проверяет выбранные planner-ом профили на реальном SEAL-прогоне.
- `bench_bootstrap_parts` измеряет `mul_plain_rescale`, `linear_transform`, `sum_slots` и `polynomial_eval`.
- `bench_bootstrap_refresh` сначала проверяет `Bootstrapper::plan_refresh_for_budget`, затем локальный scale gate для текущего prototype-refresh. Experimental refresh запускается только если следующий блок требует refresh и scale gate готов к `EvalMod`; иначе benchmark печатает blocker и завершается без падения.
- `bench_parallel_throughput` измеряет throughput при параллельной обработке независимых ciphertext.
- `demo_profile_report` печатает таблицу CKKS-параметров.

## Следующие шаги реализации

Следующий этап реализации:

1. расширить calibration layer `ops_profile -> calibrated_loss_bits` реальными sweep-измерениями;
2. сделать `BootstrapPlanner` обязательным gate перед refresh;
3. реализовать factorized FFT-like `CoeffToSlot/SlotToCoeff`, оставив dense backend как reference;
4. оптимизировать rotation keys, BSGS/hoisting и кеширование диагоналей;
5. сравнить строительные блоки с OpenFHE/Lattigo на одинаковых параметрах.
