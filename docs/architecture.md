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

## Структура реализации

Стабильный код и research-контуры разделены на уровне сборки и файлов:

- `src/bootstrap.cpp` — публичный `Bootstrapper` и guarded entry points.
- `src/bootstrap_plan.cpp` — planner/gate layer: уровни, scale, period, Mod1 depth и security budget.
- `src/bootstrap_scaling.cpp` — scalar decomposition и scale-management primitives.
- `src/bootstrap_dft.cpp` — FFT-like transform plans для bootstrap linear transforms.
- `src/bootstrap_prototype.cpp` — общий experimental prototype для dense/reference path.
- `src/bootstrap_prototype_stc_first.cpp` — experimental `SlotsToCoeffsFirst` circuit.
- `src/bootstrap_prototype_detail.cpp` — internal report/stage/scale helpers, не публичный API.

Bootstrap research apps не собираются по умолчанию. Они включаются только через
`-DM2424_BUILD_RESEARCH_APPS=ON`, чтобы default build показывал стабильную
поверхность библиотеки, а не набор диагностических harnesses.

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

Full refresh не считается стабильным API, пока `bench_bootstrap_validation` не проходит end-to-end. Старые one-off harnesses для period/scaling/prescale/layout расследований удалены: их вывод перенесён в planner objects, CTest и архитектурные ограничения. `BootstrapScaleDesign` превращает scaling decision в явный объект со статусами `period_model_blocked`, `scale_strategy_blocked`, `evalmod_capacity_blocked` и `ready_for_evalmod_p3`, а `BootstrapPeriodFeasibilityWindow` проверяет, существует ли период, который одновременно нормализует magnitude и оставляет EvalMod scale capacity. `bench_bootstrap_reference_path` оставлен как E0 oracle-harness: для маленьких slots он держит CPU/reference path рядом с ciphertext path и печатает domain, scale, level, expected magnitude, error и blocker на каждом этапе.

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

`plan_bootstrap_layout` считает circuit-level feasibility до ciphertext-прогона: period, normalization levels, scale-squash levels, EvalMod levels, transform levels, residual levels, total modulus bits и security limit. Он возвращает физический `CkksProfile`, где `profile.slots` остаётся полным CKKS slot count (`N/2`), а не числом логических bootstrap slots. Это нужно для корректного кодирования diagonal masks в rotation-based transforms. Первый расчёт показывает: dense-like layout требует около `600` bits при лимите `438` для `16384/tc128`; FFT-like `CoeffToSlot` с текущим inverse reference требует около `760` bits и помещается в `32768/tc128`; полностью factorized inverse с таким же числом слоёв требует около `1040` bits и не помещается в `32768/tc128`.

Физический gate для `ModRaiseFirst + FFT-like CoeffToSlot` показал текущий blocker: после `ModRaise` диагностическая амплитуда доминирует, и один только prescale внутри transform-а не даёт готовый вход `EvalMod` под `32768/tc128`. Поэтому ближайший рабочий путь для refresh — `SlotsToCoeffsFirst`, где сначала выполняется `SlotToCoeff`, затем `ModRaise`, затем обратный `CoeffToSlot` и нормализация перед `EvalMod P3`.

Этот путь вынесен в публичный API как `Bootstrapper::refresh_slots_to_coeffs_first(...)`. Для него нужно генерировать Galois-ключи через `Bootstrapper::scalable_refresh_rotation_steps(slots)`. На текущем этапе он является guarded refresh path: значение сохраняется в пределах tolerance и `EvalMod P3` проходит физически, но полноценное восстановление цепочки ещё требует отдельного post-refresh/modulus management слоя.

`BootstrapPrototypeReport::continuation_levels` фиксирует, сколько уровней остаётся после refresh-круга. Это текущая основная метрика прогресса к полноценному bootstrapping: нужно увеличить её и затем сделать `restore_level_criterion=true` без ухудшения ошибки. Сырой `post_refresh_mod_raise` оставлен диагностическим флагом; он показывает, что один structural raise недостаточен, потому что может поднять chain и одновременно нарушить tolerance.

После `EvalMod P3` выполняется `output_scale_repair`: plaintext-rescale с коэффициентом `1` и рассчитанным plaintext scale возвращает результат к рабочему scale около `2^40`. Это делает refresh-круг повторяемым: `bench_bootstrap_multi_cycle` уже проходит несколько циклов подряд, но остаток уровней остаётся главным ограничением.

FFT-like `CoeffToSlot` применяет попарное composition соседних layers. Для `slots=16` это снижает число ciphertext-level transform stages с `8` до `4`; multi-cycle baseline после этого проходит `3` refresh-круга подряд с минимумом `8` continuation levels. Проверочный 57-битный layout показывает, что один рост modulus bits не решает точность: ошибка остаётся порядка `1e-5`, поэтому следующий точностной blocker находится в modular-reduction/rounding semantics.

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
- `bench_bootstrap_reference_path` является oracle-harness для маленьких slots и stage-by-stage reference comparison.
- `test_bootstrap_scalable_refresh` проверяет публичный `SlotsToCoeffsFirst/FftLike/P3` путь через `Bootstrapper`, а не через отдельный benchmark.
- `demo_bootstrap_cipher_path` запускает experimental refresh-путь от существующего ciphertext.
- `bench_ckks` измеряет время операций, численную ошибку и размеры сериализованных объектов.
- `bench_chain_accuracy` калибрует связь между chain length, `scale_log2`, рабочей битностью и точностью.
- `bench_parameter_planner` проверяет выбранные planner-ом профили на реальном SEAL-прогоне.
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
