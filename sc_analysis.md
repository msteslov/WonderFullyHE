# WonderFullyHE: сценарный анализ функционала

## Граница анализа

Отчёт описывает кодовую базу на текущем состоянии рабочей ветки: публичные
заголовки `include/m2424`, их реализации в `src` и реально зарегистрированные
тесты. Под «сценарием» ниже понимается путь от входных данных API до результата
либо исключения. Это не спецификация криптографической безопасности и не
подтверждение production-готовности bootstrap.

Публичная поверхность намеренно разделена:

| Слой | Заголовок | Назначение |
| --- | --- | --- |
| Stable | `m2424/m2424.hpp` | CKKS-адаптер, вычисления, контроль точности, планирование, отчёты. |
| Research | `m2424/experimental.hpp` | Bootstrap, DFT, EvalMod, Mod1 и диагностические модели. |

## 1. Базовый жизненный цикл CKKS

### 1.1 Создание контекста и ключей

**Вход:** `CkksProfile { poly_modulus_degree, coeff_modulus_bits, scale, slots }`.

1. `SealAdapter::create(profile)` проверяет форму профиля и строит SEAL context,
   CKKS encoder и evaluator.
2. `keygen(true, true)` создаёт public/secret, relinearization и полный набор
   Galois-ключей. Альтернатива `keygen(rotation_steps, need_relin)` создаёт только
   ключи для переданных поворотов.
3. Возвращается готовый `SealAdapter`; далее его нельзя копировать, но можно
   перемещать.

**Ошибочные ветки:** невалидный профиль или параметры SEAL приводят к исключению;
операции без нужного ключа завершаются `runtime_error` (public для encrypt,
secret для decrypt, relin для умножения, Galois для rotate).

### 1.2 Кодирование, шифрование и расшифрование

**Вход:** вещественный `vector<double>` или комплексный `vector<complex<double>>`,
длина не больше физического числа CKKS slots.

`encode` / `encode_complex` → `Plain` на уровне первого modulus → `encrypt` →
`Cipher`. Обратный путь: `decrypt` → `Plain` → `decode` / `decode_complex`.

Варианты `*_like` кодируют plaintext на том же уровне и масштабе, что переданный
ciphertext; `encode_complex_at_scale_like` и `encode_scalar_at_scale_like`
принимают явный положительный конечный scale. Это обязательный путь перед
операциями с ciphertext, находящимся не на первом уровне цепочки.

**Результат:** контейнеры значений всегда имеют физическую длину CKKS encoder;
потребитель должен использовать только свой payload-size. Невалидный scale,
отсутствующий encoder или ciphertext из другого context дают исключение.

### 1.3 Гомоморфные операции

| Сценарий | Вход → обработка | Выход / ограничения |
| --- | --- | --- |
| Сложение/вычитание | `Cipher, Cipher` → SEAL `add/sub` | Cipher на совместимых уровне и scale. Выравнивание не выполняется автоматически. |
| Операция с plaintext | `Cipher, Plain` → `add_plain`, `sub_plain`, `mul_plain` | Plain обязан соответствовать `parms_id` и scale ciphertext. |
| Умножение с rescale | `Cipher, Plain` → `mul_plain_rescale`; `Cipher, Cipher` → multiply, relinearize, rescale | Потребляет один уровень; ciphertext×ciphertext требует relin keys. |
| Явный rescale | `rescale_to_next(cipher)` | Переходит к следующему modulus; на конце цепочки ошибка. |
| Переход уровня | `mod_switch_to_next_preserve_scale`, `mod_switch_to(cipher,target)` | Меняет `parms_id`, не восстанавливает потерянную глубину. |
| Совмещение | `match_level_and_scale(cipher,target)` | Сначала приводит уровень, затем допускает только малую относительную коррекцию scale; большая разница — ошибка. |
| Поворот | `rotate(cipher, steps)` | Нужен соответствующий Galois key. |

`info`, `scale`, `chain_index`, `coeff_modulus_size`, `coeff_modulus_log2` и
`serialized_size` не меняют ciphertext и возвращают его метаданные.

### 1.4 Сериализованный клиент–вычислитель

**Вход:** public key, optional relin/Galois keys и ciphertext в `SerializedBuffer`.

1. Клиент вызывает `save_public_key`, `save_relin_keys`, `save_galois_keys`,
   `save_cipher`; secret key остаётся у клиента.
2. Вычислитель создаёт совместимый adapter, вызывает `load_*`, выполняет доступные
   публичные вычисления и возвращает `save_cipher(result)`.
3. Клиент делает `load_cipher`, `decrypt`, `decode`.

`load_*` валидирует буфер через SEAL и не предоставляет secret key там, где он
не загружен. Неподходящий/повреждённый буфер либо попытка decrypt без secret key
завершаются исключением.

### 1.5 Низкоуровневый ModRaise

`mod_raise_to_first` и `bootstrap_modup_to_first` принимают CKKS ciphertext ниже
первого уровня, проверяют, что исходная modulus-chain является префиксом целевой,
и расширяют RNS-представление до первого уровня. Второй метод поддерживает
`CenteredLift` и `UncenteredLift`.

**Важно:** это structural RNS operation, не гомоморфное вычисление и не
value-preserving refresh. `unsafe_reinterpret_scale_for_diagnostics` (а также
устаревший `multiply_decoded_value`) меняет интерпретацию scale только для
диагностики; его нельзя использовать как умножение plaintext.

## 2. Проверяемый вычислительный pipeline

### 2.1 Accuracy и ABFT

`compare(expected, actual, tolerance)` требует равные размеры, вычисляет
`max_abs_error`, `mean_abs_error` и `within_tolerance`. Несовпадение размеров
или недопустимый tolerance — исключение.

ABFT-путь: `append_checksum(values)` добавляет checksum к payload;
`verify_appended_checksum(values, payload_size, tolerance)` пересчитывает её и
возвращает `ChecksumResult { expected, actual, error, ok }`. Это защита от
несогласованности вычисления, а не криптографическая аутентификация результата.

### 2.2 CheckedEvaluator

**Вход:** adapter, размер payload, tolerance, ciphertext/ plaintext и ожидаемый
открытый результат для каждой операции.

Каждый метод (`add`, `sub`, `add_plain`, `sub_plain`, `mul`,
`mul_plain_rescale`, `rescale_to_next`, `rotate`, `sum_slots`,
`linear_transform`) выполняет соответствующую операцию adapter, затем `finalize`:

1. расшифровывает результат;
2. берёт первые `payload_size` значений;
3. сравнивает их с `expected`;
4. возвращает `CheckedResult { operation, cipher, info, accuracy, ok }`;
5. записывает расход в `OperationBudgetBuilder`.

Таким образом `ok=false` — это нормальный диагностический результат при ошибке
точности; инфраструктурные ошибки (ключи, уровни, несовместимые ciphertext) всё
ещё передаются исключениями. `reset_operation_budget` очищает только статистику,
не меняет ciphertext.

### 2.3 Budget и parameter planning

`OperationBudgetBuilder` накапливает счётчики операций. Функции
`estimate_linear_transform_budget` и `estimate_sum_slots_budget` строят budget
без выполнения HE-операций.

`plan_ckks_parameters(request)` обрабатывает сценарии:

1. Проверяет `target_error`, multiplicative depth, slots и security level.
2. Выбирает loss-bits из `operation_profile` либо фактического budget
   (`use_operation_budget=true`), с возможным явным override.
3. Перебирает допустимые `poly_modulus_degree`, work-bits и уровни, сверяя сумму
   coeff moduli с ограничением SEAL для запрошенной security.
4. Возвращает `CkksPlanningResult`: профиль, оценку ошибки и `passes_target_error`.

Невозможный набор требований не маскируется: возвращается `invalid_argument`
(например security вне 128/192/256, слишком много slots, недостижимая точность).
`CheckedEvaluator::plan_refresh_for_tracked_budget` передаёт накопленный budget
в bootstrap planner, не запускает refresh.

## 3. Линейная алгебра и полиномы

### 3.1 LinearTransform и sum_slots

`LinearTransform(terms)` проверяет, что terms не пусты, rotations неотрицательны
и коэффициенты имеют согласованную длину. `apply(adapter,input)` для каждого term:

1. вращает input (либо использует input при rotation=0);
2. кодирует коэффициенты на уровне rotated ciphertext;
3. выполняет `mul_plain_rescale`;
4. выравнивает и суммирует частичные результаты.

Возвращается ciphertext transform-а. Требуются Galois keys всех ненулевых
`rotation_steps()`. `sum_slots` применяет повороты степеней двойки и складывает
слоты; `slot_count` обязан быть ненулевой степенью двойки.

### 3.2 Diagonal transform и DFT

`DiagonalLinearTransform::from_matrix` принимает квадратную конечную complex
матрицу, выделяет ненулевые cyclic diagonals и строит `DiagonalTerm`.
`apply_plain` реализует то же преобразование в открытом виде; `apply` и
`apply_at_plain_scale` выполняют rotate → complex plaintext multiply → sum.

`BootstrapDftPlan` и `FactorizedLinearTransform` строят DFT-like слои для
`HomomorphicEncode`/`HomomorphicDecode`, рассчитывают нужные rotations и
применяют слои последовательно. В encrypted-ветке каждый слой потребляет
ресурсы chain; отсутствие Galois keys, не-power-of-two slots или нехватка уровней
завершают путь исключением.

### 3.3 PolynomialEvaluator

Конструктор принимает разреженный список `(degree, coefficient)`, требует
уникальные степени, конечные коэффициенты и хотя бы один ненулевой term.
Поддерживаются только положительные степени: constant-only polynomial отклоняется.

`evaluate(adapter,input)` строит нужные степени через последовательные
ciphertext multiplications, умножает их на plaintext coefficients, выравнивает и
складывает. Результат — зашифрованное значение полинома; глубина и relin keys
должны быть предусмотрены профилем.

## 4. Профили и отчёты

`profiles::*` возвращают преднастроенные CKKS profiles; `by_name` выдаёт
`invalid_argument` для неизвестного имени. `describe_profile` формирует строковые
и числовые метрики; `to_csv_row` сериализует их в CSV.

`analyze_security(name, profile)` суммирует modulus bits, сравнивает с SEAL
лимитом для уровня security и возвращает `SecurityReport`. Ветка с неизвестным
`poly_modulus_degree` не утверждает безопасность: отчёт указывает отсутствие
поддерживаемого лимита. `project_minimum_security` агрегирует наихудший уровень.

## 5. Research: EvalMod и Mod1

### 5.1 EvalModPolynomial

Для scalar, complex и vector overload-ов вход обязан быть конечным и лежать в
интервале аппроксимации. `evaluate_plain` возвращает приближение для выбранной
степени (`P3`, `P5`, `P7`, `P3DoubleAngle`); `sine_reference` даёт открытый
reference.

Encrypted `evaluate(adapter,cipher,degree)` последовательно формирует степени,
масштабирует и rescale-ит промежуточные ciphertext. Он требует достаточной
глубины и scale capacity. При нехватке уровня/невозможном align scale выбрасывает
исключение, а не возвращает неточный «успех».

### 5.2 Mod1 approximation и Mod1Circuit

`make_mod1_approximation(model)` валидирует degree, input bound, scale и
модель. Есть две ветки:

| `BootstrapMod1Type` | Plain evaluation | Encrypted evaluation |
| --- | --- | --- |
| `LegacySineP3` | Реальный P3 sine-based путь | Доступен через `Mod1Circuit`. |
| `CosDiscrete` | Экспериментальный sine-polynomial placeholder | Не является настоящей CosDiscrete реализацией; availability зависит от модели. |

`make_wide_mod1_approximation` допускает только `CosDiscrete` degree ≥ 15,
но также строит approximation model, а не доказанную high-degree CosDiscrete
схему. `evaluate_*_plain` возвращает открытое приближение и валидирует domain.

`Mod1Circuit::evaluate_with_report` возвращает `Mod1EncryptedEvaluation` с
входным/выходным chain-index, consumed levels и output scale. Если encrypted
ветка для модели недоступна, `evaluate` выбрасывает исключение. Это важнее, чем
молчаливо подменить encrypted результат plaintext вычислением.

## 6. Research: планирование bootstrap

### 6.1 Pipeline и rotation plans

`make_research_bootstrap_plan(slots)` и `make_scalable_bootstrap_plan(slots)`
создают декларативный `BootstrapPipelinePlan`: порядок цепочки, backend,
period/scaling model, EvalMod degree и stage contracts. `bootstrap_plan_rotation_steps`
возвращает минимально необходимые rotations.

План не выполняет bootstrapping и не доказывает сохранение значения. Его
`active_gate` указывает, до какой стадии применимы текущие предпосылки.

### 6.2 Scale gate, layout и параметры

Входы scale-planner-а включают метаданные ciphertext, активные modulus bits,
ожидаемую амплитуду, target scale, минимальный остаток уровней и EvalMod margin.

Сценарий `plan_bootstrap_refresh_scale_gate`:

1. вычисляет bootstrap period и normalisation factor;
2. проверяет representability plaintext scalar;
3. планирует decomposition/scaling squash и расход уровней;
4. проверяет capacity первой EvalMod multiplication;
5. возвращает `BootstrapScaleDesign` со статусом `ReadyForEvalModP3` либо одним
   из `PeriodModelBlocked`, `ScaleStrategyBlocked`, `EvalModCapacityBlocked`.

`search_bootstrap_refresh_scale_gate` перебирает кандидаты и возвращает лучший
готовый план или blocker. `plan_bootstrap_layout` и `plan_bootstrap_parameters`
совмещают slots, security, transform-depth и scale constraints; вместо
необоснованного профиля возвращают typed blocked status.

### 6.3 Refresh decision

`plan_bootstrap_refresh(request)` сравнивает доступный chain-index с budget
следующего вычислительного блока:

- `ComputeFitsWithoutRefresh` — refresh не требуется;
- `RefreshRequired` — budget не помещается, но есть план refresh;
- `RefreshPlanBlocked` — refresh нужен, но layout/scale/security не позволяет
  безопасно запланировать его.

Возвращаемый статус — решение planner-а, а не выполнение ciphertext операции.

## 7. Research: bootstrap execution

### 7.1 Bootstrapper

`analyze_depth(input,max_steps)` шифрует вход, многократно возводит ciphertext
в квадрат до исчерпания уровней/ошибки, и возвращает `BootstrapReport` с
метриками границы и причиной остановки. Это измерительный сценарий.

`refresh*` строят `BootstrapPrototype` и запускают одну из веток:

| API | Порядок | Контроль результата |
| --- | --- | --- |
| `refresh` | legacy prototype | нет expected oracle; отчёт стадий. |
| `refresh_checked` | legacy prototype | сверяет с `ComplexVector expected`. |
| `refresh_slots_to_coeffs_first*` | scalable SlotsToCoeffsFirst + FftLike | optional expected oracle. |
| `*_guarded` | сначала `plan_refresh_for_budget` | пропускает refresh при `ComputeFitsWithoutRefresh`; иначе запускает или возвращает blocker. |

Обязательные предусловия: `slots` — power of two, expected имеет ровно slots
значений, adapter получил все rotations/relin keys, в chain достаточно уровней.
В checked варианте `BootstrapPrototypeReport` содержит `preserve_value_criterion`,
`restore_level_criterion`, остаток уровней и статусы каждой стадии. `FAIL` или
`BLOCKED` в отчёте — наблюдаемый результат исследовательского запуска, не
исключение сам по себе; ошибочные параметры и невозможные примитивы — исключения.

### 7.2 Внутренние prototype-маршруты

`BootstrapPrototype` конфигурируется normalisation mode, period mode, circuit
order, transform backend, EvalMod degree/policy и параметрами STC/ModUp.

| Circuit order | Выполняемые стадии | Статус |
| --- | --- | --- |
| `ModRaiseFirst` | ModRaise → CoeffToSlot → normalise → EvalMod → SlotToCoeff → denormalise | Prototype. |
| `SlotsToCoeffsFirst` | SlotToCoeff → scale-down/ModUp contract → CoeffToSlot → normalise → EvalMod → restore | Prototype; выполняется guarded API. |

Если normalized amplitude не попадает в EvalMod interval, путь добавляет
`BLOCKED` stages и возвращает отчёт/текущий ciphertext без запуска EvalMod.
Если post-refresh ModRaise включён, он даёт structural подъём, но критерий
сохранения значения всё равно проверяется отдельно.

### 7.3 STC/ModUp contract и reference model

`plan_stc_first_modup` строит scale-down план из metadata после SlotToCoeff.
`apply_stc_scale_down` и `bootstrap_scale_down_to_q` применяют его к ciphertext,
возвращая уровни, scales, ratio и note. Невозможный target modulus/scale,
последний уровень цепи или недопустимые числа приводят к исключению.

`bootstrap_stc_reference` работает только в plaintext/reference domain:
строит lattice sample, делает integer-lattice reduction и возвращает error/gain
метрики. Он нужен для сравнения семантики encrypted path; сам по себе не
производит ciphertext refresh.

## 8. Precision-model сценарии

Эти функции не исполняют HE-circuit. Они принимают measurements/calibration и
возвращают решения:

- `make_bootstrap_error_budget` делит total error по циклам и стадиям;
- `make_bootstrap_error_recurrence` учитывает amplification;
- `required_ciphertext_scale_log2` выводит требуемый scale из rotation noise;
- `decide_evalmod_small_signal` выбирает допустимость линейного пути;
- `fit_dft_precision_floor` оценивает quantization coefficient/noise floor;
- `estimate_bootstrap_dft_cost` считает layers, diagonals, rotations и rescale;
- `plan_bootstrap_precision` выбирает профиль и transform scale либо возвращает
  `feasible=false` с текстом blocker.

Все входные значения валидируются на конечность, положительность и непустые
candidate sets; это planning diagnostics, не runtime policy adapter-а.

## 9. Тестовое покрытие и эксплуатационный режим

По умолчанию CTest включает 12 быстрых тестов: accuracy, adapter failures, DFT,
layout, precision model, scale design, STC reference, Mod1 approximation/circuit,
operation budget, parameter planner и smoke. Каждый имеет timeout 60 секунд.

При `-DM2424_BUILD_RESEARCH_TESTS=ON` дополнительно собираются research-тесты
`test_bootstrap_scalable_refresh` и `test_bootstrap_multi_cycle_precision`, а
также diagnostic executables. Их запускают явно через
`ctest --test-dir build -L research --output-on-failure`; timeout установлен в
300 секунд. Эта граница не скрывает статус bootstrap: она предотвращает зависание
обычной CI-проверки из-за длительных исследовательских экспериментов.

### Единый цикл сравнения компонентов

`BootstrapExperimentConfig` объединяет заменяемые части одного prototype
варианта: circuit order, transform backend, EvalMod, normalisation, period,
STC/ModUp и scale settings. Цикл исследования: собрать configs → объединить
`bootstrap_experiment_rotation_steps` → создать ключи → вызвать
`run_bootstrap_experiments` на одном исходном ciphertext → сравнить
`BootstrapExperimentResult`. Outcome всегда один из `passed`, `blocked`,
`failed`, а `blocker` и stage report сохраняют причину. Так сравнение не зависит
от копий orchestration-кода в отдельных benchmark-ах.

## 10. Краткая карта результатов

| Тип результата | Где появляется | Как интерпретировать |
| --- | --- | --- |
| Value object | `Cipher`, `Plain`, vectors, plans, reports | Обычный успешный маршрут. |
| `ok=false` / criteria=false | accuracy, checked evaluator, bootstrap report | Вычисление завершилось, но не прошло проверку точности/контракта. |
| Typed blocked status | bootstrap planning | Осмысленно объяснённая невозможность плана. |
| `BLOCKED` stage | prototype bootstrap | Дальнейшая небезопасная стадия не запускалась. |
| Exception | validation, keys, SEAL state, levels/scales | Вызов не имеет корректного результата. |

Для пользовательского приложения рекомендуемый маршрут: выбрать профиль или
получить его из planner → `SealAdapter::create` → создать только нужные keys →
encode/encrypt → вычислить через adapter/CheckedEvaluator → при необходимости
сериализовать → decrypt/decode у владельца secret key. Research API следует
вызывать только с явной обработкой report/status/blocker и без предположения,
что ModRaise или prototype refresh автоматически сохраняют значение.
