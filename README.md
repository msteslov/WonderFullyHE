## Обзор

WonderFullyHE использует Microsoft SEAL в роли криптографического движка CKKS. Поверх SEAL реализован библиотечный слой `m2424` для защищённых вычислений над вещественными данными: адаптер скрывает детали SEAL, модуль точности задаёт единые метрики ошибки, ABFT-модуль проверяет численную согласованность результатов, benchmark-приложения собирают метрики времени/точности/размеров, а bootstrapping-блоки готовят вычислительную основу для `ModRaise -> CoeffToSlot -> EvalMod -> SlotToCoeff`.

## Клонирование

- Клонировать вместе с submodule: `git clone --recurse-submodules <repo-url>`
- Обновите зависимость, если уже клонировали: `git submodule update --init --recursive`

После инициализации submodule сборка полностью офлайн.

## Сборка и запуск

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/demo_basic
```

`demo_basic` выводит `max_abs_error` и `mean_abs_error` между результатом гомоморфных вычислений и эталоном на CPU. В качестве входа берём синусоиду длиной 64: данный сигнал позволяет измерить накопленную ошибку на типичном профиле, где важна точность.

Дополнительные сценарии:

```bash
./build/demo_abft
./build/bench_ckks
./build/bench_chain_accuracy
./build/bench_parameter_planner
./build/demo_noise_growth
./build/demo_secure_stats
./build/demo_client_compute_roundtrip
./build/demo_galois_key_optimization
./build/demo_bootstrap_pipeline
./build/bench_bootstrap_parts
./build/bench_bootstrap_planner
./build/bench_bootstrap_full
./build/bench_bootstrap_refresh
./build/bench_bootstrap_validation
./build/demo_bootstrap_diagonals
./build/demo_bootstrap_prototype
./build/demo_bootstrap_cipher_path
./build/demo_bootstrap_end_to_end
./build/demo_mod_raise
./build/demo_eval_mod_polynomial
./build/bench_parallel_throughput
./build/demo_precision_profiles
./build/demo_checked_pipeline
./build/demo_profile_report
./build/demo_security_report
```

`demo_abft` проверяет ABFT-инварианты: для `add/sub` полезная нагрузка дополняется checksum-слотом, для `mul` checksum произведения сравнивается с CPU-эталоном, а для `rotate` проверяется сохранение суммы по всем CKKS-слотам.

`bench_ckks` печатает CSV со временем операций, ошибкой относительно CPU-эталона и сериализованными размерами ciphertext/ключей. По умолчанию используется `basic_ckks`; профиль можно передать первым аргументом, например:

```bash
./build/bench_ckks high_precision_ckks
```

`bench_chain_accuracy` проводит контролируемый sweep точности: при фиксированных входе, `N`, операции и масштабе меняет длину chain, а отдельно проверяет влияние `scale_log2` и битности рабочих модулей. Текущий вывод: длина chain задаёт глубину, а точность в основном задаётся согласованной парой `scale_log2 ~= work_modulus_bits`.

`bench_parameter_planner` строит профиль через `CkksParameterPlanner`, запускает выбранный repeated-mul профиль на реальном SEAL-прогоне и печатает отдельные планы по explicit operation budget. Это sanity-check для калиброванной модели, а не формальное доказательство для произвольного pipeline.

`demo_noise_growth` печатает CSV по последовательным зашифрованным возведениям в квадрат. Сценарий показывает, как растёт ошибка, как меняются `scale`/`chain_index`, и где заканчивается доступная мультипликативная глубина без bootstrapping.

`demo_secure_stats` показывает прикладной сценарий защищённой обработки данных: сумма и среднее считаются над зашифрованным вектором через ротации и сложения, после расшифровки результат сравнивается с CPU-эталоном.

`demo_client_compute_roundtrip` показывает разделённый сценарий: один контекст генерирует ключи и шифрует данные, вычислительный контекст получает только публичный ключ, Galois-ключи и ciphertext, выполняет агрегацию без secret key, после чего результат расшифровывается отдельным контекстом с secret key.

`demo_galois_key_optimization` сравнивает полный набор Galois-ключей с ограниченным набором rotation keys для конкретного вычисления. Сценарий показывает размер ключей, время генерации, время `linear_transform`/`sum_slots` и численную ошибку.

`demo_bootstrap_pipeline` печатает отчёт bootstrapping-модуля: профиль `depth_ckks`, границу вычислительной глубины, параметры ciphertext и этапы конвейера `ModRaise -> CoeffToSlot -> EvalMod -> SlotToCoeff`.

`bench_bootstrap_parts` печатает CSV по строительным блокам bootstrapping: `mul_plain_rescale`, rotation-based `linear_transform`, `sum_slots` и `polynomial_eval`. В отчёт входят время, уровень ciphertext, ошибка и сериализованный размер результата.

`bench_bootstrap_planner` демонстрирует связку `CheckedEvaluator -> OperationBudgetBuilder -> plan_bootstrap_refresh`: выполняет короткий проверяемый pipeline, получает его budget и печатает, поместится ли повторный блок без refresh или нужен bootstrapping.

`bench_bootstrap_scaling` является обязательным gate перед full refresh: он выполняет только `encrypt -> lower level -> ModRaise -> CoeffToSlot -> normalization -> decrypt/check` и проверяет, представим ли normalization scalar при выбранных `period_mode` и `plain_scale_log2`. Tiny scalar больше не выполняется одним plaintext multiply: scaling layer раскладывает его на несколько `mul_plain_rescale` шагов по доступной modulus capacity.

`bench_bootstrap_prescaled_coeff_to_slot` проверяет перенос части normalization в `CoeffToSlot` diagonals. Prescaled diagonals кодируются через `apply_at_plain_scale`, чтобы избежать transparent plaintext; такой путь физически исполняется, но высокий plaintext scale поднимает ciphertext scale и пока оставляет blocker `not_enough_levels_for_scale`.
Benchmark также выполняет физический участок `CoeffToSlot -> normalization -> ScaleSquash` и печатает `p3_ready`; на текущем `boot_ckks` prescale cases пока дают `0` готовых входов для `EvalMod`.

`bench_bootstrap_full` и `bench_bootstrap_refresh` оставлены как experimental harness для разработки конвейера. `bench_bootstrap_refresh` сначала печатает operation-budget gate, затем scale gate, и не запускает experimental refresh, если следующий блок помещается без refresh или текущий prototype-refresh заблокирован. Эти harness не являются доказательством корректного bootstrapping, пока `bench_bootstrap_scaling` блокирует full validation.

`bench_bootstrap_validation` сначала запускает внутренний scaling gate. Если нет real period-mode с `evalmod_ready=true`, benchmark печатает `blocked_by_evalmod_ready_scaling` или `blocked_by_period_model` и не выполняет `EvalMod`, `SlotToCoeff` и post-refresh continuation. После прохождения gate full validation запускается только на выбранном scaling mode, а не на полном diagnostic sweep.

`bench_bootstrap_one_case` печатает узкий trace для одного фиксированного bootstrapping-сценария после scaling gate. В отчёт входят chunks/levels для decomposed normalization и denormalization, `chain_remaining_before_evalmod` и `chain_remaining_after_evalmod`; `NoBootstrapPeriod` используется только как diagnostic baseline.

`bench_bootstrap_reference_path` является E0 oracle-harness перед подключением full refresh. Для маленьких `slots=4/8/16` он строит CPU/reference path рядом с ciphertext path и печатает stage/domain/level/scale/magnitude/error/blocker. `ModRaise` в нём считается только structural stage; после него correctness проверяется по явным reference-ожиданиям каждого следующего этапа, а не по сохранению исходного значения.

`bench_bootstrap_period_model` проверяет только period/scaling model без запуска `EvalMod`: scalar correctness, попадание в интервал EvalMod, оставшиеся уровни и scale readiness перед P3. Benchmark сравнивает обычный `boot_ckks` и diagnostic `boot_deep_ckks`; `boot_deep_ckks` нужен только для проверки scale/levels budget и не означает, что full bootstrap готов. `ScaleSquash` использует обычный CKKS rescale, чтобы проверить, хватает ли уровней снизить scale до EvalMod-safe диапазона без изменения значения.

`bench_bootstrap_scale_strategy` является быстрым algebraic planner для того же blocker: он считает required modulus drop, доступные levels, минимальные scalar chunks, дополнительные `ScaleSquash` levels и capacity первого EvalMod multiplication до тяжёлого ciphertext-прогона.

`BootstrapScaleDesign` фиксирует scaling decision как кодовый gate, а не как ручную интерпретацию таблицы: `period_model_blocked`, `scale_strategy_blocked`, `evalmod_capacity_blocked` или `ready_for_evalmod_p3`. Подключать режим в `BootstrapPrototype` можно только после `ready_for_evalmod_p3`.

`BootstrapPeriodFeasibilityWindow` проверяет совместимость magnitude и EvalMod scale capacity до выбора конкретного периода. Для текущего dense/research path `bench_bootstrap_scale_strategy` показывает отрицательное окно: период, достаточный для попадания normalized slots в EvalMod interval, больше максимального периода, при котором P3 ещё имеет scale capacity.

`bench_bootstrap_prescaled_coeff_to_slot` проверяет гипотезу переноса части normalization прямо в plaintext diagonals `CoeffToSlot`. Сейчас benchmark использует `apply_at_plain_scale`: prescaled diagonals становятся представимыми, но высокий plaintext scale поднимает scale результата и пока оставляет blocker `not_enough_levels_for_scale`.

`bench_bootstrap_profile_budget` проверяет гипотетические modulus-chain layouts без SEAL keygen: security budget, оставшийся modulus под EvalMod, максимум доступного drop и capacity первого EvalMod multiplication. Его задача — отделить “плохо выбран профиль” от “period model несовместим с EvalMod scale budget”. Для P3 он также печатает `required_period_offset_log2 = start_scale_log2 + target_scale_log2 + margin`; при `target_scale_log2=60` текущий offset `44` математически недостаточен, нужен offset порядка `102`.

`demo_bootstrap_diagonals` строит комплексную матрицу канонического вложения, переводит её в диагональное разложение `sum diag_k * Rot_k(x)` и проверяет это разложение на CPU и на зашифрованном CKKS-векторе. Это первый исполняемый шаг к `CoeffToSlot`/`SlotToCoeff`.

`demo_bootstrap_prototype` связывает bootstrapping-блоки в experimental refresh-harness. Текущая архитектура считает корректным только structural/scaling gate; full refresh path не считается рабочим, пока scaling gate не проходит.

`demo_mod_raise` проверяет низкоуровневый CKKS ModRaise: ciphertext после снижения уровня расширяется обратно к первой RNS-базе без расшифрования. В отчёт входят `chain_index`, число коэффициентных модулей, масштаб, ошибка декодирования относительно ciphertext до подъёма и статус.

`demo_bootstrap_cipher_path` запускает experimental ciphertext-only путь для исследования состояния ciphertext. Его вывод не является доказательством сохранения значения.

`demo_bootstrap_end_to_end` оставлен как historical experimental demo и не включён в default CTest, потому что full refresh сейчас заблокирован scaling gate.

`demo_eval_mod_polynomial` проверяет полиномы `EvalMod` на диапазоне `[-2^-10, 2^-10]`: сначала против `sin(2*pi*u)/(2*pi)` на открытых данных, затем на зашифрованном CKKS-векторе в профиле `boot_ckks`.

`bench_parallel_throughput` измеряет масштабирование на независимых ciphertext. Benchmark разделяет `setup_ms` и `runtime_ms`: подготовка включает создание контекстов, ключей и входных ciphertext, а runtime измеряет параллельные вычисления над уже зашифрованными данными.

`demo_precision_profiles` сравнивает готовые профили `basic_ckks`, `balanced_ckks`, `depth_ckks` и `high_precision_ckks` по времени, ошибке, состоянию ciphertext и сериализованному размеру результата.

`demo_checked_pipeline` показывает контролируемый вычислительный конвейер `add -> mul -> rotate -> sum_slots`: после каждого шага печатаются `max_abs_error`, `mean_abs_error`, `scale`, `chain_index`, размер ciphertext, ABFT-статус и общий статус.

`demo_profile_report` печатает CSV-таблицу готовых CKKS-профилей из `m2424::profiles`: степень полиномиального модуля, число доступных слотов, цепочку коэффициентных модулей, суммарный размер modulus, масштаб и оценку доступной глубины умножений.

`demo_security_report` печатает CSV-таблицу проверки профилей по лимитам Microsoft SEAL для `tc128`, `tc192` и `tc256`. Общий уровень проекта определяется минимальным уровнем среди используемых профилей.

## Тесты

Default CTest запускает только тесты-инварианты. Demos и benchmarks собираются как executables, но не входят в default CTest: их нужно запускать явно, когда нужны CSV/измерения.

- `test_smoke` покрывает encode → encrypt → mul_relin_rescale → decrypt, add/sub/rotate, plaintext-операции, сериализацию ключей/ciphertext, checked evaluator, linear transform, polynomial evaluator, `sum_slots`, ABFT checksum, базовую валидацию профиля и security report.
- `test_accuracy` проверяет finite/NaN/Inf handling, worst-index и расширенные error metrics.
- `test_adapter_failures` проверяет negative API cases: некорректные inputs, отсутствие ключей, пустые/corrupt buffers и unsafe scale reinterpretation preconditions.

```bash
cmake -S . -B build -DBUILD_TESTING=ON
cmake --build build -j
ctest --test-dir build --output-on-failure
```

## Архитектура и API

`CkksProfile` описывает параметры схемы: `poly_modulus_degree`, битовые длины коэффициентов модуля, масштаб `scale` и лимит слотов. Готовые профили находятся в `m2424::profiles` и остаются пресетами для demo/tests; целевой путь развития — planner, который строит профиль из `target_error`, глубины, числа слотов и требования к security. Обёртки `Plain` и `Cipher` скрывают `seal::Plaintext`/`seal::Ciphertext`, а `SealAdapter` управляет жизненным циклом SEAL‑контекста.

Готовые профили:

| Профиль | Назначение |
|---|---|
| `profiles::fast_demo_ckks()` | Быстрая локальная проверка и простые эксперименты. |
| `profiles::basic_ckks()` | Основной профиль для demo, ABFT, защищённой статистики и базовых benchmark. |
| `profiles::balanced_ckks()` | Дополнительный уровень умножения на `N = 8192` без перехода на более тяжёлый `N = 16384`. |
| `profiles::depth_ckks()` | Анализ глубины и bootstrapping-блоки. |
| `profiles::high_precision_ckks()` | Более высокий масштаб `2^50` для сценариев, где точность важнее скорости. |
| `profiles::boot_ckks()` | Профиль для экспериментов с полным refresh: длинная цепочка `60-40-40-40-40-40-40-40-60` при `N = 16384`. |

Минимальный пример:

```cpp
#include "m2424/m2424.hpp"

auto adapter = m2424::SealAdapter::create(m2424::profiles::basic_ckks());
adapter.keygen(true, true);

auto encrypted = adapter.encrypt(adapter.encode({1.0, 2.0, 3.0}));
auto squared = adapter.mul_relin_rescale(encrypted, encrypted);
auto decoded = adapter.decode(adapter.decrypt(squared));
```

Если профиль выбирается из конфигурации или аргументов командной строки, используйте `m2424::profiles::by_name("high_precision_ckks")`.

Подробные проектные заметки:
- `docs/architecture.md` — слои библиотеки и роль каждого модуля.
- `docs/ckks_parameters.md` — выбранные CKKS-параметры, расчёт слотов, scale, modulus chain и глубины.
- `docs/m2424_math_model.tex` / `docs/m2424_math_model.pdf` — математическая модель библиотеки: входные объекты, CKKS-операции, ABFT-инварианты, линейные преобразования и refresh-harness.
- `docs/project_status.md` — текущий статус реализации, ограничения и следующие этапы.

Основные методы адаптера:
- `SealAdapter::create(profile)` — конфигурирует CKKS‑контекст и encoder под заданный профиль.
- `keygen(need_relin, need_galois)` — генерирует секретный/публичный ключи и по требованию Relin/Galois наборы.
- `keygen(rotation_steps, need_relin)` — генерирует только нужные Galois-ключи для заданных ротаций.
- `slot_count` — возвращает фактическое число CKKS-слотов для выбранного профиля.
- `encode` / `decode` — преобразуют вещественный вектор в CKKS plaintext и обратно.
- `encode_like`, `encode_scalar_like`, `encode_scalar_at_scale_like` — кодируют plaintext на уровне заданного ciphertext; для скаляра можно использовать масштаб ciphertext или явно заданный масштаб.
- `encrypt` / `decrypt` — обычные операции CKKS над plaintext/ciphertext.
- `add`, `sub`, `add_plain`, `sub_plain`, `mul_plain`, `mul_plain_rescale`, `mul_relin_rescale`, `rotate` — гомоморфные примитивы, делегирующие в `seal::Evaluator`.
- `mod_switch_to`, `match_level_and_scale` — выравнивание ciphertext перед сложением членов разных уровней.
- `mod_raise_to_first` — расширение CKKS ciphertext с текущего уровня к первой RNS-базе modulus chain без расшифрования.
- `multiply_decoded_value` — управляемое изменение CKKS scale для умножения декодируемого значения на положительный коэффициент без расхода уровня.
- `serialized_size`, `public_key_size`, `relin_keys_size`, `galois_keys_size` — вспомогательные методы для benchmark-измерений.
- `save_public_key`, `save_secret_key`, `save_relin_keys`, `save_galois_keys`, `save_cipher` — сериализация ключей и ciphertext в байтовый буфер.
- `load_public_key`, `load_secret_key`, `load_relin_keys`, `load_galois_keys`, `load_cipher` — загрузка ключей и ciphertext в новый CKKS-контекст с проверкой совместимости параметров.
- `info`, `scale`, `chain_index`, `coeff_modulus_size` — диагностика состояния ciphertext для анализа глубины и подготовки bootstrapping.

Модуль `m2424::accuracy` задаёт единые метрики точности: `max_abs_error`, `mean_abs_error` и `compare(expected, actual, tolerance)`. Demo и тесты используют этот общий код, чтобы критерии корректности не расходились между сценариями.

Для целевой точности `1e-9` библиотека должна выбирать параметры расчётом:

```text
required_result_bits = ceil(-log2(target_error))
work_bits = required_result_bits + calibrated_loss_bits(ops_profile)
scale_log2 ~= work_bits
work_levels >= multiplicative_depth
```

В текущих экспериментах для двух последовательных `mul_relin_rescale` потеря составляет около 14 бит, поэтому для `1e-9` минимальный практический режим — `scale_log2 = 45` и 45-битные рабочие модули. 50-битный режим остаётся conservative-вариантом.

Модуль `m2424::parameter_planner` реализует первую версию этого расчёта: `plan_ckks_parameters(request)` возвращает минимальный `CkksProfile`, выбранные биты, глубину, `SecurityReport`, `estimated_abs_error_bound` и флаг `passes_target_error`. Для нетривиальных программ нужно передавать `CkksOperationBudget` по всем используемым примитивам (`add/sub`, plaintext ops, `mul`, `rescale`, `mod_switch`, rotations, linear transforms, EvalMod, refresh), иначе planner использует только coarse `operation_profile`.

Модуль `m2424::CheckedEvaluator` выполняет операции через `SealAdapter` и возвращает `CheckedResult`: ciphertext, `CipherInfo`, метрики точности, tolerance и статус. Он нужен для сценариев, где после каждого шага вычисления надо контролировать ошибку, уровень ciphertext и масштаб.

Модуль `m2424::operation_budget` собирает бюджет операций для planner-а. `CheckedEvaluator` автоматически накапливает budget для `add/sub`, `mul`, `rotate`, `sum_slots` и `linear_transform`; прямые низкоуровневые вызовы `SealAdapter` пока считаются вручную через `OperationBudgetBuilder`.

Модуль `m2424::DiagonalLinearTransform` строит и применяет диагональное разложение комплексной матрицы:

```text
A*x = sum_k diag_k * Rot_k(x)
```

Для bootstrapping-блоков добавлены генераторы `canonical_embedding_matrix(slots)` и `invert_matrix(matrix)`. Они дают численные коэффициенты для прототипов `CoeffToSlot` и `SlotToCoeff`; коэффициенты не выписываются вручную, а вычисляются из корней единицы CKKS.

Модуль `m2424::EvalModPolynomial` реализует нечётные варианты `P3`, `P5` и `P7`. Для текущего интервала `|u| <= 2^-10` default должен оставаться `P3`: его ошибка аппроксимации около `1e-14`, что достаточно для target `1e-9` и дешевле по depth, чем `P5/P7`.

```text
P7(u) = u - 6.579736267393*u^3 + 12.98787880453*u^5 - 12.20811674381*u^7
```

Chebyshev/minimax-полиномы стоит использовать как offline-инструмент для генерации минимальной степени при изменении интервала или tolerance, а не как более тяжёлый default.

Модуль `m2424::Bootstrapper` пока является experimental точкой входа для refresh. Рабочим regression gate считается `bench_bootstrap_scaling`; full refresh blocked, если normalization scalar не представим после `ModRaise -> CoeffToSlot`.

Модуль `m2424::BootstrapPipelinePlan` задаёт целевую архитектуру bootstrapping: домены значений между стадиями, transform backend (`DenseDiagonal` для текущего research path или `FftLike` для будущего масштабируемого path), scaling strategy и active gate. Bootstrapping должен выполняться только по рассчитанному плану, который заранее проверяет `chain_index`, `scale_log2`, остаток modulus, interval для `EvalMod` и error budget.

`plan_bootstrap_refresh` связывает operation budget и bootstrapping decision: по текущему `CipherInfo` и бюджету следующего блока он возвращает `compute_fits_without_refresh`, `refresh_required` или `refresh_plan_blocked`.

`Bootstrapper::plan_refresh_for_budget` выполняет тот же gate для реального ciphertext: берёт `CipherInfo` из адаптера и возвращает решение до запуска experimental refresh.

`plan_bootstrap_refresh_scale_gate` проверяет уже собранные stage facts для prototype-refresh: period, active modulus bits, normalization levels и capacity первого `EvalMod` multiplication.

`search_bootstrap_refresh_scale_gate` перебирает допустимые `period_mode`, manual period, plaintext scale и target scale, возвращая лучший feasible design или ближайший blocker. Scale diagnostics включают `missing_drop_log2`, `missing_scalar_levels` и `missing_total_levels`, чтобы выбирать следующий scale-management шаг без слепого запуска full refresh.

Search учитывает `coeff_to_slot_prescale_log2` и `coeff_to_slot_plain_scale_log2`: prescale без физической цены больше не считается готовым design, потому что plaintext scale transform-а напрямую влияет на capacity первого `EvalMod` multiplication.

Модуль `m2424::BootstrapPrototype` оставлен как внутренний проверочный harness. Текущий `CoeffToSlot/SlotToCoeff` строится через dense diagonal transform и годится как reference для маленьких `slots`; production path должен перейти на настоящий factorized FFT-like backend с sparse layers, BSGS/hoisting и ограниченным набором rotation keys.

Модуль `m2424::abft` содержит checksum-инструменты: `append_checksum`, `checksum`, `verify_appended_checksum`, `verify_checksum_value`.

Модуль `m2424::LinearTransform` применяет линейные преобразования вида `sum_i a_i * rotate(ct, k_i)`. Он нужен для rotation-based блоков `CoeffToSlot` и `SlotToCoeff`. Отдельная функция `sum_slots` считает сумму заданного числа слотов и помещает результат в первый слот.

Модуль `m2424::PolynomialEvaluator` вычисляет полином от ciphertext по степенному базису. Он используется как программная основа для этапа `EvalMod`.

Модуль `m2424::Bootstrapper` выделяет bootstrapping как отдельный компонент библиотеки. Реализация связывает диагностику вычислительной глубины с этапами CKKS bootstrapping-конвейера и фиксирует параметры ciphertext: `scale`, `chain_index`, размер ciphertext и критерии `Dec(c') ≈ Dec(c)`, `level(c') > level(c)`.

Модуль `m2424::profile_report` формирует табличное описание выбранных CKKS-параметров.

Модуль `m2424::security_report` проверяет суммарный размер коэффициентного модуля относительно лимитов Microsoft SEAL для `tc128`, `tc192` и `tc256`.

Строгие математические формулировки для каждого метода вынесены в `api.tex`.

Функция `m2424::version()` отдаёт семантическую версию библиотеки и используется в демо как sanity‑check линковки.

## Структура каталога

- `include/m2424/` — публичные C++-заголовки WonderFullyHE.
- `src/` — реализация адаптера, accuracy, ABFT, linear transform, polynomial evaluator, bootstrapping-отчёта, profile/security reports и версии.
- `apps/` — демо.
- `tests/` — компактные проверки корректности.
- `.github/workflows/ci.yml` и `.gitlab-ci.yml` — CI-проверки сборки и тестов.
- `extern/seal/` — git submodule Microsoft SEAL; все операции делегируются туда.

## Планируемое развитие

Закрыто в текущей версии:

- [x] Базовый CSV benchmark для encode/encrypt/decrypt/add/mul/rotate и размеров ключей.
- [x] Расширить ABFT checksum с add/sub на mul и rotate.
- [x] Добавить noise-growth demo перед проектированием bootstrapping-прототипа.
- [x] Добавить прикладной сценарий защищённой статистики.
- [x] Выделить начальный bootstrapping-модуль с диагностикой глубины и статусом этапов.
- [x] Добавить отчёт по CKKS-профилям и расчётным параметрам.
- [x] Добавить программный отчёт по криптостойкости CKKS-профилей.
- [x] Добавить CI-конфигурации, лицензию, публичное описание проекта и `api.pdf`.
- [x] Добавить plaintext-операции для bootstrapping-блоков.
- [x] Добавить ограниченную генерацию Galois-ключей под заданные ротации.
- [x] Добавить rotation-based `LinearTransform`, `sum_slots` и `PolynomialEvaluator`.
- [x] Добавить benchmark строительных блоков bootstrapping.
- [x] Добавить демонстрацию уменьшения размера Galois-ключей при генерации только нужных ротаций.
- [x] Добавить benchmark параллельной обработки независимых ciphertext.
- [x] Вынести готовые CKKS-профили в публичный API `m2424::profiles`.
- [x] Добавить CMake alias target `m2424::m2424` для подключения библиотеки через `add_subdirectory`.
- [x] Добавить сериализацию публичного/секретного/evaluation-ключей и ciphertext.
- [x] Добавить разделённый roundtrip-сценарий: шифрование, вычисление без secret key, расшифрование результата.
- [x] Добавить низкоуровневый CKKS `ModRaise` для восстановления первой RNS-базы ciphertext.
- [x] Добавить scaling gate для проверки представимости normalization scalar после `ModRaise -> CoeffToSlot`.
- [ ] Исправить совместимость `ModRaise -> CoeffToSlot`: текущий scaling gate показывает взрыв амплитуды до порядка `1e80`.
- [ ] Вернуть full refresh validation только после прохождения scaling gate.

Следующие инженерные задачи:

- [ ] Расширить benchmark до sweep-режима: несколько `poly_modulus_degree`, разные размеры payload и несколько повторов для усреднения времени.
- [ ] Добавить отдельные измерения для `rotate` при разных шагах и наборах Galois-ключей.
- [ ] Расширить ABFT checksum на цепочки операций, а не только на одиночные add/sub/mul/rotate.
- [ ] Снять отдельные метрики памяти для SEAL-контекста, наборов ключей и промежуточных ciphertext; текущий benchmark уже фиксирует сериализованные размеры.
- [x] Разделить контексты выполнения для демонстрации: один контекст генерирует ключи и шифрует данные, вычислительный контекст получает только публичные/evaluation-ключи и ciphertext.
- [x] Подставить математически рассчитанные матрицы `CoeffToSlot` и `SlotToCoeff`.
- [x] Подставить коэффициенты полинома `EvalMod` из математической модели.
- [ ] Собрать end-to-end refresh-путь поверх готовых строительных блоков после прохождения scaling gate.
- [ ] Перенести refresh-путь из experimental prototype в стабильный публичный API после прохождения scaling gate.
- [x] Добавить benchmark scaling gate перед full refresh validation.
- [x] Добавить контролируемый benchmark `bench_chain_accuracy` для связи chain/scale/work bits и точности.
- [x] Добавить первую версию `CkksParameterPlanner`: `target_error`, depth, slots и security -> минимальный `CkksProfile`.
- [ ] Добавить строгий validation-benchmark для сохранения значения только после прохождения scaling gate.
- [ ] Добавить сравнение с OpenFHE для тех же строительных блоков и для end-to-end refresh после завершения bootstrapping.
- [ ] Добавить calibration layer `ops_profile -> calibrated_loss_bits`.
- [ ] Сделать `BootstrapPlanner` обязательным gate перед refresh.
- [ ] Реализовать factorized FFT-like backend для `CoeffToSlot/SlotToCoeff`; dense backend оставить reference-path.
