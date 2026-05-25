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
./build/demo_noise_growth
./build/demo_secure_stats
./build/demo_client_compute_roundtrip
./build/demo_galois_key_optimization
./build/demo_bootstrap_pipeline
./build/bench_bootstrap_parts
./build/bench_bootstrap_full
./build/demo_bootstrap_diagonals
./build/demo_bootstrap_prototype
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

`demo_noise_growth` печатает CSV по последовательным зашифрованным возведениям в квадрат. Сценарий показывает, как растёт ошибка, как меняются `scale`/`chain_index`, и где заканчивается доступная мультипликативная глубина без bootstrapping.

`demo_secure_stats` показывает прикладной сценарий защищённой обработки данных: сумма и среднее считаются над зашифрованным вектором через ротации и сложения, после расшифровки результат сравнивается с CPU-эталоном.

`demo_client_compute_roundtrip` показывает разделённый сценарий: один контекст генерирует ключи и шифрует данные, вычислительный контекст получает только публичный ключ, Galois-ключи и ciphertext, выполняет агрегацию без secret key, после чего результат расшифровывается отдельным контекстом с secret key.

`demo_galois_key_optimization` сравнивает полный набор Galois-ключей с ограниченным набором rotation keys для конкретного вычисления. Сценарий показывает размер ключей, время генерации, время `linear_transform`/`sum_slots` и численную ошибку.

`demo_bootstrap_pipeline` печатает отчёт bootstrapping-модуля: профиль `depth_ckks`, границу вычислительной глубины, параметры ciphertext и этапы конвейера `ModRaise -> CoeffToSlot -> EvalMod -> SlotToCoeff`.

`bench_bootstrap_parts` печатает CSV по строительным блокам bootstrapping: `mul_plain_rescale`, rotation-based `linear_transform`, `sum_slots` и `polynomial_eval`. В отчёт входят время, уровень ciphertext, ошибка и сериализованный размер результата.

`bench_bootstrap_full` измеряет полный refresh-harness в профиле `boot_ckks`: подготовку rotation steps, генерацию ключей, полный runtime, время `CoeffToSlot`, `EvalMod`, `SlotToCoeff`, финальную ошибку и статус. Диагональные преобразования выполняются через baby-step/giant-step, поэтому для dense-transform на 16 логических slots требуется 6 rotation keys вместо 15.

`demo_bootstrap_diagonals` строит комплексную матрицу канонического вложения, переводит её в диагональное разложение `sum diag_k * Rot_k(x)` и проверяет это разложение на CPU и на зашифрованном CKKS-векторе. Это первый исполняемый шаг к `CoeffToSlot`/`SlotToCoeff`.

`demo_bootstrap_prototype` связывает bootstrapping-блоки в один проверяемый refresh-harness: `mod_raise_harness -> CoeffToSlot -> EvalMod -> SlotToCoeff -> refresh_result`. Для каждого этапа печатаются уровень ciphertext, масштаб, максимальная ошибка, runtime и статус относительно tolerance `2e-5`. В текущей версии `mod_raise_harness` использует повторное шифрование в `boot_ckks`, а остальные этапы выполняются над ciphertext.

`demo_eval_mod_polynomial` проверяет полином `EvalMod` степени 7 на диапазоне `[-2^-10, 2^-10]`: сначала против `sin(2*pi*u)/(2*pi)` на открытых данных, затем на зашифрованном CKKS-векторе в профиле `boot_ckks`.

`bench_parallel_throughput` измеряет масштабирование на независимых ciphertext. Benchmark разделяет `setup_ms` и `runtime_ms`: подготовка включает создание контекстов, ключей и входных ciphertext, а runtime измеряет параллельные вычисления над уже зашифрованными данными.

`demo_precision_profiles` сравнивает готовые профили `basic_ckks`, `balanced_ckks`, `depth_ckks` и `high_precision_ckks` по времени, ошибке, состоянию ciphertext и сериализованному размеру результата.

`demo_checked_pipeline` показывает контролируемый вычислительный конвейер `add -> mul -> rotate -> sum_slots`: после каждого шага печатаются `max_abs_error`, `mean_abs_error`, `scale`, `chain_index`, размер ciphertext, ABFT-статус и общий статус.

`demo_profile_report` печатает CSV-таблицу готовых CKKS-профилей из `m2424::profiles`: степень полиномиального модуля, число доступных слотов, цепочку коэффициентных модулей, суммарный размер modulus, масштаб и оценку доступной глубины умножений.

`demo_security_report` печатает CSV-таблицу проверки профилей по лимитам Microsoft SEAL для `tc128`, `tc192` и `tc256`. Общий уровень проекта определяется минимальным уровнем среди используемых профилей.

## Тесты

`test_smoke` — покрывает encode → encrypt → mul_relin_rescale → decrypt, add/sub/rotate, plaintext-операции, сериализацию ключей/ciphertext, checked evaluator, linear transform, polynomial evaluator, `sum_slots`, ABFT checksum, ошибки без нужных ключей, базовую валидацию профиля и security report.

```bash
cmake -S . -B build -DBUILD_TESTING=ON
cmake --build build -j
ctest --test-dir build --output-on-failure
```

## Архитектура и API

`CkksProfile` описывает параметры схемы: `poly_modulus_degree`, битовые длины коэффициентов модуля, масштаб `scale` и лимит слотов. Готовые профили находятся в `m2424::profiles`, поэтому обычному пользователю не нужно вручную подбирать modulus chain для первого запуска. Обёртки `Plain` и `Cipher` скрывают `seal::Plaintext`/`seal::Ciphertext`, а `SealAdapter` управляет жизненным циклом SEAL‑контекста. В текущей C++-реализации публичные типы и функции находятся в пространстве имён `m2424`.

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
- `docs/project_status.md` — текущий статус реализации, ограничения и следующие этапы.

Основные методы адаптера:
- `SealAdapter::create(profile)` — конфигурирует CKKS‑контекст и encoder под заданный профиль.
- `keygen(need_relin, need_galois)` — генерирует секретный/публичный ключи и по требованию Relin/Galois наборы.
- `keygen(rotation_steps, need_relin)` — генерирует только нужные Galois-ключи для заданных ротаций.
- `slot_count` — возвращает фактическое число CKKS-слотов для выбранного профиля.
- `encode` / `decode` — преобразуют вещественный вектор в CKKS plaintext и обратно.
- `encode_like`, `encode_scalar_like` — кодируют plaintext на уровне и масштабе заданного ciphertext.
- `encrypt` / `decrypt` — обычные операции CKKS над plaintext/ciphertext.
- `add`, `sub`, `add_plain`, `sub_plain`, `mul_plain_rescale`, `mul_relin_rescale`, `rotate` — гомоморфные примитивы, делегирующие в `seal::Evaluator`.
- `mod_switch_to`, `match_level_and_scale` — выравнивание ciphertext перед сложением членов разных уровней.
- `serialized_size`, `public_key_size`, `relin_keys_size`, `galois_keys_size` — вспомогательные методы для benchmark-измерений.
- `save_public_key`, `save_secret_key`, `save_relin_keys`, `save_galois_keys`, `save_cipher` — сериализация ключей и ciphertext в байтовый буфер.
- `load_public_key`, `load_secret_key`, `load_relin_keys`, `load_galois_keys`, `load_cipher` — загрузка ключей и ciphertext в новый CKKS-контекст с проверкой совместимости параметров.
- `info`, `scale`, `chain_index`, `coeff_modulus_size` — диагностика состояния ciphertext для анализа глубины и подготовки bootstrapping.

Модуль `m2424::accuracy` задаёт единые метрики точности: `max_abs_error`, `mean_abs_error` и `compare(expected, actual, tolerance)`. Demo и тесты используют этот общий код, чтобы критерии корректности не расходились между сценариями.

Модуль `m2424::CheckedEvaluator` выполняет операции через `SealAdapter` и возвращает `CheckedResult`: ciphertext, `CipherInfo`, метрики точности, tolerance и статус. Он нужен для сценариев, где после каждого шага вычисления надо контролировать ошибку, уровень ciphertext и масштаб.

Модуль `m2424::DiagonalLinearTransform` строит и применяет диагональное разложение комплексной матрицы:

```text
A*x = sum_k diag_k * Rot_k(x)
```

Для bootstrapping-блоков добавлены генераторы `canonical_embedding_matrix(slots)` и `invert_matrix(matrix)`. Они дают численные коэффициенты для прототипов `CoeffToSlot` и `SlotToCoeff`; коэффициенты не выписываются вручную, а вычисляются из корней единицы CKKS.

Модуль `m2424::EvalModPolynomial` реализует полином:

```text
P7(u) = u - 6.579736267393*u^3 + 12.98787880453*u^5 - 12.20811674381*u^7
```

Рабочий диапазон первой версии: `|u| <= 2^-10`. Ciphertext-версия считает степени `u^2`, `u^3`, `u^5`, `u^7` отдельной схемой, без общего последовательного подъёма степени.

Модуль `m2424::BootstrapPrototype` собирает строительные блоки в refresh-harness. Он генерирует минимальный набор rotation steps для заданного числа slots, применяет диагональные `CoeffToSlot`/`SlotToCoeff`, вызывает `EvalModPolynomial` и возвращает отчёт по этапам: chain index, scale, max error, runtime и статус. Внутри `DiagonalLinearTransform` используется baby-step/giant-step-разложение, которое уменьшает число CKKS-ротаций при плотных диагональных матрицах.

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

Следующие инженерные задачи:

- [ ] Расширить benchmark до sweep-режима: несколько `poly_modulus_degree`, разные размеры payload и несколько повторов для усреднения времени.
- [ ] Добавить отдельные измерения для `rotate` при разных шагах и наборах Galois-ключей.
- [ ] Расширить ABFT checksum на цепочки операций, а не только на одиночные add/sub/mul/rotate.
- [ ] Снять отдельные метрики памяти для SEAL-контекста, наборов ключей и промежуточных ciphertext; текущий benchmark уже фиксирует сериализованные размеры.
- [x] Разделить контексты выполнения для демонстрации: один контекст генерирует ключи и шифрует данные, вычислительный контекст получает только публичные/evaluation-ключи и ciphertext.
- [ ] Подставить математически рассчитанные матрицы `CoeffToSlot` и `SlotToCoeff`.
- [ ] Подставить коэффициенты полинома `EvalMod` из математической модели.
- [ ] Собрать end-to-end `Bootstrapper::refresh(cipher)` поверх готовых строительных блоков.
- [ ] Добавить сравнение с OpenFHE для тех же строительных блоков и для end-to-end refresh после завершения bootstrapping.
