## Обзор

WonderFullyHE использует Microsoft SEAL в роли криптографического движка CKKS. Поверх SEAL реализуется библиотечный слой для защищённых вычислений: адаптер скрывает детали SEAL, модуль точности задаёт единые метрики ошибки, ABFT-модуль проверяет численную согласованность результатов, а benchmark-приложение собирает первичные метрики времени, точности и размера сериализованных объектов.

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
./build/demo_bootstrap_pipeline
./build/demo_profile_report
```

`demo_abft` проверяет ABFT-инварианты: для `add/sub` полезная нагрузка дополняется checksum-слотом, для `mul` checksum произведения сравнивается с CPU-эталоном, а для `rotate` проверяется сохранение суммы по всем CKKS-слотам.

`bench_ckks` печатает CSV со временем операций, ошибкой относительно CPU-эталона и сериализованными размерами ciphertext/ключей.

`demo_noise_growth` печатает CSV по последовательным зашифрованным возведениям в квадрат. Сценарий показывает, как растёт ошибка, как меняются `scale`/`chain_index`, и где заканчивается доступная мультипликативная глубина без bootstrapping.

`demo_secure_stats` показывает прикладной сценарий защищённой обработки данных: сумма и среднее считаются над зашифрованным вектором через ротации и сложения, после расшифровки результат сравнивается с CPU-эталоном.

`demo_bootstrap_pipeline` печатает отчёт bootstrapping-модуля: профиль `depth_ckks`, границу вычислительной глубины, параметры ciphertext и этапы конвейера `ModRaise -> CoeffToSlot -> EvalMod -> SlotToCoeff`.

`demo_profile_report` печатает CSV-таблицу CKKS-профилей, используемых в демо: степень полиномиального модуля, число слотов, цепочку коэффициентных модулей, суммарный размер modulus, масштаб и оценку доступной глубины умножений.

## Тесты

`test_smoke` — покрывает encode → encrypt → mul_relin_rescale → decrypt, а также add/sub/rotate, ABFT checksum, ошибки без нужных ключей и базовую валидацию профиля.

```bash
cmake -S . -B build -DBUILD_TESTING=ON
cmake --build build -j
ctest --test-dir build --output-on-failure
```

## Архитектура и API

`CkksProfile` описывает параметры схемы: `poly_modulus_degree`, битовые длины коэффициентов модуля, масштаб `scale` и лимит слотов. Обёртки `Plain` и `Cipher` скрывают `seal::Plaintext`/`seal::Ciphertext`, а `SealAdapter` управляет жизненным циклом SEAL‑контекста. В текущей C++-реализации публичные типы и функции находятся в пространстве имён `m2424`.

Подробные проектные заметки:
- `docs/architecture.md` — слои библиотеки и роль каждого модуля.
- `docs/ckks_parameters.md` — выбранные CKKS-параметры, расчёт слотов, scale, modulus chain и глубины.

Основные методы адаптера:
- `SealAdapter::create(profile)` — конфигурирует CKKS‑контекст и encoder под заданный профиль.
- `keygen(need_relin, need_galois)` — генерирует секретный/публичный ключи и по требованию Relin/Galois наборы.
- `slot_count` — возвращает фактическое число CKKS-слотов для выбранного профиля.
- `encode` / `decode` — преобразуют вещественный вектор в CKKS plaintext и обратно.
- `encrypt` / `decrypt` — обычные операции CKKS над plaintext/ciphertext.
- `add`, `sub`, `mul_relin_rescale`, `rotate` — гомоморфные примитивы, делегирующие в `seal::Evaluator`.
- `serialized_size`, `public_key_size`, `relin_keys_size`, `galois_keys_size` — вспомогательные методы для benchmark-измерений.
- `info`, `scale`, `chain_index`, `coeff_modulus_size` — диагностика состояния ciphertext для анализа глубины и подготовки bootstrapping.

Модуль `m2424::accuracy` задаёт единые метрики точности: `max_abs_error`, `mean_abs_error` и `compare(expected, actual, tolerance)`. Demo и тесты используют этот общий код, чтобы критерии корректности не расходились между сценариями.

Модуль `m2424::abft` содержит checksum-инструменты: `append_checksum`, `checksum`, `verify_appended_checksum`, `verify_checksum_value`.

Модуль `m2424::Bootstrapper` выделяет bootstrapping как отдельный компонент библиотеки. Реализация связывает диагностику вычислительной глубины с этапами CKKS bootstrapping-конвейера и фиксирует параметры ciphertext: `scale`, `chain_index`, размер ciphertext и критерии `Dec(c') ≈ Dec(c)`, `level(c') > level(c)`.

Модуль `m2424::profile_report` формирует табличное описание выбранных CKKS-параметров.

Строгие математические формулировки для каждого метода вынесены в `api.tex`.

Функция `m2424::version()` отдаёт семантическую версию библиотеки и используется в демо как sanity‑check линковки.

## Структура каталога

- `include/m2424/` — публичные C++-заголовки WonderFullyHE.
- `src/` — реализация `SealAdapter` и версии.
- `apps/` — демо.
- `tests/` — компактные проверки корректности.
- `extern/seal/` — git submodule Microsoft SEAL; все операции делегируются туда.

## Планируемое развитие

- [x] Базовый CSV benchmark для encode/encrypt/decrypt/add/mul/rotate и размеров ключей.
- [ ] Throughput для `mul_relin_rescale` в зависимости от `poly_modulus_degree` и глубины цепочки.
- [ ] Латентность `rotate`/`rotate_inverse` для разных наборов Galois‑ключей.
- [ ] Профилировка `encode`/`decode` на батчах (длина сигнала, влияние `scale`).
- [ ] Снять отдельные метрики памяти для SEAL‑контекста, наборов ключей и промежуточных ciphertext.
- [x] Расширить ABFT checksum с add/sub на mul и rotate.
- [ ] Расширить ABFT checksum на цепочки операций.
- [x] Добавить noise-growth demo перед проектированием bootstrapping-прототипа.
- [x] Добавить прикладной сценарий защищённой статистики.
- [x] Выделить начальный bootstrapping-модуль с диагностикой глубины и статусом этапов.
- [x] Добавить отчёт по CKKS-профилям и расчётным параметрам.
- [ ] Реализовать строительные блоки bootstrapping: полиномиальную аппроксимацию и rotation-based linear transforms.
