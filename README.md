## Обзор

Проект использует Microsoft SEAL в роли криптографического движка CKKS. Данная программа демонстрирует базовый функционал, который будет использоваться при добавлении новых технологии(бутстрапинг, ABFT...)

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

## Тесты

`test_smoke` — покрывает encode → encrypt → mul_relin_rescale → decrypt.

```bash
cmake -S . -B build -DBUILD_TESTING=ON
cmake --build build -j
ctest --test-dir build --output-on-failure
```

## Архитектура и API

`CkksProfile` описывает параметры схемы: `poly_modulus_degree`, битовые длины коэффициентов модуля, масштаб `scale` и лимит слотов. Обёртки `Plain` и `Cipher` скрывают `seal::Plaintext`/`seal::Ciphertext`, а `SealAdapter` управляет жизненным циклом SEAL‑контекста.

Основные методы адаптера:
- `SealAdapter::create(profile)` — конфигурирует CKKS‑контекст и encoder под заданный профиль.
- `keygen(need_relin, need_galois)` — генерирует секретный/публичный ключи и по требованию Relin/Galois наборы.
- `encode` / `decode` — преобразуют вещественный вектор в CKKS plaintext и обратно.
- `encrypt` / `decrypt` — обычные операции CKKS над plaintext/ciphertext.
- `add`, `sub`, `mul_relin_rescale`, `rotate` — гомоморфные примитивы, делегирующие в `seal::Evaluator`.

Строгие математические формулировки для каждого метода вынесены в `api.tex`.

Функция `m2424::version()` отдаёт семантическую версию библиотеки и используется в демо как sanity‑check линковки.

## Структура каталога

- `include/m2424/` — публичные заголовки адаптера.
- `src/` — реализация `SealAdapter` и версии.
- `apps/` — демо.
- `tests/` — компактные проверки корректности.
- `extern/seal/` — git submodule Microsoft SEAL; все операции делегируются туда.

## TODO: Бенчмарки, которые надо реализовать

- [ ] Throughput для `mul_relin_rescale` в зависимости от `poly_modulus_degree` и глубины цепочки.
- [ ] Латентность `rotate`/`rotate_inverse` для разных наборов Galois‑ключей.
- [ ] Профилировка `encode`/`decode` на батчах (длина сигнала, влияние `scale`).
- [ ] Снять отдельные метрики памяти для SEAL‑контекста, наборов ключей и промежуточных ciphertext.
