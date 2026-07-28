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
