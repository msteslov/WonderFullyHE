# WonderFullyHE

WonderFullyHE — C++17-библиотека для приближённых защищённых вычислений на CKKS поверх Microsoft SEAL.

Стабильный API подключается через:

```cpp
#include "m2424/m2424.hpp"
```

Bootstrap, EvalMod и DFT остаются исследовательскими API и подключаются отдельно:

```cpp
#include "m2424/experimental.hpp"
```

## Сборка

Клонирование с Microsoft SEAL:

```bash
git clone --recurse-submodules <repo-url>
```

Сборка библиотеки, приложений и тестов:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Минимальный пример:

```bash
./build/demo_basic
```

## Приложения

- `demo_basic` — минимальный цикл шифрования и вычисления.
- `demo_secure_stats` — зашифрованные сумма и среднее.
- `demo_checked_pipeline` — pipeline с контролем бюджета операций.
- `bench_ckks` — измерение базовых CKKS-операций.

Остальные demos и benchmarks собираются по умолчанию и служат примерами API или измерительными утилитами.

## Тесты

По умолчанию CTest запускает быстрые детерминированные инвариантные тесты.

Исследовательские bootstrap-тесты и диагностику можно включить отдельно:

```bash
cmake -S . -B build -DM2424_BUILD_RESEARCH_APPS=ON -DM2424_BUILD_RESEARCH_TESTS=ON
cmake --build build -j
ctest --test-dir build -L research --output-on-failure
```

Для сравнения bootstrap-вариантов используйте единый research runner:

```bash
cmake -S . -B build -DM2424_BUILD_RESEARCH_APPS=ON
cmake --build build --target bench_bootstrap_experiments -j
./build/bench_bootstrap_experiments
```

Варианты задаются через `BootstrapExperimentConfig`; результат каждого содержит
общий outcome, blocker и отчёт стадий.

## Пример API

```cpp
#include "m2424/m2424.hpp"

auto adapter = m2424::SealAdapter::create(m2424::profiles::basic_ckks());
adapter.keygen(true, true);

auto encrypted = adapter.encrypt(adapter.encode({1.0, 2.0, 3.0}));
auto squared = adapter.mul_relin_rescale(encrypted, encrypted);
auto decoded = adapter.decode(adapter.decrypt(squared));
```

## Документация

- `docs/architecture.md`
- `docs/ckks_parameters.md`
