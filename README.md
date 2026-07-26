# WonderFullyHE

WonderFullyHE — C++17-библиотека для приближённых защищённых вычислений на CKKS поверх Microsoft SEAL.

Стабильный API подключается через:

```cpp
#include "m2424/m2424.hpp"
```

DFT, EvalMod и parameter-planning API подключаются отдельно:

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
