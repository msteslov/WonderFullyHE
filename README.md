# WonderFullyHE

WonderFullyHE — C++17-библиотека для приближённых защищённых вычислений на CKKS поверх Microsoft SEAL.

Стабильный API подключается через:

```cpp
#include "m2424/m2424.hpp"
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
./build/bin/demo_basic
```

## Приложения

- `demo_basic` — минимальный цикл шифрования и вычисления.
- `demo_secure_stats` — зашифрованные сумма и среднее.
- `bench_ckks` — измерение базовых CKKS-операций.

Остальные demos и benchmarks собираются по умолчанию и служат примерами API или измерительными утилитами.

## Тесты

По умолчанию CTest запускает быстрые детерминированные инвариантные тесты.

## Кандидаты bootstrap

Библиотека содержит каталог проектных кандидатов `bootstrapCandidates()`. Это не
готовые реализации bootstrap: кандидат становится применимым только после
полной оценки подтверждённых границ ошибок через `forecastBootstrapFeasibility`.

Первый измерительный фильтр кандидатов запускается отдельно и не реализует
bootstrap-алгоритмы:

```bash
./build/bin/bench_bootstrap_candidates --all
./build/bin/bench_bootstrap_candidates balanced_8192_s55
```

Последняя базовая калибровка: `docs/bootstrap_candidate_baseline.md`.

Для CoeffToSlot сейчас реализованы FFT- и BSGS-кандидаты. До их benchmark
сравнения ни один не считается оптимальным; критерии и подготовительные работы
зафиксированы в `docs/coeff_to_slot_comparison_plan.md`.

## Пример API

```cpp
#include "m2424/m2424.hpp"

auto adapter = m2424::SealAdapter::create(m2424::profiles::basic_ckks());
adapter.generateKeys(true, true);

auto encrypted = adapter.encrypt(adapter.encode({1.0, 2.0, 3.0}));
auto product = adapter.multiply(encrypted, encrypted);
auto reduced = adapter.relinearize(product);
auto squared = adapter.rescaleToNext(reduced);
auto decoded = adapter.decode(adapter.decrypt(squared));
```

## Документация

- `docs/architecture.md`
- `docs/bootstrap_candidate_baseline.md`
- `docs/coeff_to_slot_comparison_plan.md`

Публичные классы и методы документированы в заголовочных файлах через Doxygen. После установки Doxygen HTML-документацию можно сгенерировать так:

```bash
cmake -S . -B build
cmake --build build --target docs
```

Главная страница будет создана в `build/docs/html/index.html`.

Структура каталога сборки:

```text
build/
  bin/       # demos, benchmarks и тесты
  lib/       # библиотека m2424 и Microsoft SEAL
  docs/html/ # HTML-документация Doxygen
```
