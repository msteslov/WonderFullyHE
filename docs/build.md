# Сборка и проверка

Команды выполняются из корня репозитория. Пути к исполняемым файлам приведены
для одноконфигурационного генератора CMake: Unix Makefiles или Ninja.

## Зависимости

Основная конфигурация требует компилятора C++17, CMake 3.16+ и Microsoft SEAL
из подмодуля `extern/seal`. После клонирования необходимо выполнить:

```bash
git submodule update --init --recursive
```

Экспериментальная конфигурация дополнительно требует `pkg-config`, GMP с C++
интерфейсом (`gmpxx`) и MPFR. Используемые имена пакетов:

| Среда | Пакеты анализа |
|---|---|
| Ubuntu | `pkg-config`, `libgmp-dev`, `libmpfr-dev` |
| macOS, Homebrew | `pkg-config`, `gmp`, `mpfr` |

Doxygen необходим только для генерации справочника.

## Основная конфигурация

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build -j
ctest --test-dir build --output-on-failure
./build/bin/demo_basic
```

| Опция CMake | Значение по умолчанию | Назначение |
|---|---|---|
| `BUILD_TESTING` | `ON` | Тесты и диагностические интерфейсы тестовых сценариев |
| `M2424_BUILD_APPS` | `ON` | Демонстрации и измерительные приложения |
| `M2424_ENABLE_EVALMOD_ANALYSIS` | `OFF` | Анализ, сертифицирующая подготовка и дополнительные тесты |

Для сборки основной библиотеки без приложений и тестов:

```bash
cmake -S . -B build-core-only -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_TESTING=OFF -DM2424_BUILD_APPS=OFF \
  -DM2424_ENABLE_EVALMOD_ANALYSIS=OFF
cmake --build build-core-only -j
```

Значения опций сохраняются в CMake cache. При повторном использовании каталога
нужную конфигурацию следует задавать явно.

## Экспериментальная конфигурация

```bash
cmake -S . -B build-analysis -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_TESTING=ON -DM2424_ENABLE_EVALMOD_ANALYSIS=ON
cmake --build build-analysis -j
ctest --test-dir build-analysis --output-on-failure
```

Основная цель — `m2424::m2424`, дополнительная — `m2424::evalmod_analysis`.
GMP/MPFR требуются только второй цели.

Некоторые методы объявлены в публичных заголовках, но их подготовка реализована
в analysis-цели: `Bootstrapper::prepare`, `SlotToCoeffPlan::prepare`,
`EvalRoundPlusCoeffToSlot::prepareCertified`. Наличие объявления само по себе
не обеспечивает линковку без соответствующей библиотеки.

## Тесты

Список тестов определяется конфигурацией:

```bash
ctest --test-dir build -N
ctest --test-dir build-analysis -N
ctest --test-dir build-analysis -R '^test_evalround_' --output-on-failure
ctest --test-dir build-analysis -R '^test_bootstrap_end_to_end$' --output-on-failure
```

Состав и таймауты заданы в `tests/CMakeLists.txt`.
Для части analysis-тестов установлен таймаут 1800 секунд. Количество тестов
и длительность из архивных отчётов относятся к соответствующим ревизиям.

## Приложения

| Исполняемый файл | Назначение |
|---|---|
| `demo_basic`, `demo_secure_stats` | CKKS-операции и зашифрованная статистика |
| `demo_abft` | Контроль численной согласованности |
| `demo_profile_report`, `demo_security_report` | Параметры профилей и отчёт безопасности |
| `demo_noise_growth`, `demo_precision_profiles` | Диагностика точности |
| `bench_ckks`, `bench_chain_accuracy`, `bench_parallel_throughput` | Измерения операций, цепочек и пропускной способности |
| `bench_coeff_to_slot` | Подготовка и исполнение CoeffToSlot |
| `bench_bootstrap_candidates` | Начальная калибровка конфигураций bootstrap |
| `report_evalmod_synthesis`, `report_evalmod_feasibility` | Синтез и анализ EvalMod; требуется analysis-цель |
| `report_evalround_k64_synthesis` | Синтез полиномов и математическое планирование K=64; требуется analysis-цель |

Примеры запуска:

```bash
./build/bin/demo_profile_report
./build/bin/bench_bootstrap_candidates --all
./build-analysis/bin/report_evalmod_synthesis --csv
./build-analysis/bin/report_evalmod_feasibility
./build-analysis/bin/report_evalround_k64_synthesis
```

Отчёт синтеза K=64 не выполняет полную проверку backend. Соответствующий
регрессионный сценарий запускается отдельно:

```bash
ctest --test-dir build-analysis -R '^test_evalround_k64_synthesis$' --output-on-failure
```

Калибровка кандидатов не исполняет полный bootstrap и не выдаёт его сертификат.

## Справочник API

После установки Doxygen необходимо повторить конфигурацию:

```bash
cmake -S . -B build
cmake --build build --target docs
```

Входы — заголовки `include/m2424/` и страницы `docs/`.
Стартовая страница — оглавление документации. Архив доступен как отдельный
раздел; его страницы явно обозначены как исторические отчёты.

`Doxyfile` задаёт фиксированный путь `build/docs/html/index.html` относительно
корня репозитория, в том числе при вызове цели из другого каталога сборки.
Без установленного Doxygen цель `docs` недоступна.

Исполняемые файлы находятся в `<build-dir>/bin/`, библиотеки —
в `<build-dir>/lib/`. Многоконфигурационные генераторы могут добавлять
подкаталог конфигурации; сборка выполняется с `--config Release`, тесты — с `-C Release`.
