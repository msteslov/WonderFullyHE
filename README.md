# WonderFullyHE

C++17-библиотека приближённых гомоморфных вычислений на CKKS поверх Microsoft
SEAL. Библиотека предоставляет операции над шифротекстами, линейные
преобразования, компоненты CKKS bootstrapping и контроль границ ошибок.

Базовый интерфейс — `m2424::SealAdapter`; общий заголовок — `m2424/m2424.hpp`.
Реализованы исполнители стадий bootstrap, подготовка планов, проверка контрактов,
интеграционные тесты и воспроизводимые бенчмарки.

Документация собрана в **[docs/](docs/README.md)**:

- [Сборка и проверка](docs/build.md).
- [Архитектура](docs/architecture.md).
- [Bootstrap: маршрут и контракты](docs/bootstrap.md).
- [Состояние реализации и ограничения](docs/status.md).

## Быстрый запуск

Требуются CMake 3.16+, компилятор C++17 и Git. Из корня репозитория:

```bash
git submodule update --init --recursive
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build -j
ctest --test-dir build --output-on-failure
./build/bin/demo_basic
```

Расширенная конфигурация анализа требует GMP/MPFR и включается отдельно.
Условия сборки и полный перечень команд приведены в документации.

## Бенчмарки и анализ на Python

В каталоге [`python/`](python/README.md) находится отдельный пакет для запуска
C++-бенчмарков, валидации конфигураций, сбора CSV/JSON, статистического анализа и
построения графиков. Для быстрой проверки после сборки проекта:

```bash
python3 -m venv python/.venv
python/.venv/bin/pip install -e 'python[dev]'
python/.venv/bin/wfhe-experiments run python/configs/smoke.json
```

Методика измерений: [docs/research_protocol.md](docs/research_protocol.md).
