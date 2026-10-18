# Структура репозитория

## Слои и направление зависимостей

```text
public API (include/m2424)
    ↓
planning / certification (src/planning; optional exact analysis uses GMP/MPFR)
    ↓
prepared immutable plan (public value/PIMPL types)
    ↓
runtime execution (src/ckks, src/math)
    ↓
SealAdapter backend boundary (src/ckks/seal_adapter*.cpp)
    ↓
Microsoft SEAL
```

Зависимость направлена вниз. Runtime не вызывает synthesis/search и не использует
reference/oracle для принятия решения. Наблюдатели runtime получают диагностику,
но не участвуют в арифметическом расписании.

## Ответственность директорий

| Путь | Ответственность |
|---|---|
| `include/m2424/` | Стабильные типы, контракты, планы и runtime API. `m2424.hpp` — стабильный umbrella header. |
| `include/m2424/experimental/evalmod_analysis/` | Явно opt-in интерфейсы точного анализа, synthesis/search, reference/oracle и cost/feasibility models. |
| `src/core/` | Контракты, общие типы, каталоги профилей/кандидатов и внутренние структуры данных. |
| `src/math/` | Математические преобразования CtS/StC и plaintext/reference реализация. Файлы exact/reference, включённые только в analysis target, не являются runtime oracle. |
| `src/planning/` | Построение и сертификация подготовленных планов. `evalround_execution_compiler.cpp` компилирует execution DAG; он не исполняет ciphertext. |
| `src/ckks/` | Исполнение подготовленных планов и адаптация к SEAL. `evalround_execution.cpp` содержит runtime DAG executor и liveness. |
| `tests/` | Тесты и test-only fixtures/oracles. Native SEAL-доступ тестов сосредоточен в `seal_adapter_test_support.cpp`. |
| `apps/` | Demos, benchmarks и reports, сгруппированные списками CMake без изменения target names. |

Внутренние заголовки имеют суффикс `_internal.hpp`, находятся в `src/` и не
устанавливаются. `src/ckks/seal_adapter_internal.hpp` хранит PIMPL/native SEAL
состояние, общее для implementation units; это не публичный API.

## Public / experimental / internal

- Stable public: заголовки непосредственно в `include/m2424/`. Они не включают
  GMP/MPFR. Reference-функции canonical embedding публичны, но документированы
  как plaintext `O(N^2)` reference, не production runtime.
- Experimental public: всё под `include/m2424/experimental/evalmod_analysis/`.
  Эти интерфейсы доступны только при `M2424_ENABLE_EVALMOD_ANALYSIS=ON` и могут
  зависеть от GMP/MPFR.
- Internal: все `_internal.hpp` в `src/`, а также test fixtures/oracles в `tests/`.
  Они не должны включаться клиентским кодом.

Часть методов `prepare` стабильных планов реализована в optional analysis target.
Это существующая link-time граница, сохранённая ради совместимости; клиент,
вызывающий эти методы, должен линковать `m2424::evalmod_analysis`.

## Bootstrap paths

Production execution имеет один маршрут:

```text
Cipher input
  → SealAdapter::modRaiseToTop / modRaiseSparse
  → prepared certified CoeffToSlot HP+LP
  → executeEvalRound on the certified LP DAG
  → executeCertifiedArithmetic for HP - gamma * I_hat
  → prepared certified SlotToCoeff
  → Cipher output
```

`Bootstrapper::prepare` строит этот маршрут в `src/planning/bootstrap_planner.cpp`.
`Bootstrapper::preflight` проверяет контекст, вход, сертификаты и ключи до
исполнения. `Bootstrapper::apply` в `src/ckks/bootstrap.cpp` только исполняет
сохранённые стадии. Legacy sparse restore доступен лишь при test diagnostics;
production restore выполняется внутри сертифицированного первого фактора CtS.

## CtS и StC

Обе стороны используют порядок `plan → requirements → prepare/certificate →
preflight → apply → metrics` и неизменяемые prepared plans.

CtS намеренно сложнее: он принимает raised ciphertext, строит две ветви HP/LP,
совместно использует first-factor sparse restoration и выпускает domain evidence
для EvalRound. StC принимает две уже очищенные половины и объединяет их в один
ciphertext. Поэтому CtS certification хранится рядом с dual-branch preparation,
а StC certificate — в `PreparedSlotToCoeffPlan`; это не два execution path.

## EvalRound

| Сущность | Реализация |
|---|---|
| Математическая модель/polynomial representation | `include/m2424/evalround.hpp`, `src/planning/evalround_polynomial_internal.hpp` |
| Synthesis/search | `src/planning/evalround_synthesis.cpp`, experimental headers |
| Компиляция сертифицированного DAG | `src/planning/evalround_execution_compiler.cpp` |
| Runtime исполнение DAG | `src/ckks/evalround_execution.cpp` |

`executeCertifiedArithmetic` сохраняет порядок узлов и освобождает значения по
последнему использованию. Planning builder находится в отдельном
`certified_arithmetic_planning_internal.hpp`; runtime-декларация — в
`certified_arithmetic_runtime_internal.hpp`.

## CMake

Production target `m2424::m2424` состоит из явно названных групп core,
backend-SEAL, math и planning и зависит только от SEAL. Optional target
`m2424::evalmod_analysis` состоит из analysis math/planning, зависит от
`m2424::m2424`, GMP и MPFR. Test-only native helpers собраны в
`m2424_test_support` и не входят в production-библиотеку.

## Будущие изменения

- Изменения SEAL primitives, сериализации и key handling размещаются в
  `seal_adapter*.cpp`, не в planners.
- Новые backend implementations должны реализовывать отдельную backend boundary;
  runtime plan/execution не должен включать vendor headers.
- Performance/memory изменения prepared execution размещаются в runtime units
  после отдельного измеримого baseline PR.
- CPU/GPU выбор должен появиться ниже runtime executor и выше конкретного backend.
- Boolean DAG является отдельным верхним execution layer и не добавляется в
  arithmetic/certification types EvalRound.
