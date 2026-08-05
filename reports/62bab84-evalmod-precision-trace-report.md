# 62bab84 — отчёт по EvalMod precision validation и differential trace

## Цель коммита

Коммит `62bab84` выполняет один узкий этап: перестраивает контракт проверки точности
EvalMod и добавляет per-node MPFR-versus-CKKS trace. Approximation generator, Remez,
SlotToCoeff, domain model и interval certificate в этом коммите намеренно не
переписывались.

Исходная проблема была в том, что backend validation считалась успешной при одном
лишь факте выполнения ciphertext DAG:

```cpp
result.passed = result.executionSucceeded;
```

Из-за этого тест мог одновременно требовать `passed=true` и
`matchesPolynomialReference=false`. Такой контракт скрывал численную проблему.

## Новый precision contract

В `EvalModProblem` добавлен `EvalModPrecisionBudget` с независимыми полями:

```cpp
double implementation;
double approximation;
```

Validator теперь независимо проверяет:

```text
implementationError = backend − MPFR compiled polynomial
approximationError  = MPFR polynomial − exact EvalMod target
totalMeasuredError  = backend − exact EvalMod target
```

Итоговый `passed` равен true только при одновременном выполнении четырёх условий:

```text
executionSucceeded
&& implementationError <= implementationBudget
&& approximationError <= approximationBudget
&& totalMeasuredError <= targetAbsoluteError
```

Failure classification стала явной:

- `implementation_budget_exceeded`;
- `approximation_budget_exceeded`;
- `total_budget_exceeded`.

`executionSucceeded` и `BackendMeasured` по-прежнему показывают, что схема физически
исполнилась, но больше не подменяют verdict точности.

## MPFR compiled-DAG interpreter

Для differential trace реализовано отдельное MPFR-исполнение каждого узла compiled
DAG. Оно не вызывает прямой polynomial Horner вместо графа, а повторяет структуру
узлов:

- Input;
- EncodeConstant;
- Add/AddPlain;
- MultiplyPlain/MultiplyCipher;
- Relinearize;
- Rescale;
- ModSwitch;
- AlignScale.

Semantic identity operations сохраняют MPFR value, арифметические узлы выполняются с
512-bit precision. Decimal constants импортируются напрямую в MPFR, без `stod` в
reference path.

## Differential ciphertext trace

При validation executor после каждого ciphertext node выполняет decrypt/decode и
сравнивает результат с соответствующим MPFR DAG node. Для каждого измеренного узла
сохраняются:

```text
node index
operation
chain index
actual CKKS scale
MPFR semantic value в худшем slot
decrypted CKKS value в худшем slot
maximum absolute error
error increase относительно предыдущего ciphertext node
```

Также сохраняется `firstImplementationBudgetExceedingNode`. Differential mode
предназначен только для validation с secret key; production executor без запроса
trace не выполняет decrypt.

## Эксперимент degree 9

Для точностной диагностики зафиксирован degree-9 `PeriodicSineBaseline`, профиль
N=32768, scale около `2^59` и пять симметричных входов в центральном допустимом
интервале:

```text
-0.12, -0.06, 0, 0.06, 0.12
```

Normalization gain остаётся неединичным. Budgets заданы явно:

```text
implementationBudget = 1e-10
approximationBudget   = 1e-10
targetAbsoluteError   = 1e-10
```

Локальные запуски показали implementation error порядка `1e-12` или ниже. Ни один
узел не стал `firstImplementationBudgetExceedingNode`; MPFR polynomial reference
достигается с запасом.

При этом approximation error degree-9 periodic baseline существенно превышает
`1e-10`. Поэтому окончательный verdict кандидата:

```text
executionSucceeded          = true
matchesPolynomialReference  = true
matchesEvalModTarget        = false
passed                      = false
failure                     = approximation_budget_exceeded
```

CTest для этого сценария зелёный не потому, что кандидат объявлен точным. Он зелёный,
потому что validator правильно обнаружил и классифицировал ожидаемый approximation
failure. Старое требование `!matchesPolynomialReference` удалено.

## Почему в этом коммите не исправлялась CKKS-операция

Первоначальная гипотеза состояла в том, что превышение `1e-10` возникает на encoding,
MultiplyCipher, rescale или metadata alignment. Differential trace на фиксированном
degree-9 эксперименте эту гипотезу не подтвердил: implementation path укладывается в
budget, а первым падает независимый approximation gate.

Поэтому не было оснований менять конкретную CKKS operation только ради получения
зелёного total-error теста. Это исказило бы диагноз. Следующий шаг относится к
target-driven approximation solver, но он сознательно не включён в текущий узкий
коммит.

На более широких значениях domain metadata alignment остаётся отдельным известным
риском и должен устраняться настоящим prime-aware scale planner. Текущий результат не
объявляет этот вопрос решённым; он только локализует причину провала фиксированного
precision experiment.

## Изменения тестового контракта

Backend test теперь проверяет одновременно:

- ciphertext execution действительно состоялось;
- implementation reference достигнут;
- implementation error не превышает `1e-10`;
- approximation и total gates не достигнуты degree-9 baseline;
- `passed=false`;
- failure reason равен `approximation_budget_exceeded`;
- отсутствует implementation-budget-breaking node;
- differential trace непустой;
- short-chain и least-squares rejection продолжают работать;
- minimax executions всё ещё используются для backend/calibration coverage, но не
  выдаются за достижение total target.

## Проверка

После изменений выполнено:

```text
cmake --build build -j4
ctest --test-dir build --output-on-failure
report_evalmod_synthesis JSON parse
report_evalmod_synthesis CSV column consistency
git diff --check
```

Результат локального прогона:

```text
12/12 tests passed
test_evalmod_backend:
  implementation error около 1e-12
  verdict approximation_budget_exceeded
  test contract PASS
```

Статус GitHub Actions этим локальным прогоном не подтверждается.

## Что остаётся дальше

Следующий отдельный milestone должен заниматься approximation path:

1. отказаться от PeriodicSineBaseline как accuracy candidate;
2. ввести target-driven degree search;
3. использовать устойчивый basis, предпочтительно Chebyshev;
4. проверять MPFR direct polynomial против MPFR compiled-DAG на большой сетке;
5. принимать кандидата только при certified approximation bound ниже budget.

Scale planner, exact coefficient representation и operation-isolated calibration
также остаются необходимыми, но их не следует смешивать с уже выполненной
диагностикой в одном коммите.
