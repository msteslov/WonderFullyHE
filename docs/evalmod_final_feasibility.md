# Финальный feasibility-результат nonlinear EvalMod

## Результат

Для исходного профиля и для явно ограниченного класса ниже получен конечный
отрицательный результат:

```text
status = CertificateClassInfeasible
global EvalMod impossibility = false
```

Это вариант B из финального ТЗ, а не `NoCertifiedPlanInSearchSpace`.
Машиночитаемый отчёт строится командой:

```bash
./build-analysis/bin/report_evalmod_feasibility
```

## Область доказательства

Результат относится ко всем plan, которые генерирует текущий explicit-alpha
stable-scale polynomial compiler, включая direct и возможную polynomial
Composite-композицию, при следующих ограничениях:

* стандартные SEAL `N = 1024..32768`;
* `tc128` и полный coefficient-modulus не выше таблицы SEAL;
* coefficient primes не более 60 bit;
* special prime того же scale-класса;
* `S_out <= qSource`;
* явная нормализация `alpha = 1/2` через `MultiplyPlain + Rescale` до
  polynomial DAG;
* deterministic coefficient-support/canonical-embedding certificate,
  реализованный в `certifyEvalModDagArithmetic()`.

Это **не** доказательство глобальной невозможности EvalMod. В частности, оно
не покрывает новый compiler, который поглощает `alpha` в коэффициенты, новый
no-rescale/non-stationary layout или backend с иной доказанной моделью
rounding.

## Исходный профиль

```text
N                          = 32768
K                          = 1
rho                        = 0.08
qSource                    = 2^60
CtS scale                  = 2^59
output scale               = 2^57
output gain                = 8
target final               = 1.0e-10
target normalized          = 1.25e-11
actual modulus bits        = 533
tc128 modulus limit        = 881
max same-scale data primes = 13
max levels                 = 12
max dense odd degree upper = 41
```

`max levels` здесь не выбранная вручную граница. Это
`floor((881 - 60) / 59) - 1 = 12` после резервирования actual 60-bit special
prime.

## Полный numerical budget

Таблица ниже получена для реального prepared degree-9 nonlinear plan. Все
`current` arithmetic-строки являются deterministic rigorous bounds. Колонка
`final` включает exact denormalization gain 8 и SlotToCoeff operator norm 1.
Нулевое «best» нигде не выдаётся за достигнутый rigorous bound: если более
сильного доказанного механизма нет, повторяется current bound. Для Rescale
отдельно указан минимальный обязательный operation bound внутри класса.

| Source | current normalized | best justified rigorous | fail log2 | count | propagated final |
|---|---:|---:|---:|---:|---:|
| Input/CtS | 0 | 0 | `-inf` | 1 | 0 |
| Approximation | 3.0519165e-2 | 3.0519165e-2 | `-inf` | 1 | 2.4415332e-1 |
| Coefficient quantization | 1.7788540e-8 | 1.7788540e-8 | `-inf` | 22 | 1.4230832e-7 |
| MultiplyPlain local error | 0 | 0 | `-inf` | 19 | 0 |
| Cipher×Cipher cross propagation | 5.9741232e-15 | 5.9741232e-15 | `-inf` | 6 | 4.7792986e-14 |
| Rescale | 2.1306879e-6 | 9.3135100e-10 operation floor | `-inf` | 14 | 1.7045503e-5 |
| Relinearization/KeySwitch | 2.1096577e-22 | 2.1096577e-22 | `-inf` | 6 | 1.6877262e-21 |
| ModSwitch | 0 | 0 | `-inf` | 9 | 0 |
| Scale representation | 0 | 0 | `-inf` | 0 | 0 |
| Period mismatch | 0 | 0 | deterministic | 1 | 0 |
| Denormalization local error | 0 | 0 | `-inf` | 1 | 0 |
| SlotToCoeff additive | 0 | 0 | `-inf` | 1 | 0 |
| Final additive | 0 | 0 | `-inf` | 1 | 0 |

Arithmetic total representative plan:

```text
2.148476447249207e-6 normalized
```

Dominant source and требуемый выигрыш:

```text
dominant source             = Rescale
current propagated bound    = 2.1306879011734705e-6
normalized budget           = 1.2499999999999989e-11
required improvement factor = 1.7045503209387779e5
```

## Строгий stop для исходного профиля

Для size-2 ciphertext divide-and-round оставляет два residual polynomial с
coefficient support `1/2`. Для ternary secret `|s_i| <= 1` текущий certificate
даёт

```text
B_RS = (N/2) * (1 + N) / Delta.
```

На actual первом normalization Rescale это
`9.3135099619435707e-10`, уже в 74.508 раза больше всего normalized target.

Нельзя просто утверждать, что local error обязан без изменения дойти до
выхода. Поэтому proof дополнительно использует необходимую чувствительность.
Если polynomial `p` аппроксимирует identity на центральном интервале
`[-rho,rho]` с ошибкой не выше `epsilon`, то

```text
sup |p'| >= |p(rho)-p(-rho)|/(2 rho) >= 1 - epsilon/rho.
```

Absolute DAG propagation certificate не может иметь input-error amplification
меньше этой величины: сумма absolute path sensitivities всегда не меньше
`sup |p'|`. При target `epsilon <= 1.25e-11` фактор равен как минимум
`0.99999999984375`. Следовательно, обязательный normalization Rescale не может
быть ослаблен downstream coefficients до target внутри указанного класса.

Этот stop не зависит от degree, baby-step или latency. Resource-only enumeration
текущего dense odd Paterson–Stockmeyer compiler даёт верхнюю степень 41 при
security-valid depth 12 (legacy alignment делает это именно верхней, а не
оптимистичной executable границей). Поэтому все direct Minimax/Chebyshev degrees
до этой границы исчерпаны строгим pruning
до дорогостоящего Remez search; arbitrary `maximumDegree=15` в доказательстве
не используется.

Для полноты отчёт сохраняет bounded diagnostics, но не использует их как stop:

| family | diagnostic degrees | best diagnostic approximation | resource degree upper | final status |
|---|---:|---:|---:|---|
| Direct Minimax | 7–15 | 2.1816321340093734e-2 | 41 | pruned by mandatory normalization Rescale |
| Direct Chebyshev | 7–15 | 2.1816321340093734e-2 | 41 | pruned by mandatory normalization Rescale |
| Composite arcsin/sine | not generated | not claimed | 41 | pruned by feasibility before implementation |

Таким образом, degree 15 — только последняя diagnostic точка. Конечный результат
покрывает всю resource-bound область до 41 единым degree-independent proof.

## Почему probabilistic arithmetic не реализован

Фактический `CBD(eta=21)` удовлетворяет MGF-bound

```text
E exp(tX) <= exp(10.5 t^2 / 2).
```

Даже если оптимистично разрешить non-adaptive first-use decomposition weights
и отдать всем relinearization union budget `2^-128`, оценка KeySwitch при его
реальном pre-rescale scale составляет около `1.17e-26`, а вместе с его ModDown
около `2.14e-26`. Это подтверждает, что KeySwitch не является bottleneck.

Для главного источника — Rescale — rounding residual является функцией
ciphertext coefficients. Независимость и zero mean относительно адаптивного
ciphertext не доказаны. Поэтому Hoeffding/subgaussian concentration для него
применять нельзя. Probabilistic KeySwitch implementation улучшил бы уже
несущественную строку и не приблизил бы план к target, поэтому feasibility gate
его отклоняет до кодинга.

## Composite feasibility

Для

```text
arcsin(sin(2*pi*z))/(2*pi)
```

при `rho=0.08` branch margin до ближайшего rounding discontinuity равен `0.42`.
То есть математическая ветвь корректна и approximation-side Composite имеет
смысл исследовать. Но в текущем explicit-alpha compiler Composite проходит тот
же normalization Rescale до sine DAG и тот же sensitivity stop. Double-angle
меняет degree/depth, но не устраняет этот первый gate. Поэтому ТЗ запрещает
тратить код на Composite implementation, который заранее не способен пройти
arithmetic budget.

## Security-valid profile exhaustion

Для каждого `N` scale выбран максимально благоприятно: не больше 60 bit и не
больше трети `tc128` modulus budget, чтобы остались два data prime и same-scale
special prime.

| N | tc128 bits | max stable scale bits | max data primes | max levels | one-RS floor | result |
|---:|---:|---:|---:|---:|---:|---|
| 1024 | 27 | 9 | 2 | 1 | 1.025e3 | arithmetic stop |
| 2048 | 54 | 18 | 2 | 1 | 8.00390625 | arithmetic stop |
| 4096 | 109 | 36 | 2 | 1 | 1.2210011e-4 | arithmetic stop |
| 8192 | 218 | 60 | 2 | 1 | 2.9107383e-11 | approximation stop |
| 16384 | 438 | 60 | 6 | 5 | 1.1642243e-10 | arithmetic stop |
| 32768 | 881 | 60 | 13 | 12 | 4.6567550e-10 | arithmetic stop |

У `N=8192` operation floor меньше максимально возможного normalized budget
`1e-10`, но доступен только один level. Даже если бесплатно поглотить часть
линейной работы, это не даёт polynomial degree выше 2. Для любого quadratic,
который имеет error не выше `epsilon` в точках `-1,0,1`, Lagrange interpolation
в точке `rho` даёт

```text
epsilon >= rho / (2 + rho - rho^2)
        = 0.03858024691358025.
```

Это на восемь порядков выше target. Таким образом, все standard tc128 profiles
в области доказательства имеют явный arithmetic или approximation stop.

## Оставшиеся классы

Результат оставляет открытыми и явно рекомендует следующие направления:

1. Compiler, который поглощает `alpha` в polynomial coefficients и строит новый
   error-layout proof без обязательного normalization Rescale.
2. No-rescale/non-stationary exact-scale layout с отдельным modulus-headroom
   доказательством.
3. Backend, где stochastic rounding residual действительно независим от
   adaptive ciphertext и это доказано.
4. Неполиномиальная bootstrapping-конструкция или другой HE backend.

Отдельная regression `test_evalmod_feasibility` проверяет exact primes/security,
полную таблицу sources, sensitivity/profile stops, Composite stop и отсутствие
глобального claim.
