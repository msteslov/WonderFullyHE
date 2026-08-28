# EvalMod: математический контракт анализа

Компонент разделяет experimental synthesis/measurement и prepared execution.
`prepareEvalMod`, `preflightEvalMod` и `applyEvalMod` образуют production-shaped API:
последняя функция не синтезирует polynomial, не кодирует coefficients, не decrypt-ит
и не выбирает fallback. План получает `rigorous=true` только после всех реализованных
certificate gates; отрицательный план сохраняет структурированный status.

Текущий nonlinear `N=32768, K=1, rho=0.08` regression ещё не является полным
production certificate: degree-9 periodic approximation не проходит approximation
budget, а deterministic-support arithmetic bound порядка `2.2e-6` превышает
normalized target `1.25e-11`. Это отрицательный результат для конкретного candidate,
не утверждение глобальной невозможности. Нужны proven probabilistic norm bounds и
target-driven поиск более сильной approximation/composition.

## Нормировка

После ModRaise и CoeffToSlot анализируется

```text
x = (m_tilde + e + q_src I) / S_CtS
z = (S_CtS / q_src) x = I + (m_tilde + e) / q_src
```

Контракт домена хранит две независимые величины:

- `integerBound = K`, определяющую число допустимых интервалов;
- `normalizedResidualBound = rho < 1/2`, определяющую расстояние до разрывов.

`estimateEvalModDomain` принимает уже явно построенную tail-модель состояния
ciphertext. `integerNoiseSubgaussianSigma` — доказанный subgaussian parameter, а не
обычная standard deviation. Его необходимо вывести из распределений secret и
encryption error, key-switching noise и фактического состояния ciphertext. Тип модели
фиксируется как deterministic, Gaussian или subgaussian. Функция не пытается вывести
его только из `q`, `Q` и message bound. Для `K` используется MPFR union bound с
направленным вверх округлением по всем коэффициентам и заданной failure probability;
`rho` включает bounds на encoded message, существующую encoding error, нормированную
ошибку CoeffToSlot и additive normalization error, заранее
консервативно поделённых на exact `q_src`. API специально не принимает `q_src` как
`double`: exact division должна происходить в high-precision слое вызывающей стороны.
Relative period mismatch умножается внутри анализа на `K + residual`, поэтому caller
не может забыть зависимость period error от найденного integer bound. Derived
discontinuity margin не хранится и вычисляется с округлением вниз.

## Oracle

`exactCoefficientOracle` выполняет CRT в GMP integer, затем exact centered reduction
modulo полного `q_src`. Только после этого выполняется высокоточное деление на
semantic output scale. Это исключает потерю младших битов большого `Q` через
промежуточный `double`. Scale хранится как exact rational; `fromBinaryDouble`
импортирует точные mantissa/exponent фактического SEAL scale. Результат сохраняет
exact numerator/scale, rounded decimal и MPFR rounding-error bound.

## Аппроксимация

`diagnoseEvalModPolynomialOnGrid` измеряет polynomial на объединении интервалов
`[k-rho,k+rho]`, а не на сплошном интервале с разрывами. Он сообщает maximum error,
real/complex derivative maxima и ошибку на всей границе комплексной окрестности.
Коэффициенты, radius и координаты grid передаются/строятся в MPFR без промежуточного
`double`; NaN и infinity отклоняются. Исполняемый evaluator остаётся monomial.
Chebyshev synthesis функционален: MPFR Remez строит coefficients в базисе `T_k`,
после чего exact-rational recurrence конвертирует их в monomial coefficients без
ошибки basis conversion. Composite basis пока не является certified family.

Grid diagnostic остаётся измерением. Отдельный certificate разбивает каждый real
interval на `subdivisions` ячеек, вычисляет center через outward MPFR intervals,
исполняет interval Horner и добавляет Lipschitz half-cell gap с `RNDU`. Коэффициенты
импортируются парой `RNDD/RNDU`, поэтому отрицательные значения не занижаются.
Complex rectangle пока использует более грубый maximum-modulus majorant.
`proved=true` выставляется только для конечного полного покрытия.

## Распространение ошибки и стоимость

`propagatedBootstrapError` реализует

```text
||T_StC|| G_out (L_F e_CtS + e_approx + e_poly + e_scale + e_period)
  + e_StC + e_final.
```

где `G_out = q_src / S_out` — консервативная денормировка.

`estimateEvalModCost` отдельно оценивает latency, working set и нижнюю оценку data
modulus bits по maximum из multiplicative depth и level consumption; это учитывает
минимальный scale headroom lazy evaluation. Cost строится по
скомпилированному circuit description и измеренной backend cost model. Это ещё не
optimizer modulus chain и не security estimator. Исследовательские компоненты собраны
в optional target `m2424::evalmod_analysis`; стабильная библиотека `m2424::m2424` не
зависит от GMP/MPFR. API находится в `m2424/experimental/evalmod_analysis/` и не входит
в стабильный umbrella header.

### Exact prime/scale schedule

`SealAdapter` публикует фактические `dataModulusValues`,
`specialKeyModulusValue` и активные `coeffModulusValues(cipher)`. Certified planner
получает primes активной базы на входе EvalMod, а не только запрошенные
битовые размеры из профиля.

`buildExactEvalModScaleSchedule` хранит ideal rational scale, битовый образ actual
binary64 scale и exact rational representation этого binary64 для каждого DAG-узла.
Для `Rescale` используется конкретный сбрасываемый prime:

```text
Delta_runtime,out = fl64(Delta_runtime,in / double(p_dropped)).
```

Certified compiler выравнивает ветви настоящим prepared `MultiplyPlain(1)` и при
необходимости `Rescale`; nonlinear certified DAG не содержит `AlignScale`. Legacy
circuit сохраняет metadata-only узлы только для diagnostics. `ModSwitch` сохраняет
scale и не добавляет CKKS arithmetic error,
но `certifyEvalModModSwitchHeadroom` доказывает до его исполнения:

```text
2 * ceil((M + E) * Delta) < Q_next.
```

Нестрогий value/error bound может дать diagnostic headroom, но не может
перевести `headroomCertificate.rigorous` в `true`.

### Prepared exact constants

`prepareEvalModConstants` разбирает decimal coefficient как exact rational, умножает
его на exact binary64 scale, который будет записан в SEAL plaintext metadata, и
округляет scaled value в integer без `std::stod`. Integer раскладывается по
активным RNS primes и заранее записывается в NTT scalar plaintext нужного
уровня. Для каждой константы хранятся rounded integer и exact encoding-error
certificate.

`executeEvalModCircuitDiagnostic` оставляет parsing/encoding только для исследований.
`executePreparedEvalMod` требует `PreparedEvalModConstants`, проверяет modulus base,
полноту и привязку constants к DAG-узлам и не вызывает `std::stod` или CKKS encoder.
`executeEvalModAfterCoeffToSlot` также принимает только prepared constants.

### Certified arithmetic bounds

`certifyEvalModDagArithmetic` не использует empirical constants. Для Add применяется
`M_a+M_b`, `E_a+E_b`; для multiplication —
`M_a E_b + M_b E_a + E_a E_b`; coefficient quantization приходит из prepared
constant certificate. `BoundKind::Unknown` никогда не заменяется нулём.

Текущая vendored SEAL сборка использует ternary secret support `[-1,1]` и
centered-binomial evaluation-key noise support `[-21,21]`. Из конечных support,
RNS decomposition и special-prime ModDown выводится deterministic key-switch bound.
Для rescale coefficientwise rounding residual имеет `|r_i|<=1/2`; переход к slot
norm использует явную canonical-embedding triangle inequality.
`SealAdapter::analyzeCkksRescale` read-only извлекает actual dropped-prime limb,
делает inverse NTT и строит residual bound. Differential measurements используются
только для проверки `observed <= certified`, не для построения certificate.

`prepareEvalModSearch` перекомпилирует каждый реализованный candidate под exact
runtime context и возвращает либо prepared plan, либо ограниченный отрицательный
результат. `NoCertifiedPlanInSearchSpace` содержит families, degree/depth range,
лучшие approximation/arithmetic bounds, первый failing gate и всегда явно сообщает
`globalImpossibilityProved=false`. Resource infeasibility и отсутствие rigorous
operation bound имеют отдельные статусы.

## Approximation synthesis milestone

`synthesizeEvalMod` связывает problem, rigorous domain estimate, parameter search,
polynomial DAG, scale schedule, MPFR diagnostic, propagation и backend
cost. Поиск перебирает degree и Paterson–Stockmeyer baby step для odd multi-interval
MPFR Remez в monomial и Chebyshev bases. Дополнительно строятся periodic-sine baseline, исключённый из выбора
least-squares diagnostic prototype и inverse-sine composition.
Remez-кандидат допускается к сравнению только при упорядоченных extrema, чередовании
знаков ошибки, сбалансированных амплитудах, устойчивом exchange set и стабилизации
dense-grid maximum. Исключение одного кандидата превращается в
`ApproximationNotConverged` и не прерывает synthesis.

Operation counts, depth, level consumption и peak liveness выводятся из узлов DAG.
`ModSwitch`, `AlignScale` и `AddPlain` представлены отдельными узлами и входят в
стоимостную модель, а не вставляются backend-исполнителем скрыто.
Нормировка `S_CtS/q_src` и денормировка `q_src/S_out` хранятся exact rational и входят
в scale propagation. Legacy arithmetic estimate остаётся только diagnostic. Для
каждого узла сохраняются value/error, relative-scale и noise estimates.
`EvalModArithmeticErrorModel` задаёт разные empirical constants
для encoding, add, multiply, relinearize, rescale, ModSwitch и metadata correction;
серия backend measurements может быть сведена в conservative calibrated model через
`calibrateEvalModArithmeticModel`. Такая калибровка всегда остаётся
`rigorous=false`; даже исторически выставленный caller-ом flag не способен создать
certificate. Строгий результат приходит только из `EvalModArithmeticCertificate`.

Кандидаты отдельно хранят `circuitValid`, `backendRunnable`, `backendMeasured`,
`approximationCertified`, `arithmeticErrorCertified` и `rigorouslyValidated`.
`executable` означает возможность экспериментального запуска, а не математическое
доказательство. `estimatedBootstrapError` использует grid diagnostic только для
эвристического ranking, а `certifiedBootstrapError` появляется лишь при строгих
approximation и arithmetic bounds. JSON/CSV сохраняют пути раздельно.

## Backend validation

`executeEvalModDagPlaintext` является smoke-проверкой семантики compiled DAG.
Независимый reference вычисляет polynomial напрямую Horner-схемой в MPFR с exact
rational normalization/denormalization. Diagnostic executor принимает существующие
adapter/ciphertext, буквально исполняет DAG и возвращает trace chain/scale каждого
узла. Legacy metadata alignment ограничен явным малым бюджетом
`maxMetadataScaleCorrectionLog2`; отдельный `maxPlannedScaleDriftLog2` допускает
заранее заявленное естественное отклонение после rescale простыми около `2^60`, не
разрешая подменять его metadata correction. Обёртка
`executeEvalModAfterCoeffToSlot` применяет тот же circuit к обеим половинам настоящего
`CoeffToSlotResult`.

Prepared-plan regression проверяет context fingerprint, exact active primes,
chain level, exact binary64 input scale, ciphertext size и наличие evaluation keys.
`applyEvalMod` также исполняется на public evaluator, в котором secret key отсутствует.

Validator только создаёт тестовый ciphertext, вызывает production executor и отдельно
сообщает implementation, approximation и total measured errors, совпадение с
polynomial/точным EvalMod target и покрытие prediction. Успешное выполнение даёт
`BackendMeasured`; `BackendValidated` требует также строгой arithmetic bound.

`EvalModPrecisionBudget` задаёт независимые implementation и approximation budgets.
Legacy `passed` теперь истинно только если одновременно выполнены execution,
implementation, approximation и total gates; оно больше не является синонимом
`executionSucceeded`. При differential validation executor расшифровывает каждый
ciphertext node и сравнивает его с MPFR-исполнением того же compiled DAG. Trace хранит
operation, chain index, actual scale, MPFR/CKKS values, absolute error и её приращение.
На фиксированном degree-9 periodic circuit implementation error находится около
`1e-12` или ниже и проходит `1e-10`, тогда как approximation gate ожидаемо падает:
baseline не выдаётся за точную EvalMod approximation.
Backend-тест использует неединичные normalization/denormalization gains, boundary и
random точки всех интервалов, scale около `2^59`, несколько encryption trials и
PeriodicSine/MultiIntervalMinimax 7/9/11 с разными baby-step при `N=32768`.
Least-squares prototype остаётся запрещённым независимо от его grid error.

На `N=16384` end-to-end test исполняет настоящий линейный `K=0` EvalMod DAG после
четырёхуровневого CoeffToSlot и сравнивает обе половины с exact MPFR target. При scale
около `2^59.5` в доступный security budget не помещается одновременно этот CoeffToSlot
и нелинейный minimax DAG; полный нелинейный профиль остаётся тестом `N=32768`.
