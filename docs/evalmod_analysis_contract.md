# EvalMod: математический контракт анализа

Этот milestone намеренно не содержит homomorphic `EvalMod::apply()`. Он отделяет
семантику, вероятностную параметризацию, численную аппроксимацию и backend cost.

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
`double`; NaN и infinity отклоняются. Сейчас реализован
только monomial evaluator; Chebyshev и Composite являются явными, но пока отклоняемыми
вариантами basis, а не неявно смешанными схемами.

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
минимальный scale headroom lazy evaluation. Полный scale/prime schedule всё ещё должен
заменить эту временную нижнюю оценку. Cost строится по
скомпилированному circuit description и измеренной backend cost model. Это ещё не
optimizer modulus chain и не security estimator. Исследовательские компоненты собраны
в optional target `m2424::evalmod_analysis`; стабильная библиотека `m2424::m2424` не
зависит от GMP/MPFR. API находится в `m2424/experimental/evalmod_analysis/` и не входит
в стабильный umbrella header.

## Approximation synthesis milestone

`synthesizeEvalMod` связывает problem, rigorous domain estimate, parameter search,
polynomial DAG, scale schedule, MPFR diagnostic, propagation и backend
cost. Поиск перебирает degree и Paterson–Stockmeyer baby step для odd multi-interval
MPFR Remez. Дополнительно строятся periodic-sine baseline, исключённый из выбора
least-squares diagnostic prototype и inverse-sine composition.
Remez-кандидат допускается к сравнению только при упорядоченных extrema, чередовании
знаков ошибки, сбалансированных амплитудах, устойчивом exchange set и стабилизации
dense-grid maximum. Исключение одного кандидата превращается в
`ApproximationNotConverged` и не прерывает synthesis.

Operation counts, depth, level consumption и peak liveness выводятся из узлов DAG.
`ModSwitch`, `AlignScale` и `AddPlain` представлены отдельными узлами и входят в
стоимостную модель, а не вставляются backend-исполнителем скрыто.
Нормировка `S_CtS/q_src` и денормировка `q_src/S_out` хранятся exact rational и входят
в scale propagation. Arithmetic error оценивается по операциям DAG; неизвестное
значение никогда не трактуется как ноль. Для каждого узла сохраняются value/error,
relative-scale и noise bounds. `EvalModArithmeticErrorModel` задаёт разные constants
для encoding, add, multiply, relinearize, rescale, ModSwitch и metadata correction;
серия backend measurements может быть сведена в conservative calibrated model через
`calibrateEvalModArithmeticModel`. Такая калибровка остаётся empirical, поэтому
`rigorous=false` до отдельного вывода CKKS/SEAL bounds.

Кандидаты отдельно хранят `circuitValid`, `backendRunnable`, `backendMeasured`,
`approximationCertified`, `arithmeticErrorCertified` и `rigorouslyValidated`.
`executable` означает возможность экспериментального запуска, а не математическое
доказательство. `estimatedBootstrapError` использует grid diagnostic только для
эвристического ranking, а `certifiedBootstrapError` появляется лишь при строгих
approximation и arithmetic bounds. JSON/CSV сохраняют пути раздельно.

## Backend validation

`executeEvalModDagPlaintext` является smoke-проверкой семантики compiled DAG.
Независимый reference вычисляет polynomial напрямую Horner-схемой в MPFR с exact
rational normalization/denormalization. Production API `executeEvalModCircuit`
принимает существующие adapter/ciphertext, буквально исполняет DAG и возвращает trace
chain/scale каждого узла. Metadata alignment ограничен явным малым бюджетом
`maxMetadataScaleCorrectionLog2`; отдельный `maxPlannedScaleDriftLog2` допускает
заранее заявленное естественное отклонение после rescale простыми около `2^60`, не
разрешая подменять его metadata correction. Обёртка
`executeEvalModAfterCoeffToSlot` применяет тот же circuit к обеим половинам настоящего
`CoeffToSlotResult`.

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
