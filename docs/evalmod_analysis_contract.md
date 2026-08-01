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

Текущая реализация является детерминированным grid diagnostic, но не формальным
interval-arithmetic proof. Будущий Arb API получит отдельный тип interval certificate.

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

## Следующий эксперимент

1. Снять empirical распределение `I` и residual на реальном выходе ModRaise/CtS.
2. Построить аналитическую secret/error model и сравнить quantiles, не подменяя bound
   эмпирическим максимумом.
3. Добавить generator periodic/minimax candidates и interval-certified derivative.
4. После этого формировать единый CSV synthesis report и выбирать circuit/chain.
