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

`analyzeEvalModDomain` принимает уже явно построенную субгауссову модель состояния
ciphertext. Он не пытается вывести её только из `q`, `Q` и message bound. Для `K`
используется union bound по всем коэффициентам с заданной failure probability;
`rho` выводится из bounds на encoded message и существующую ошибку, заранее
консервативно поделённых на exact `q_src`. API специально не принимает `q_src` как
`double`: exact division должна происходить в high-precision слое вызывающей стороны.

## Oracle

`exactCoefficientOracle` выполняет CRT в GMP integer, затем exact centered reduction
modulo полного `q_src`. Только после этого выполняется высокоточное деление на
semantic output scale. Это исключает потерю младших битов большого `Q` через
промежуточный `double`.

## Аппроксимация

`certifyEvalModPolynomial` измеряет monomial polynomial на объединении интервалов
`[k-rho,k+rho]`, а не на сплошном интервале с разрывами. Он сообщает maximum error,
maximum derivative и ошибку на границе заданной комплексной окрестности.

Текущая реализация является детерминированным grid diagnostic, но не формальным
interval-arithmetic proof. Поле называется certificate на уровне API, однако перед
использованием как доказательства его необходимо подкрепить interval/Arb backend.

## Распространение ошибки и стоимость

`propagatedBootstrapError` реализует

```text
||T_StC|| (L_F e_CtS + e_approx + e_poly + e_scale + e_period)
  + e_StC + e_final.
```

`estimateEvalModCost` отдельно оценивает latency, working set и modulus bits по
скомпилированному circuit description и измеренной backend cost model. Это ещё не
optimizer modulus chain и не security estimator.

## Следующий эксперимент

1. Снять empirical распределение `I` и residual на реальном выходе ModRaise/CtS.
2. Построить аналитическую secret/error model и сравнить quantiles, не подменяя bound
   эмпирическим максимумом.
3. Добавить generator periodic/minimax candidates и interval-certified derivative.
4. После этого формировать единый CSV synthesis report и выбирать circuit/chain.
