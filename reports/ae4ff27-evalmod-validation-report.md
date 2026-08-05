# ae4ff27 — отчёт по разделению estimated/certified EvalMod validation

## Контекст

Коммит `ae4ff27` продолжает reusable EvalMod executor из `5d23a75`. Целью было
устранить смешение grid-based оценки и доказанного error path, заменить глобальный
почти вакуумный approximation bound локальным interval certificate, выполнить
настоящий EvalMod DAG после CoeffToSlot и сформировать практически применимую модель
арифметической ошибки без ложного статуса rigorous.

## Главный обнаруженный дефект

До этого `intervalCertified=true` мог сосуществовать с `predictedBootstrapError`,
построенным из `diagnostic.approximationMaxError`. Grid diagnostic полезен для ranking,
но не является доказательством. Если бы arithmetic model позднее получила флаг
rigorous, numerical gate мог бы использовать недоказанную approximation error.

Контракт разделён на два пути:

- `estimatedBootstrapError` использует grid diagnostic и служит только для поиска и
  provisional ranking;
- `certifiedBootstrapError` является `optional` и появляется только при строгих
  approximation и arithmetic bounds;
- `propagationBounds.approximation` всегда получает сертифицированную границу, если
  certificate существует;
- `satisfiesNumericalTarget` проверяет только `certifiedBootstrapError`.

Legacy-поле `predictedBootstrapError` оставлено для совместимости отчётов, но больше
не используется для provisional selection или rigorous numerical gate.

## Локальный approximation certificate

Прежний certificate ограничивал полином глобальной суммой абсолютных коэффициентов.
Граница была конечной и безопасной, но параметр `subdivisions` практически не влиял
на real bound.

Новая real-сертификация работает так:

1. Каждый допустимый интервал `[k-rho,k+rho]` делится на `subdivisions` ячеек.
2. Центр каждой ячейки строится как MPFR interval с направленным округлением.
3. Каждый decimal coefficient импортируется парой `RNDD/RNDU`; это отдельно важно
   для отрицательных коэффициентов.
4. Полином вычисляется interval Horner. Умножение интервалов рассматривает четыре
   комбинации границ с outward rounding.
5. Для `p(x)-x+k` берётся строгий upper bound в центре.
6. На всю ячейку bound переносится добавлением `(sup|p'|+1)*halfCell` с `RNDU`.
7. `proved=true` выставляется только при конечных bounds и полном покрытии.

Для complex rectangle пока сохранён более грубый maximum-modulus majorant. Это
осознанный компромисс: real approximation bound стал полезным локальным
сертификатом, а complex path остаётся безопасным, но консервативным.

На текущем report-профиле periodic baseline имеет grid error порядка `3.4e-3`, а
сертифицированный real bound порядка `4e-2`. Это заметно плотнее прежней глобальной
границы порядка десятков единиц и при этом не подменяет доказательство измерением.

## End-to-end CoeffToSlot → EvalMod на N=16384

Identity circuit был удалён из интеграционной проверки. Новый тест выполняет:

```text
ModRaise
→ реальный четырёхуровневый CoeffToSlot
→ compiled linear EvalMod DAG
→ independent exact MPFR EvalMod target
```

DAG содержит явную normalization multiplication и rescale, потребляет последний
доступный уровень и применяется к обеим половинам `CoeffToSlotResult` через
`executeEvalModAfterCoeffToSlot`.

Область теста намеренно задаётся как `K=0`: source modulus выбран так, чтобы все
фактические CoeffToSlot coordinates после normalization находились внутри
`(-1/2,1/2)`. Поэтому линейный polynomial `p(x)=x` является настоящим EvalMod для
этого domain, а не identity passthrough: ciphertext физически проходит normalization,
уровень и scale меняются, результат независимо сравнивается с exact target.

Измеренная end-to-end EvalMod error в локальном прогоне была порядка `2e-13` для
обеих ciphertext-половин; CoeffToSlot error осталась порядка `1e-10`.

### Почему не был оставлен кубический кандидат на N=16384

Сначала проверялся нелинейный cubic DAG. Чтобы освободить достаточно уровней, был
опробован CoeffToSlot depth 2. Такая факторизация резко увеличивает prepared diagonal
state и оказалась непригодной по памяти для обычного тестового запуска. При штатном
depth 4 security-compatible modulus chain N=16384 оставляет только один уровень,
чего недостаточно для нелинейного polynomial DAG при scale около `2^59.5`.

Вместо скрытого ослабления параметров принято следующее разделение:

- N=16384 проверяет настоящий, но domain-limited `K=0` linear EvalMod end-to-end;
- N=32768 продолжает ciphertext matrix для nonlinear PeriodicSine и
  MultiIntervalMinimax degree 7/9/11.

Это ограничение связано с level/security budget, а не с ownership API.

## Scale handling на простых около 2^60

Тест N=16384 использует initial scale `2^59.5` и 60-bit primes. После rescale
естественный scale отличается от номинального примерно на половину бита. Старый
единый threshold `1e-6` ошибочно трактовал бы это как metadata correction failure.

Теперь разделены два бюджета:

- `maxMetadataScaleCorrectionLog2=1e-6` — только для фактического изменения metadata;
- `maxPlannedScaleDriftLog2=0.55` — заявленный допуск естественного prime-induced
  drift, который сохраняется в trace и не исправляется молча.

`AlignScale` по-прежнему не может менять metadata при расхождении больше `1e-6`.
Больший естественный drift допустим только для плановой проверки отдельной ветви;
слияние несовместимых ветвей всё равно fail-fast.

## Улучшения polynomial compiler

End-to-end работа обнаружила две лишние траты уровней:

- compiler заранее строил все baby powers даже для degree 1;
- коэффициент `1` умножался как обычный plaintext coefficient с последующим rescale.

Исправления:

- powers строятся только до `min(babyStep, degree)`;
- monomial coefficient `1` переиспользует готовую power node;
- unit normalization/denormalization не создают фиктивные операции;
- exact integer denormalization может выполняться plaintext multiplication без
  дополнительного rescale.

Эти оптимизации следуют семантике DAG и уменьшают level consumption, не меняя
полиномиальную функцию.

## Arithmetic error model и калибровка

Универсальная добавка `2^-targetPrecisionBits` заменена структурой
`EvalModArithmeticErrorModel` с отдельными параметрами для:

- encoding;
- addition;
- multiplication;
- relinearization;
- rescale;
- ModSwitch;
- metadata scale correction.

`calibrateEvalModArithmeticModel` принимает серию успешных backend measurements и
строит conservative empirical model по максимальной implementation error с safety
factor. Backend test калибрует модель на periodic и minimax executions с несколькими
шифрованиями.

Калиброванная модель намеренно имеет `calibrated=true`, `rigorous=false`. Измерение с
safety factor полезно для инженерного planner, но не доказывает worst-case CKKS noise
для всех keys, messages и modulus chains.

## Что было проверено

Локально выполнено:

```text
cmake --build build -j4
ctest --test-dir build --output-on-failure
report_evalmod_synthesis JSON parse
report_evalmod_synthesis CSV column consistency
git diff --check
```

Результат: 12 из 12 тестов прошли. Отдельно прошли:

- N=16384 ModRaise → CoeffToSlot → real K=0 EvalMod;
- N=32768 periodic backend validation;
- MultiIntervalMinimax degree 7/9/11 с разными babyStep и несколькими encryption
  trials;
- JSON/CSV machine-readable reports.

GitHub Actions status в рамках этого локального коммита не подтверждался; утверждение
о прохождении относится к локальной сборке и `ctest`.

## Оставшиеся ограничения

- Complex certificate всё ещё использует глобальный majorant.
- Empirical arithmetic calibration не является CKKS proof.
- Полный nonlinear N=16384 pipeline не помещается в текущий level/security budget.
- Remez остаётся MPFR candidate generator с discretized exchange, а не формально
  certified minimax solver.
- Полный SlotToCoeff после EvalMod ещё не включён в этот end-to-end тест.

Следующий содержательный milestone — либо nonlinear end-to-end профиль N=32768 с
SlotToCoeff, либо новый depth/scale design, который докажет выполнимость nonlinear
N=16384 без нарушения security budget.
