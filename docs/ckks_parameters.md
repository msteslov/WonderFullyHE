# CKKS-параметры

WonderFullyHE использует Microsoft SEAL как CKKS-движок. Готовые параметры находятся в `m2424::profiles`, а подбор параметров для конкретного вычисления выполняет `plan_ckks_parameters`.

## Профили

Актуальная таблица профилей формируется командой:

```bash
./build/demo_profile_report
```

| Профиль | N | Слоты | Биты коэффициентных модулей | Суммарный размер, бит | Масштаб | log2(scale) | Оценка глубины умножений |
|---|---:|---:|---|---:|---:|---:|---:|
| `fast_demo_ckks` | 4096 | 2048 | `40-30-30` | 100 | `1.07374e+09` | 30 | 1 |
| `basic_ckks` | 8192 | 4096 | `60-40-40-60` | 200 | `1.09951e+12` | 40 | 2 |
| `balanced_ckks` | 8192 | 4096 | `50-40-40-40-48` | 218 | `1.09951e+12` | 40 | 3 |
| `depth_ckks` | 16384 | 8192 | `60-40-40-40-40-60` | 280 | `1.09951e+12` | 40 | 4 |
| `high_precision_ckks` | 16384 | 8192 | `60-50-50-50-60` | 270 | `1.1259e+15` | 50 | 3 |
| `boot_ckks` | 16384 | 8192 | `60-40-40-40-40-40-40-40-60` | 400 | `1.09951e+12` | 40 | 7 |
| `boot_deep_ckks` | 32768 | 16384 | `60-40-40-40-40-40-40-40-40-40-40-40-40-40-40-40-40-40-40-40-60` | 880 | `1.09951e+12` | 40 | 19 |
| `precision_boot_deep_ckks` | 32768 | 16384 | `60-50-50-50-50-50-50-50-50-50-50-50-50-50-50-50-60` | 870 | `1.1259e+15` | 50 | 15 |
| `precision_boot_ultra_ckks_59` | 32768 | 16384 | `60-60-60-60-60-60-60-60-60-60-60-60-60-60` | 840 | `5.76461e+17` | 59 | 12 |

## Назначение профилей

| Задача | Профиль |
|---|---|
| Быстрая локальная проверка сборки и простых операций | `fast_demo_ckks` |
| Основные demo, ABFT и базовые benchmark-и | `basic_ckks` |
| Дополнительный уровень умножения без перехода на `N = 16384` | `balanced_ckks` |
| Анализ глубины leveled-вычислений | `depth_ckks` |
| Более высокая точность обычных CKKS-цепочек | `high_precision_ckks` |
| Bootstrap-диагностика с 40-битным рабочим масштабом | `boot_ckks`, `boot_deep_ckks` |
| Bootstrap-диагностика с повышенным масштабом | `precision_boot_deep_ckks`, `precision_boot_ultra_ckks_59` |

## Масштаб и уровни

В основных профилях используется масштаб:

```text
scale = 2^40 ~= 1.0995 * 10^12
```

После умножения двух ciphertext масштаб примерно становится `2^80`. Операция `rescale_to_next` делит результат на ближайший рабочий модуль и возвращает масштаб к диапазону около `2^40`, расходуя один уровень modulus chain.

`fast_demo_ckks` использует `2^30`, потому что его рабочие модули 30-битные. `high_precision_ckks` использует `2^50`, что снижает ошибку округления, но увеличивает стоимость операций и требует большего `N`.

Длина chain задаёт доступную глубину вычисления. При фиксированном масштабе и одинаковых рабочих модулях увеличение chain позволяет выполнить больше последовательных `mul + rescale`, но само по себе не улучшает точность одной и той же операции.

## Целевая точность `1e-9`

Для `target_error = 1e-9` нужно около 30 бит финальной точности:

```text
ceil(-log2(1e-9)) = 30
```

CKKS теряет биты на encoding, encryption noise, multiplication, rescale, key switching и линейные трансформации. Для двух последовательных `mul_relin_rescale` текущая калибровка planner-а использует около 14 бит запаса:

```text
work_bits = 30 + 14 ~= 44
```

Практические режимы planner-а:

```text
speed:        scale_log2 = 45, рабочие модули 45 бит
conservative: scale_log2 = 50, рабочие модули 50 бит
```

Проверка `bench_parameter_planner` на реальном SEAL-прогоне показывает, что выбранные профили проходят цель `1e-9` для leveled-операций глубины 1, 2 и 3, а также для бюджетов `mul + linear_transform` и `mul + EvalMod P3`.

## Подбор параметров

`plan_ckks_parameters` принимает:

```text
target_error
multiplicative_depth
slots
ops_profile или CkksOperationBudget
security_bits
```

Результат содержит:

```text
CkksProfile
selected_work_bits
selected_scale_log2
estimated_precision_bits
estimated_abs_error_bound
passes_target_error
SecurityReport
```

Для простых сценариев можно использовать `ops_profile`. Для реальных pipeline лучше передавать `CkksOperationBudget`, где отдельно указаны:

- ciphertext multiplication;
- plaintext multiplication;
- plaintext-rescale;
- explicit rescale;
- mod-switch;
- rotations;
- linear transforms;
- EvalMod P3;
- refresh.

## Слоты

В CKKS Microsoft SEAL число слотов равно:

```text
slots = N / 2
```

Примеры:

```text
4096 / 2 = 2048 slots
8192 / 2 = 4096 slots
16384 / 2 = 8192 slots
32768 / 2 = 16384 slots
```

В bootstrap layout логическое число bootstrap slots не должно сужать физический `CkksProfile::slots`: diagonal masks и rotation-based transforms должны кодироваться на полном CKKS slot count.

## Проверка безопасности

Профили проверяются через `SealAdapter::create` и `SecurityReport` относительно лимитов Microsoft SEAL для `tc128`, `tc192` и `tc256`.

Подбор параметров балансирует:

- доступную мультипликативную глубину;
- численную точность после rescale;
- суммарный размер коэффициентного модуля;
- `poly_modulus_degree`, необходимый для выбранного уровня безопасности.

## Bootstrapping-точность

Обычная модель `work_bits = result_bits + calibrated_loss_bits` описывает leveled CKKS arithmetic, но не доказывает точность bootstrapping. Для refresh-круга используется recurrence:

```text
e_{k+1} <= A * e_k + b
```

Для итоговой точности `1e-9` на нескольких refresh-кругах вклад одного круга должен быть существенно ниже `1e-9`. Поэтому bootstrap-профили с `N = 32768`, длинной chain и масштабом `2^50`/`2^59` используются как диагностические профили для проверки DFT, rotation noise, ModUp и EvalMod по отдельности.
