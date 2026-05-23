# Архитектура

WonderFullyHE реализует библиотеку защищённых вычислений поверх Microsoft SEAL.

```text
Демонстрационные приложения и эксперименты
        |
        v
Публичный API WonderFullyHE
        |
        +-- SealAdapter: CKKS-контекст, ключи, encode/encrypt/evaluate/decrypt
        +-- accuracy: отчёты по максимальной и средней ошибке
        +-- abft: checksum-проверки корректности вычислений
        +-- LinearTransform: линейные преобразования через ротации
        +-- PolynomialEvaluator: вычисление полиномов от ciphertext
        +-- Bootstrapper: диагностика глубины и отчёт по bootstrapping-конвейеру
        +-- profile_report: воспроизводимые таблицы CKKS-параметров
        |
        v
Microsoft SEAL как CKKS-движок
```

## Слои публичного API

### SealAdapter

`SealAdapter` отделяет публичный API библиотеки от Microsoft SEAL. Прямые SEAL-типы скрыты за `Plain` и `Cipher`, а наружу вынесены CKKS-операции, используемые библиотекой:

- создание CKKS-контекста из `CkksProfile`;
- генерация ключей;
- encode/decode;
- encrypt/decrypt;
- сложение и вычитание;
- сложение и вычитание ciphertext с plaintext;
- умножение ciphertext на plaintext с rescale;
- умножение с relinearize и rescale;
- ротация вектора;
- генерация только нужных Galois-ключей для заданных ротаций;
- выравнивание уровня и масштаба ciphertext;
- размеры сериализованных объектов;
- диагностика ciphertext: `scale`, `chain_index`, `coeff_modulus_size`.

Microsoft SEAL остаётся внутренней зависимостью. Текущая C++-реализация размещает публичные типы и функции в пространстве имён `m2424`.

### accuracy

Модуль `accuracy` задаёт общие критерии численной точности:

```text
max_abs_error
mean_abs_error
compare(expected, actual, tolerance)
```

Демонстрационные приложения и тесты используют этот модуль как единый критерий точности.

### abft

Модуль `abft` реализует checksum-проверки корректности. Текущие проверки покрывают:

- добавленный checksum для сложения;
- добавленный checksum для вычитания;
- эталонный checksum для покомпонентного умножения;
- сохранение суммы при ротации.

ABFT-слой используется для проверки численной согласованности защищённых вычислений.

### LinearTransform

`LinearTransform` реализует преобразования вида:

```text
sum_i a_i * rotate(ct, k_i)
```

Этот блок нужен для этапов `CoeffToSlot` и `SlotToCoeff`, где вычисление строится из ротаций CKKS-слотов, умножений на plaintext-коэффициенты и сложений.

Отдельная функция `sum_slots` считает сумму заданного числа слотов и помещает результат в первый слот ciphertext. Для неё используется последовательность ротаций по степеням двойки.

### PolynomialEvaluator

`PolynomialEvaluator` вычисляет полином от ciphertext по заданным степеням и коэффициентам. Блок нужен для этапа `EvalMod`, где модульная редукция приближается полиномом.

### Bootstrapper

`Bootstrapper` является точкой входа для bootstrapping-модуля. Текущая реализация предоставляет:

- диагностику мультипликативной глубины;
- отчёт по CKKS bootstrapping-конвейеру;
- явные этапы: `ModRaise`, `CoeffToSlot`, `EvalMod`, `SlotToCoeff`.
- параметры ciphertext на границе глубины: `scale`, `chain_index`, `coeff_modulus_size`, `serialized_bytes`;
- критерии из модели: `Dec(c') ≈ Dec(c)` и `level(c') > level(c)`.

Модуль отделяет анализ глубины и состояние bootstrapping-конвейера от низкоуровневого адаптера SEAL.

### profile_report

`profile_report` формирует воспроизводимый отчёт по CKKS-параметрам:

```bash
./build/demo_profile_report
```

Вывод содержит `N`, число слотов, цепочку модулей, суммарный размер коэффициентного модуля, масштаб и оценку мультипликативной глубины.

## Демонстрационные приложения

- `demo_secure_stats` показывает защищённую облачную агрегацию: сумму и среднее над зашифрованными данными.
- `demo_abft` проверяет корректность гомоморфных операций через ABFT-инварианты.
- `demo_noise_growth` показывает расход глубины и остановку вычисления без bootstrapping.
- `demo_bootstrap_pipeline` печатает отчёт bootstrapping-модуля.
- `bench_ckks` измеряет время операций, численную ошибку и размеры сериализованных объектов.
- `bench_bootstrap_parts` измеряет `mul_plain_rescale`, `linear_transform`, `sum_slots` и `polynomial_eval`.
- `demo_profile_report` печатает таблицу CKKS-параметров.

## Следующие шаги реализации

Следующий этап реализации — развитие вычислительных блоков bootstrapping:

1. подстановка рассчитанных матриц `CoeffToSlot` и `SlotToCoeff` в `LinearTransform`;
2. подстановка коэффициентов полинома `EvalMod` в `PolynomialEvaluator`;
3. сборка `Bootstrapper::refresh(cipher)` и проверка критериев `Dec(c') ≈ Dec(c)`, `level(c') > level(c)` на end-to-end сценарии.
