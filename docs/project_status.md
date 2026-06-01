# Статус проекта

WonderFullyHE — учебно-исследовательский прототип библиотеки для защищённых вычислений на CKKS поверх Microsoft SEAL.

## Реализовано

- Адаптер `SealAdapter`, скрывающий низкоуровневые типы Microsoft SEAL.
- Готовые CKKS-профили в `m2424::profiles`: `fast_demo_ckks`, `basic_ckks`, `balanced_ckks`, `depth_ckks`, `high_precision_ckks`.
- Первая версия `CkksParameterPlanner`: подбор `scale`, рабочих модулей, длины chain и `N` из `target_error`, depth, slots и security.
- `OperationBudgetBuilder` и автоматический tracking бюджета в `CheckedEvaluator` для базовых проверяемых pipeline.
- `plan_bootstrap_refresh`: первый gate, который связывает текущий `CipherInfo`, бюджет следующего блока и решение `compute_fits_without_refresh` / `refresh_required` / `refresh_plan_blocked`.
- `Bootstrapper::plan_refresh_for_budget`: применение refresh gate к реальному ciphertext перед запуском experimental refresh.
- Операции `encode`, `encrypt`, `decrypt`, `decode`.
- Гомоморфные операции `add`, `sub`, `mul_relin_rescale`, `rotate`.
- Plaintext-операции `add_plain`, `sub_plain`, `mul_plain`, `mul_plain_rescale`.
- Кодирование plaintext на уровне ciphertext через `encode_like`, `encode_scalar_like` и `encode_scalar_at_scale_like`.
- Генерация только нужных Galois-ключей для заданных rotation steps.
- Сериализация и загрузка публичного ключа, секретного ключа, Relin/Galois-ключей и ciphertext.
- Разделённый roundtrip-сценарий: шифрование данных, вычисление без secret key, расшифрование результата.
- Выравнивание ciphertext по уровню и масштабу через `mod_switch_to` и `match_level_and_scale`.
- Метрики точности `max_abs_error`, `mean_abs_error`, `tolerance`.
- ABFT-проверки для `add`, `sub`, `mul` и `rotate`.
- Rotation-based `LinearTransform`.
- `sum_slots` для суммы слотов в первом CKKS-слоте.
- `PolynomialEvaluator` для вычисления полиномов от ciphertext.
- Демонстрация защищённой статистической агрегации.
- Демонстрация оптимизации Galois-ключей: полный набор против ограниченного набора ротаций.
- Анализ расхода вычислительной глубины до остановки на `x^32`.
- Низкоуровневый CKKS `ModRaise` для расширения ciphertext к первой RNS-базе.
- Experimental bootstrapping-конвейер `ModRaise -> CoeffToSlot -> eval_mod_normalization -> EvalMod -> SlotToCoeff` с отчётом по `scale`, `chain_index` и размеру ciphertext.
- Public architecture layer `BootstrapPipelinePlan`, который отделяет research backend `DenseDiagonal` от целевого scalable backend `FftLike` и фиксирует stage contracts.
- `plan_bootstrap_layout` строит физический `CkksProfile` из рассчитанного layout. Важно: логические bootstrap slots не сужают `CkksProfile::slots`; профиль оставляет полный CKKS slot count, чтобы diagonal masks и rotation transforms кодировались физически корректно.
- Scaling gate перед full refresh validation. Найденный leak был в normalization scalar: после `ModRaise -> CoeffToSlot` tiny scalar нельзя применять одним plaintext multiply. Scaling layer теперь поддерживает decomposition на несколько `mul_plain_rescale` шагов; gate имеет PASS-режимы, а следующий blocker находится в full EvalMod/denormalization path.
- One-case trace отделяет `scalar_pass` от `evalmod_ready`: `NoBootstrapPeriod` является diagnostic baseline, а full validation не запускается без реального period-mode, который одновременно scalar-correct и попадает в интервал EvalMod.
- Period-model gate показал следующий blocker: реальные period-mode могут попасть в EvalMod interval, но пока не дают одновременно достаточный уровень и безопасный scale перед P3.
- Experimental API `Bootstrapper::refresh(cipher, slots, tolerance)` и стабильный helper `Bootstrapper::refresh_rotation_steps(slots)`.
- Scalable refresh API `Bootstrapper::refresh_slots_to_coeffs_first(...)` и `Bootstrapper::scalable_refresh_rotation_steps(slots)`: публичный путь использует `SlotsToCoeffsFirst + FftLike + EvalMod P3` и покрыт отдельным CTest.
- `output_scale_repair` после `EvalMod P3`: value-preserving plaintext-rescale возвращает output scale к рабочему диапазону около `2^40`, что позволяет запускать несколько refresh-кругов подряд.
- `bench_bootstrap_multi_cycle`: текущий baseline после сжатия FFT-like `CoeffToSlot` — `slots=16, cycles=3` проходит с `max_abs_error < 2e-5`, минимум `8` continuation levels. Проверочный 57-битный layout `60 + 13*57 + 60` тоже проходит `slots=16, cycles=3`, но ошибка остаётся около `1e-5`.
- Нормализация входа `EvalMod` по амплитуде после `CoeffToSlot`.
- Historical end-to-end demo для refresh; не входит в default CTest, пока scaling gate не проходит.
- Benchmark публичного refresh-пути.
- Benchmark времени операций, ошибок и размеров ciphertext/ключей.
- Benchmark `bench_chain_accuracy` для контролируемой проверки влияния длины chain, `scale_log2` и рабочей битности на точность.
- Benchmark `bench_parameter_planner` для проверки выбранных planner-ом профилей на реальном SEAL-прогоне.
- Benchmark строительных блоков bootstrapping: plaintext multiplication, linear transform, slot summation, polynomial evaluation.
- Benchmark параллельной обработки независимых ciphertext с разделением setup/runtime.
- Calibrated bootstrap precision model: recurrence `e_{k+1} <= A e_k + b`, rotation scale model и diagnostic profile `precision_boot_ultra_ckks_59` для проверки `scale = 2^59` на Microsoft SEAL.
- Manual diagnostics для slots=4 DFT precision: `precision_boot_ultra_ckks_59` снижает rotation/key-switch noise ниже `1e-10`, а current `FftLike`/`DenseDiagonal` DFT roundtrip проходит `2e-10` gate в диагностическом прогоне. `SmallSlots4Butterfly` семантически совпадает с reference STC/CTS, но не выбран для refresh, потому что его больший rescale depth даёт худшую encrypted roundtrip error.
- Security report по CKKS-профилям относительно лимитов Microsoft SEAL `tc128`, `tc192` и `tc256`.
- CMake alias target `m2424::m2424` для подключения библиотеки через `add_subdirectory`.
- CMake-сборка, smoke-тесты, GitLab CI и GitHub Actions CI.

## В работе

- Расширение калибровочной модели `ops_profile -> calibrated_loss_bits`; для target `1e-9` текущий быстрый ориентир — 45-битные рабочие модули при `scale_log2 = 45`.
- `BootstrapPlanner` как обязательный gate перед refresh: проверка уровней, scale, period window, EvalMod interval и error budget до ciphertext-прогона.
- Factorized FFT-like backend для `CoeffToSlot/SlotToCoeff`; текущий dense diagonal backend остаётся reference для малых `slots`.
- Довести `SlotsToCoeffsFirst` от guarded refresh до полноценного level-refresh: сейчас путь сохраняет значение, проходит `EvalMod P3` на `boot_deep_ckks` и поддерживает повторные циклы; для `slots=4` следующий gate — перенос successful DFT precision result в refresh path без переключения на `SmallSlots4Butterfly`, пока этот backend не станет lowest-error candidate.
- Следующий blocker: сырой `post_refresh_mod_raise` поднимает chain, но нарушает tolerance, поэтому нужен корректный output scaling/denormalization перед post-raise, а не простой structural raise.
- Оптимизация rotation keys, BSGS/hoisting и кеширования plaintext-диагоналей.

## Ограничения

- ABFT-модуль контролирует вычислительную согласованность, но не заменяет криптографическую аутентификацию результата.
- Для реального внешнего применения требуется отдельное управление ключами, формат обмена данными и анализ side-channel рисков.
- Full bootstrapping пока не является стабильным API: текущий код имеет stage-by-stage gates и oracle-harness, но не гарантирует бесконечное число операций с ошибкой `1e-9`.
- `precision_boot_ultra_ckks_59` является diagnostic scale/noise profile, а не production bootstrap profile. Финальный профиль должен выводиться из segmented bootstrap layout после прохождения DFT, ModRaise, EvalMod и multi-cycle gates.
- Повышать степень `EvalMod` выше `P3` не является текущим способом улучшить точность: при `|u| <= 2^-10` ошибка `P3` уже меньше целевого порога, а bottleneck находится в scale/noise/transforms.
